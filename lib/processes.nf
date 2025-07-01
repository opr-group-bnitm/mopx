nextflow.enable.dsl = 2


process concat_fastq_files {
    label "common"
    cpus 1
    input:
        tuple val(meta), path('input_dir')
    output:
        tuple val(meta), path('concat.fastq')
    """
    set -euo pipefail

    touch concat.fastq

    # Check if there are any .fastq or .fastq.gz files
    n_fastq=\$(find input_dir/ -type f \\( -name '*.fastq' -o -name '*.fq' \\) | wc -l)
    n_fastq_gz=\$(find input_dir/ -type f \\( -name '*.fastq.gz' -o -name '*.fq.gz' \\) | wc -l)
    total=\$(( n_fastq + n_fastq_gz ))

    if [[ "\${total}" -eq 0 ]]
    then
        echo "ERROR: No .fastq or .fastq.gz files found" >&2
        exit 1
    fi

    # Concatenate uncompressed .fastq files
    find input_dir/ -type f \\( -name '*.fastq' -o -name '*.fq' \\) | while read -r fq
    do
        cat "\${fq}" >> concat.fastq
    done

    # Concatenate compressed .fastq.gz files
    find input_dir/ -type f \\( -name '*.fastq.gz' -o -name '*.fq.gz' \\) | while read -r fqz
    do
        gzip -dc "\${fqz}" >> concat.fastq
    done
    """
}


process filter_reads {
    label "common"
    cpus 4
    input:
        tuple val(meta), path('reads.fastq'), path('ref.fasta')
    output:
        tuple val(meta), path('filtered.fastq')
    """
    minimap2 \\
        -ax map-ont \\
        --secondary=no \\
        ref.fasta \\
        -t ${task.cpus} \\
        reads.fastq \\
    | samtools fastq \\
        -F 4 > filtered.fastq
    """
}


process canu {
    label "canu"
    cpus 16
    memory "30 GB"
    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('canu_contigs.fasta')
    """
    canu \\
        -nanopore-raw ${reads} \\
        -p asm \\
        -d out \\
        genomeSize=${params.canu_genome_size} \\
        minReadLength=${params.canu_min_read_length} \\
        minOverlapLength=${params.canu_min_overlap_length} \\
        stopOnLowCoverage=${params.canu_stop_on_low_coverage} \\
        maxThreads=${task.cpus} \\
        minInputCoverage=${params.canu_min_input_coverage} \\
        maxMemory=${task.memory.toGiga()}g

    mv out/asm.contigs.fasta canu_contigs.fasta
    """
}


process pop_canu_bubbles {
    label "common"
    cpus 1
    input:
        tuple val(meta), path('canu_contigs.fasta')
    output:
        tuple val(meta), path('canu_contigs.nobbubles.fasta')
    """
    #!/usr/bin/env python

    from Bio import SeqIO

    with open('canu_contigs.nobbubles.fasta', 'w') as out_handle:
        for record in SeqIO.parse('canu_contigs.fasta', 'fasta'):
            if 'suggestBubble=yes' not in record.description:
                SeqIO.write(record, out_handle, 'fasta')
    """
}


process flye {
    label "flye"
    cpus 16
    memory "30 GB"
    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path('flye_contigs.fasta')
    """
    flye \\
        --nano-raw ${reads} \\
        --out-dir out \\
        --genome-size ${params.flye_genome_size} \\
        --threads ${task.cpus}

    mv out/assembly.fasta flye_contigs.fasta
    """
}


process merge_and_cluster {
    label "common"
    input:
        tuple val(meta), path('contigs_*.fasta')
    output:
        tuple val(meta), path('clustered.fasta')
    """
    cat contigs_*.fasta > merged.fasta
    cd-hit-est \\
        -n ${params.contigs_cluster_cdhit_wordlen} \\
        -c ${params.contigs_cluster_cdhit_min_similarity} \\
        -i merged.fasta -o clustered.fasta
    """
}


process split_ref {
    label "common"
    cpus 1
    input:
        tuple val(meta), val(itr_len), path('ref.fasta')
    output:
        tuple val(meta), path('itr.fasta'), emit: itr
        tuple val(meta), path('body.fasta'), emit: core
    """
    #!/usr/bin/env python

    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    ref = SeqIO.read('ref.fasta', 'fasta')

    itr_record = SeqRecord(
        id='ITR',
        description=f"ref={ref.id}",
        seq=ref.seq[:${itr_len}],
    )
    body_record = SeqRecord(
        id='core_region',
        description=f"ref={ref.id}",
        seq=ref.seq[${itr_len}:-${itr_len}],
    )

    SeqIO.write(itr_record, 'itr.fasta', 'fasta')
    SeqIO.write(body_record, 'body.fasta', 'fasta')
    """
}


process minimap {
    label "common"
    cpus 2
    input:
        tuple val(meta), path(ref), path(reads)
    output:
        tuple val(meta), path(ref), path('sorted.bam'), path('sorted.bam.bai')
    """
    minimap2 -Y -ax map-ont $ref $reads -p $params.readmap_min_secondary2primary_ratio -N $params.readmap_max_secondary \
    | samtools view -b -o mapped.bam
    samtools sort -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}


process minimap_eqx {
    label "common"
    cpus 2
    input:
        tuple val(meta), path(ref), path(reads)
    output:
        tuple val(meta), path(ref), path('sorted.bam'), path('sorted.bam.bai')
    """
    minimap2 -eqx -Y -ax map-ont $ref $reads -p $params.readmap_min_secondary2primary_ratio -N $params.readmap_max_secondary \
    | samtools view -b -o mapped.bam
    samtools sort -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}


process sniffles {
    label "structural_variants"
    cpus 2
    memory "10 GB"
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("mapped_to_ref.bam"),
            path("mapped_to_ref.bam.bai")
    output:
        tuple val(meta), path("ref.fasta"), path("sv.filtered.vcf")
    """
    if [[ \$(samtools view mapped_to_ref.bam | wc -l) -ne 0 ]]
    then
        set +e
        sniffles \\
            -t ${task.cpus} \\
            --input mapped_to_ref.bam \\
            --reference ref.fasta \\
            --vcf sv.vcf \\
            --minsupport ${params.sniffles_predraft_min_support} \\
            --minsvlen ${params.sniffles_predraft_min_sv_len}
        set -e

        bcftools view -i 'INFO/AF>=${params.sniffles_predraft_min_variant_allele_fraction}' sv.vcf -Ov -o sv.filtered.vcf
    else
        # create an empty vcf file
        refid=\$(head -n 1 ref.fasta | sed 's/^>//' | awk '{print \$1}')
        reflen=\$(grep -v '^>' ref.fasta | tr -d '\\n' | wc -c)

        # Create the minimal VCF
        {
            echo "##fileformat=VCFv4.2"
            echo "##contig=<ID=\$refid,length=\$reflen>"
            echo -e "#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO"
        } > sv.filtered.vcf
    fi
    """
}


process structural_variant_consensus {
    label "structural_variants"
    input:
        tuple val(meta),
            path("ref.fasta"),
            path("sv.filtered.vcf")
    output:
        tuple val(meta),
            path("sv_consensus.fasta")
    """
    if [[ \$(bcftools view -H sv.filtered.vcf | wc -l) -ne 0 ]]
    then
      bcftools sort sv.filtered.vcf -Oz -o structural_variants.vcf.gz
      tabix structural_variants.vcf.gz
      bcftools consensus \\
        --fasta-ref ref.fasta \\
        -o sv_consensus.fasta structural_variants.vcf.gz
    else
      cp ref.fasta sv_consensus.fasta
    fi
    """
}


process secondary_to_primary_alignments {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("draft.fasta"), path("in.bam"), path("in.bam.bai")
    output:
        tuple val(meta), path("draft.fasta"), path("out.bam"), path("out.bam.bai")
    """
    transform_bam.py --input in.bam --output unsorted.bam
    samtools sort unsorted.bam -o out.bam
    samtools index out.bam
    """
}


process medaka_polish {
    label "medaka"
    cpus 2
    input:
        tuple val(meta), path("draft.fasta"), path("mapped.bam"), path("mapped.bam.bai"), path(model)
    output:
        tuple val(meta), path("polished.fasta")
    """
    medaka inference mapped.bam consensus_probs.hdf \
        --model "$model" --threads $task.cpus \
        || { echo "Medaka inference failed."; exit 1; }

    medaka sequence consensus_probs.hdf draft.fasta polished.fasta \
        --threads $task.cpus \
        || { echo "Consensus stitching failed."; exit 1; }
    """
}


process merge_predraft {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("segment_*.fasta")
    output:
        tuple val(meta), path("merged.fasta")
    """
    #!/usr/bin/env python3

    import glob
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    sequences = [SeqIO.read(fn, 'fasta') for fn in glob.glob("segment_*.fasta")]
    itr_seqs = [s for s in sequences if s.id == "ITR"]
    core_seqs = [s for s in sequences if s.id == "core_region"]

    assert len(itr_seqs) == 1
    assert len(core_seqs) == 1

    itr = itr_seqs[0]
    core = core_seqs[0]
    full_seq = itr.seq + core.seq + itr.seq.reverse_complement()

    outseq = SeqRecord(
        id="PreDraft",
        description=itr.description,
        seq=full_seq
    )

    SeqIO.write(outseq, 'merged.fasta', 'fasta')
    """
}


process place_contigs {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("ref.fasta"), path("contigs.fasta")
    output:
        tuple val(meta), path("draft.fasta"), path("draft.log")
    """
    place_contigs.py \
        --ref ref.fasta \
        --contigs contigs.fasta \
        --outfile draft.fasta \
        --log draft.log \
        --draftid $meta.draft_name
    """
}


process align_with_terminal_repeats {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("draft.fasta"), path(reads), path("tr_range.txt")
    output:
        tuple val(meta), path("draft.fasta"), path("mapped.bam"), path("mapped.bam.bai")
    """
    seq_id=\$(grep '^>' "draft.fasta" | head -n 1 | cut -d ' ' -f 1 | sed 's/^>//')

    non_repeat_start=\$(head -n 1 tr_range.txt)
    non_repeat_end=\$(tail -n 1 tr_range.txt)

    samtools faidx draft.fasta
    samtools faidx draft.fasta "\$seq_id:1-\${non_repeat_start}" > prefix.fasta
    samtools faidx draft.fasta "\$seq_id:\${non_repeat_end}-" > suffix.fasta

    for ref in draft prefix suffix
    do
        minimap2 -Y -ax map-ont \${ref}.fasta $reads --secondary=no \
        | samtools view -b -o unsorted_\${ref}.bam
        samtools sort -o \${ref}.bam unsorted_\${ref}.bam
        samtools index \${ref}.bam
    done

    merge_draft_alignments.py \
        --max-score-diff $params.terminal_repeat_max_map_diff \
        --all draft.bam \
        --left prefix.bam \
        --right suffix.bam \
        --out merged_unsorted.bam

    samtools sort merged_unsorted.bam -o mapped.bam
    samtools index mapped.bam
    """
}


process assign_terminal_repeat_regions {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("ref.fasta"), path("draft.fasta"), val(trl_end), val(trr_start)
    output:
        tuple val(meta), path("draft.fasta"), path("alignment.fasta"), path("tr_range.txt")
    """
    cat ref.fasta draft.fasta > concat.fasta
    mafft --auto concat.fasta > alignment.fasta

    transfer_terminal_repeat_range.py \
        --alignment alignment.fasta \
        --trl-end $trl_end \
        --trr-start $trr_start \
        --out tr_range.txt
    """
}


process correct_and_mask_consensus {
    label "common"
    cpus 1
    input:
        tuple val(meta), path('draft.fasta'), path('sorted.bam'), path('sorted.bam.bai')
    output:
        tuple val(meta), path("consensus_masked.fasta"), path("nucleotide_counts.tsv")
    """
    mask_consensus.py \
        --bam sorted.bam \
        --reference draft.fasta \
        --outdir . \
        --min-coverage $params.min_coverage
    """
}


process evaluate_contigs {
    label "common"
    cpus 1
    input:
        tuple val(meta), path('mapped_contigs.bam'), path('mapped_contigs.bai')
    output:
        tuple val(meta), path('contig_mapping_stats.tsv')
    """
    mapped_contigs_evaluation.py --bam mapped_contigs.bam --out contig_mapping_stats.tsv
    """
}


process mafft_align {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("seq_*.fasta")
    output:
        tuple val(meta), path("alignment.fasta")
    """
    cat seq_*.fasta > concat.fasta
    mafft --auto concat.fasta > alignment.fasta
    """
}


process shift_gaps {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("alignment_in.fasta")
    output:
        tuple val(meta), path("alignment_out.fasta")
    """
    shift_gaps.py --align alignment_in.fasta --out alignment_out.fasta
    """
}


process variants_from_alignment {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("alignment.fasta"), path("basecounts.tsv")
    output:
        tuple val(meta), path("variants.vcf.gz"), path("variants.vcf.gz.tbi")
    """
    variants_from_alignment.py --align alignment.fasta --basecounts basecounts.tsv --out variants.vcf    
    bgzip variants.vcf
    tabix -p vcf variants.vcf.gz
    """
}


process variants_report {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("annotations.gb"), path("variants.vcf.gz"), path("variants.vcf.gz.tbi")
    output:
        tuple val(meta), path("report.txt")
    """
    variants_report.py --variants variants.vcf.gz --gb annotations.gb > report.txt
    """
}


process annotate_cds_changes {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("alignment.fasta"), path("annotations.gb")
    output:
        tuple val(meta), path("alignments.txt"), path("transfers.tsv")
    """
    annotate.py --alignment alignment.fasta --gb annotations.gb --out out
    mv out/transfers.tsv transfers.tsv
    mv out/alignments.txt alignments.txt
    """
}


process mutational_spectrum {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("ref.fasta"), path("variants.vcf.gz"), path("variants.vcf.gz.tbi")
    output:
        tuple val(meta), path("pairs.tsv"), path("triples.tsv")
    """
    mutational_spectrum.py --vcf variants.vcf.gz --ref ref.fasta  --out .
    """
}


process annotate_aa_changes {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("reference.gbk"), path("variants.vcf.gz"), path("variants.vcf.gz.tbi")
    output:
        tuple val(meta), path("annotated.vcf.gz"), path("annotated.vcf.gz.tbi"), path("snpEff_genes.txt"), path("snpEff_summary.html")
    """
    ref_id=\$(awk '/^VERSION/ {print \$2}' reference.gbk)

    if [[ -z "\$ref_id" ]]
    then
        echo "Error: Could not get reference ID from genbank file." >&2 
        exit 1
    fi

    mkdir -p data/\${ref_id}
    cp reference.gbk data/\${ref_id}/genes.gbk

    if [[ ! -f data/\${ref_id}/genes.gbk ]]
    then
        echo "Error: genes.gbk file not found." >&2
        exit 1
    fi

    snpEff build \
        -configOption "\${ref_id}.genome=\${ref_id}" \
        -dataDir \$PWD/data \
        -genbank -v \
        \${ref_id}

    snpEff \
        -configOption "\${ref_id}.genome=\${ref_id}" \
        -dataDir \$PWD/data \
        -v -no-downstream -no-upstream -no-utr \
        \${ref_id} \
        variants.vcf.gz > annotated.vcf

    bgzip annotated.vcf
    tabix -p vcf annotated.vcf.gz
    """
}


process output {
    // publish inputs to output directory
    cpus 1
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: {
            dirname ? (fname_out ? "$dirname/$fname_out" : "$dirname/$fname_in") : fname_in
        }
    )
    input:
        tuple path(fname_in), val(dirname), val(fname_out)
    output:
        path(fname_in)
    """
    """
}
