// Copyright (c) 2025 Outbreak Preparedness and Response Group at BNITM
// This file is part of MOPX and is licensed under the MIT License.
// See the LICENSE file in the root of this repository for full license details.

nextflow.enable.dsl = 2

process concat_fastq_files {
    label "common"
    cpus 1
    input:
        tuple val(meta), path('input_dir')
    output:
        tuple val(meta), path('concat.fastq')
    """
    #!/usr/bin/env python3

    from pathlib import Path
    import gzip

    fcount = 0

    with open('concat.fastq', 'w') as f_out:
        for fname in Path('input_dir').iterdir():
            if fname.name.endswith('.fq') or fname.name.endswith('.fastq'):
                fcount += 1
                with open(fname) as f:
                    f_out.write(f.read())
            elif fname.name.endswith('.fq.gz') or fname.name.endswith('.fastq.gz'):
                fcount += 1
                with gzip.open(fname, 'rt') as f:
                    f_out.write(f.read())

    if fcount == 0:
        raise RuntimeError('No fastq files found!')
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
    set -e

    if canu \\
        -nanopore-raw ${reads} \\
        -p asm \\
        -d out \\
        genomeSize=${params.canu_genome_size} \\
        minReadLength=${params.canu_min_read_length} \\
        minOverlapLength=${params.canu_min_overlap_length} \\
        stopOnLowCoverage=${params.canu_stop_on_low_coverage} \\
        maxThreads=${task.cpus} \\
        minInputCoverage=${params.canu_min_input_coverage} \\
        maxMemory=${task.memory.toGiga()}g; then

        mv out/asm.contigs.fasta canu_contigs.fasta

    else
        echo "Canu failed; creating empty contigs file."
        touch canu_contigs.fasta
    fi
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


process assign_terminal_hairpin_region {
    label "common"
    input:
        tuple val(meta), path("polished_draft.fasta"), path("hairpin.fasta")
    output:
        tuple val(meta), path("hairpin_region_to_refine.fasta"), path("hairpin_matching.json")
    """
    match_terminal_hairpin.py \\
        --hairpin hairpin.fasta \\
        --genome polished_draft.fasta \\
        --margin 300 \\
        --flanking_region_diff_tolerance 50 \\
        --out .
    """
}


process minimap_nosecondary {
    label "common"
    cpus 2
    input:
        tuple val(meta), path(ref), path(reads)
    output:
        tuple val(meta), path(ref), path('sorted.bam'), path('sorted.bam.bai')
    """
    minimap2 -Y -ax map-ont $ref $reads --secondary=no \
    | samtools view -b -o mapped.bam
    samtools sort -o sorted.bam mapped.bam
    samtools index sorted.bam
    """
}


process medaka_variants {
    label "medaka"
    cpus 2
    input:
        tuple val(meta),
            path("draft.fasta"),
            path("mapped.bam"),
            path("mapped.bam.bai"),
            path(model)
    output:
        tuple val(meta), path("annotated.vcf")
    """
    if [[ -s "draft.fasta" ]]
    then
        medaka inference mapped.bam consensus_probs.hdf --model "$model" --threads $task.cpus

        medaka vcf consensus_probs.hdf draft.fasta variants.vcf --gvcf

        bcftools sort variants.vcf -o sorted.vcf

        medaka tools annotate sorted.vcf draft.fasta mapped.bam annotated.vcf
    else
        touch annotated.vcf
    fi
    """
}


process make_hairpin_consensus {
    label "common"
    input:
        tuple val(meta),
            path("hairpin_matching.json"),
            path("variants.vcf")
    output:
        tuple val(meta),
            path("hairpin_consensus.fasta")
    """
    if [[ -s "variants.vcf" ]]
    then
        build_refined_terminal_hairpin.py \\
            --hairpin-matches hairpin_matching.json \\
            --gvcf variants.vcf \\
            --min-depth ${params.min_coverage} \\
            --min-qual ${params.medaka_variant_minqual} \\
            --out hairpin_consensus.fasta
    fi
    touch hairpin_consensus.fasta
    """
}


process integrate_fixed_hairpin {
    label "common"
    input:
        tuple val(meta),
            path("corrected_draft.fasta"),
            path("basecounts.tsv"),
            path("hairpin_matching.json"),
            path("hairpin_consensus.fasta")
    output:
        tuple val(meta), path("fixed_hairpin_draft.fasta"), path("fixed_hairpin_basecounts.tsv") 
    """
    if [[ -s hairpin_consensus.fasta ]]
    then
        integrate_fixed_hairpin.py \\
            --draft corrected_draft.fasta \\
            --basecounts basecounts.tsv \\
            --hairpin-consensus hairpin_consensus.fasta \\
            --hairpin-matches hairpin_matching.json \\
            --out .
    else
        cp corrected_draft.fasta fixed_hairpin_draft.fasta
        cp basecounts.tsv fixed_hairpin_basecounts.tsv
    fi
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
        tuple val(meta), path("draft.fasta"), path("sorted.bam"), path("sorted.bam.bai")
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


process gb2fasta {
    label "common"
    cpus 1
    input:
        tuple val(meta), path("reference.gb")
    output:
        tuple val(meta), path("reference.fasta")
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    gb = SeqIO.parse('reference.gb', 'genbank')
    SeqIO.write(gb, 'reference.fasta', 'fasta')
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

    snpEff build \\
        -configOption "\${ref_id}.genome=\${ref_id}" \\
        -dataDir \$PWD/data \\
        -genbank -v \\
        \${ref_id}

    snpEff \\
        -configOption "\${ref_id}.genome=\${ref_id}" \\
        -dataDir \$PWD/data \\
        -v -no-downstream -no-upstream -no-utr \\
        \${ref_id} \\
        variants.vcf.gz > annotated.vcf

    bgzip annotated.vcf
    tabix -p vcf annotated.vcf.gz
    """
}


process report {
    label "common"
    cpus 1
    input:
        tuple val(samplename), path("final_draft.fasta"), path("variants.vcf.gz")
    output:
        tuple val(samplename), path("report.html")
    """
    report.py \\
        --samplename ${samplename} \\
        --mopx-version ${workflow.manifest.version} \\
        --draft final_draft.fasta \\
        --vcf variants.vcf.gz \\
        --out report.html
    """
}


process output {
    label "common"
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
