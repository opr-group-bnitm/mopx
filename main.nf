nextflow.enable.dsl = 2

include {
    concat_fastq_files;
    filter_reads;
    canu;
    pop_canu_bubbles;
    flye;
    merge_and_cluster;
    split_ref;
    sniffles;
    structural_variant_consensus;
    merge_predraft;
    place_contigs;
    assign_terminal_repeat_regions;
    align_with_terminal_repeats;
    medaka_polish as medaka_polish_predraft;
    medaka_polish as medaka_polish_draft;
    minimap as minimap_to_separated_ref;
    minimap as minimap_to_separated_ref_sv;
    minimap as minimap_reads_to_draft;
    minimap_eqx as minimap_contigs_to_draft;
    minimap as minimap_reads_to_ref;
    minimap as minimap_contigs_to_ref;
    medaka_polish;
    correct_and_mask_consensus;
    evaluate_contigs;
    assign_terminal_hairpin_region;
    minimap_nosecondary;
    medaka_variants;
    make_hairpin_consensus;
    integrate_fixed_hairpin;
    // annotation
    gb2fasta;
    mafft_align;
    shift_gaps;
    variants_from_alignment;
    variants_report;
    annotate_cds_changes;
    annotate_aa_changes;
    // general
    output
} from './lib/processes.nf'


workflow de_novo_assemblies {
    take:
        reads
    main:
        canu_assemblies = reads
        | canu
        | pop_canu_bubbles
    
        flye_assemblies = reads
        | flye

        assemblies_merged = canu_assemblies
        | mix(flye_assemblies)
        | groupTuple(by: 0)
        | merge_and_cluster

    emit:
        contigs = assemblies_merged
}


workflow draft_genome {
    take:
        meta
    main:

        // part 1: build pre-draft
        ref_split = meta
        | map {meta -> [meta, params.reference_terminal_repeat_left_end, meta.reference]}
        | split_ref

        refs_in_pieces = ref_split.itr
        | mix(ref_split.core)

        predraft = refs_in_pieces
        | map {meta, ref -> [meta, ref, meta.reads]}
        | minimap_to_separated_ref
        | sniffles
        | structural_variant_consensus
        | map {meta, sv_consensus -> [meta, sv_consensus, meta.reads]}
        | minimap_to_separated_ref_sv
        | map {meta, draft, bam, bai -> [meta, draft, bam, bai, params.medaka_consensus_model]}
        | medaka_polish_predraft
        | map {meta, segment -> [meta.sample, segment]}
        | groupTuple(by: 0)
        | merge_predraft

        predraft_with_meta = meta
        | map {meta -> [meta.sample, meta]}
        | join(predraft, by: 0)
        | map {sampleid, meta, predraft -> [meta, predraft]}

        // part 2: place contigs
        raw_draft = predraft_with_meta
        | map {meta, predraft -> [meta, predraft, meta.contigs]}
        | place_contigs

        tr_annotated_draft = raw_draft
        | map { meta, draft, log -> [meta, meta.reference, draft, meta.trl_end, meta.trr_start]}
        | assign_terminal_repeat_regions

        polished_draft = tr_annotated_draft
        | map {meta, draft, alignment, tr_range -> [meta, draft, meta.reads, tr_range]}
        | align_with_terminal_repeats
        | map {meta, draft, bam, bai -> [meta, draft, bam, bai, params.medaka_consensus_model]}
        | medaka_polish

        alignments_reads_to_polished = polished_draft 
        | map {meta, draft -> [meta, draft, meta.reads]}
        | minimap_reads_to_draft

        corrected_draft = alignments_reads_to_polished
        | correct_and_mask_consensus
 
        // refine and integracte the hairpin

        hairpin_match = polished_draft
        | map { meta, polished -> [meta, polished, params.hairpin]}
        | assign_terminal_hairpin_region

        hairpin_alignments = hairpin_match
        | map { meta, hp_ref, hp_match -> [meta + [hairpin_match_info: hp_match], hp_ref, meta.reads]}
        | minimap_nosecondary

        hairpin_variants = hairpin_alignments
        | map { meta, ref, bam, bai -> [meta, ref, bam, bai, params.medaka_consensus_model] }
        | medaka_variants

        hairpin_consensus = hairpin_variants
        | map { meta, vcf -> [meta, meta.hairpin_match_info, vcf] }
        | make_hairpin_consensus

        draft_by_sample = corrected_draft
        | map { meta, draft, basecounts -> [meta.sample, draft, basecounts] }

        draft_with_corrected_hairpin = hairpin_consensus
        | map { meta, hp_cons -> [meta.sample, meta, hp_cons] }
        | join(draft_by_sample, by: 0)
        | map { sampleid, meta, hp_cons, draft, basecounts -> [meta, draft, basecounts, meta.hairpin_match_info, hp_cons] }
        | integrate_fixed_hairpin

        // Alignments for troubleshooting/further analysis
        // 1 - map contigs to polished draft
        contigs_mapped_to_polished = polished_draft
        | map {meta, draft -> [meta, draft, meta.contigs]}
        | minimap_contigs_to_draft

        contigs_eval = contigs_mapped_to_polished
        | map {meta, ref, bam, bai -> [meta, bam, bai]}
        | evaluate_contigs

        // 2 - map contigs to reference sequence (output directory: ref.fasta, mapped.bam, mapped.bam.bai)
        contigs_mapped_to_ref = meta
        | map {meta -> [meta, meta.reference, meta.contigs]}
        | minimap_contigs_to_ref

        // 3 - map reads to reference sequence (output directory: ref.fasta, mapped.bam, mapped.bam.bai)
        reads_mapped_to_ref = meta
        | map {meta -> [meta, meta.reference, meta.reads]}
        | minimap_reads_to_ref

        write_this = Channel.empty()
        | mix(
            raw_draft | map {meta, draft, log -> [draft, 'draft', 'draft.fasta']},
            raw_draft | map {meta, draft, log -> [log, 'logs', 'draft.log']},
            tr_annotated_draft | map {meta, draft, alignment, tr_range -> [alignment, 'repeat_transfer', 'alignment.fasta']},
            tr_annotated_draft | map {meta, draft, alignment, tr_range -> [alignment, 'repeat_transfer', 'range.txt']},
            polished_draft | map {meta, draft -> [draft, 'draft', 'polished.fasta']},
            draft_with_corrected_hairpin | map {meta, draft, basecounts -> [draft, 'draft', 'draft_final.fasta']},
            draft_with_corrected_hairpin | map {meta, draft, basecounts -> [basecounts, 'draft', 'basecounts_final.tsv']},
            alignments_reads_to_polished | map {meta, draft, bam, bai -> [bam, 'alignments', 'reads_vs_polished.bam']},
            alignments_reads_to_polished | map {meta, draft, bam, bai -> [bai, 'alignments', 'reads_vs_polished.bam.bai']},
            contigs_mapped_to_polished | map {meta, draft, bam, bai -> [bam, 'alignments', 'contigs_vs_polished.bam']},
            contigs_mapped_to_polished | map {meta, draft, bam, bai -> [bai, 'alignments', 'contigs_vs_polished.bam.bai']},
            contigs_mapped_to_ref | map {meta, draft, bam, bai -> [bam, 'alignments', 'contigs_vs_ref.bam']},
            contigs_mapped_to_ref | map {meta, draft, bam, bai -> [bai, 'alignments', 'contigs_vs_ref.bam.bai']},
            reads_mapped_to_ref | map {meta, draft, bam, bai -> [bam, 'alignments', 'reads_vs_ref.bam']},
            reads_mapped_to_ref | map {meta, draft, bam, bai -> [bai, 'alignments', 'reads_vs_ref.bam.bai']},
            contigs_eval | map {meta, table -> [table, 'stats', 'contigs_eval.tsv']}
        )
    emit:
        // emit everything that is used next
        draft = draft_with_corrected_hairpin
        write_this = write_this
}


workflow annotate_draft {
    take:
        draft_ref
    main:

        alignment = draft_ref
        | map {meta -> [meta, meta.genbank]}
        | gb2fasta
        | map {meta, ref_fasta -> [meta, [ref_fasta, meta.draft]]}
        | mafft_align
        | shift_gaps

        variants = alignment
        | map {meta, alignment -> [meta, alignment, meta.basecounts]}
        | variants_from_alignment

        // TODO: remove variants_report? and aa_report?
        // - replace with annotation based on ORF to AA alignments as done for
        //   the genome annotations? 
        variant_report = variants
        | map {meta, vcf, tbi -> [meta, meta.genbank, vcf, tbi]}
        | variants_report

        aa_report = alignment
        | map {meta, alignment -> [meta, alignment, meta.genbank]}
        | annotate_cds_changes

        vcf_variants = variants
        | map {meta, vcf, tbi -> [meta, meta.genbank, vcf, tbi]}
        | annotate_aa_changes

        write_this = Channel.empty()
        | mix(
            alignment | map {meta, align -> [align, 'annotation', 'alignment.fasta']},
            variant_report | map {meta, report -> [report, 'annotation', 'variant_report.txt']},
            vcf_variants | map {meta, vcf, tbi, genetable, var_summary -> [vcf, 'annotation', 'variants.vcf.gz']},
            vcf_variants | map {meta, vcf, tbi, genetable, var_summary -> [tbi, 'annotation', 'variants.vcf.gz.tbi']},
            vcf_variants | map {meta, vcf, tbi, genetable, var_summary -> [genetable, 'annotation', 'gene_variant_summary.txt']},
            vcf_variants | map {meta, vcf, tbi, genetable, var_summary -> [var_summary, 'annotation', 'variant_summary.html']},
            aa_report | map {meta, aa_align, transfers -> [aa_align, 'annotation', 'protein_alignments.txt']},
            aa_report | map {meta, aa_align, transfers -> [transfers, 'annotation', 'cds_transfer.tsv']}
        )
    emit:
        write_this = write_this
}


workflow {

    // TODO: implement multi-sample
    // fastqs = Channel.fromList(FastqScan.findFastqDirs(params.reads))
    // | concat_fastq_files

    fastqs = Channel.of(["sample", params.reads])
    | concat_fastq_files
    | map { sampleid, reads -> [sampleid, reads, params.reference]}
    | filter_reads

    // Assembly
    if(params.contigs == null ) {
        assembly_out = fastqs
        | de_novo_assemblies
        assembly = assembly_out.contigs
    } else {
        assembly = Channel.of(["sample", params.contigs])
    }

    // Consensus
    if(params.draft == null) {
        draft_input = assembly
        | join(fastqs, by: 0)
        | map { sampleid, contigs, reads -> [
            "sample": sampleid,
            "draft_name": params.draft_name,
            "contigs": contigs,
            "reference": params.reference,
            "reads": reads,
            "trl_end": params.reference_terminal_repeat_left_end,
            "trr_start": params.reference_terminal_repeat_right_start
        ]}
        // TODO: pass on the sample-ID!
        draft_results = draft_input | draft_genome
        draft = draft_results.draft | map { meta, draft, basecounts -> [draft, basecounts] }
    } else {
        draft_results = null
        draft = Channel.of(params.draft)
    }

    // Annotation
    annotation = draft
    | map { draft, basecounts -> [
        "draft": draft,
        "genbank": params.reference_genbank,
        "basecounts": basecounts
    ] }
    | annotate_draft

    Channel.empty()
    | mix(
        (draft_results ? draft_results.write_this : Channel.empty()),
        annotation.write_this
    )
    | output
}


// TODO:
// - chose reference
// - report + output
