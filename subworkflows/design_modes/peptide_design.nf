include { BINDER_DESIGN } from '../../modules/local/binder_design/main'
include { BINDER_SUMMARY } from '../../modules/local/binder_summary/main'
include { LIGANDMPNN_DESIGN } from '../../modules/local/ligandmpnn_design/main'
include { MERGE_BINDER_RESULTS } from '../../modules/local/merge_binder_results/main'
include { BIOPHYSICAL_FILTER } from '../../modules/local/biophysical_filter/main'
include { FINALIZE_PEPTIDE_RESULTS } from '../../modules/local/finalize_peptide_results/main'
workflow PEPTIDE_DESIGN {

    take:
    target_pdb

    main:

    BINDER_DESIGN(target_pdb)

    binder_pdb_ch = BINDER_DESIGN.out[0].flatten()

    BINDER_SUMMARY(BINDER_DESIGN.out[0])

    LIGANDMPNN_DESIGN(binder_pdb_ch)

    MERGE_BINDER_RESULTS(
        BINDER_SUMMARY.out,
        LIGANDMPNN_DESIGN.out[0].collect()
    )

    BIOPHYSICAL_FILTER(MERGE_BINDER_RESULTS.out)

    FINALIZE_PEPTIDE_RESULTS(BIOPHYSICAL_FILTER.out)
    emit:
    final_peptide_candidates = FINALIZE_PEPTIDE_RESULTS.out
}
