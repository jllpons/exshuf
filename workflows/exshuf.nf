include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_HITS_W_EXONS     } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_HITS_W_INTRONS   } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_HITS_W_GENE_OPP  } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_HITS_DEACTIVATED } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_SUBTRACT                                         } from '../modules/nf-core/bedtools/subtract/main'
include { GAWK_GET_IDS                                              } from '../modules/local/gawk/main'
include { GAWK_MATCH_KWORDS_AND_ID                                  } from '../modules/local/gawk/main'
include { AGGR_EXSHUF_GFFS                                          } from '../modules/local/python/main'


workflow EXSHUF {

    take:
    features
    hits
    kw_type_gene
    kw_type_exon
    kw_type_case_sensitive
    ch_versions


    main:

    ch_features = Channel.fromPath(features)
    ch_hits = Channel.fromPath(hits)


    GAWK_GET_IDS(
        ch_features,
        kw_type_gene,
        kw_type_case_sensitive
    )
    // TODO: Here we could add a GREP module that would allow to filter the IDs based
    //       on a provided
    // The following will create a channel where each item is a single line of the file (ID in this case)
    ch_ids = GAWK_GET_IDS.out.ids_file
        .splitText()
        .map { it.trim() }



    GAWK_MATCH_KWORDS_AND_ID(
        ch_ids,
        ch_features,
        kw_type_gene,
        kw_type_exon,
        kw_type_case_sensitive
    )

    // In this case, this Channel will contain an array of 3 elements:
    // 1. A Groovy map containing the ID of the gene as a meta field for the nf-core process
    // 2. The GFF line containing the coordinates of the gene
    // 3. The GFF line containing the coordinates of the exons
    // <https://nf-co.re/modules/bedtools_subtract>
    ch_subtract_input = GAWK_MATCH_KWORDS_AND_ID.out.results
        .map { id, gene_gff, exons_gff ->
            def meta = [id: id]
            return [meta, gene_gff, exons_gff]
        }

    BEDTOOLS_SUBTRACT(
        ch_subtract_input,
    )


    ch_intersect_w_exons_input = GAWK_MATCH_KWORDS_AND_ID.out.results
        .map { id, gene_gff, exons_gff ->
            def meta = [id: id]
            return [meta, exons_gff]
        }
        .combine(ch_hits)
        .map { meta, exons_gff, hits ->
            return [meta, hits, exons_gff]
        }

    ch_intersect_w_introns_input = BEDTOOLS_SUBTRACT.out.bed.combine(ch_hits)
        .map { meta, bed, hits ->
            return [meta, hits, bed]
        }

    ch_intersect_w_gene_opp_input = GAWK_MATCH_KWORDS_AND_ID.out.results
        .map { id, gene_gff, exons_gff ->
            def meta = [id: id]
            return [meta, gene_gff]
        }
        .combine(ch_hits)
        .map { meta, gene_gff, hits ->
            return [meta, hits, gene_gff]
        }

    ch_intersect_deactivated_input = ch_features.combine(ch_hits)
        .map { features, hits ->
            def meta = [id: "Deactivated Hits"]
            return [meta, hits, features]
        }


    BEDTOOLS_INTERSECT_HITS_W_EXONS(
        ch_intersect_w_exons_input,
        [[], []] // Empty array to avoid error for not passing in the chromosome size file
    )

    BEDTOOLS_INTERSECT_HITS_W_INTRONS(
        ch_intersect_w_introns_input,
        [[], []]
    )

    BEDTOOLS_INTERSECT_HITS_W_GENE_OPP(
        ch_intersect_w_gene_opp_input,
        [[], []]
    )

    BEDTOOLS_INTERSECT_HITS_DEACTIVATED(
        ch_intersect_deactivated_input,
        [[], []]
    )

    ch_aggr_input = BEDTOOLS_INTERSECT_HITS_W_EXONS.out.intersect
        .join(BEDTOOLS_INTERSECT_HITS_W_INTRONS.out.intersect, by: [0])
        .join(BEDTOOLS_INTERSECT_HITS_W_GENE_OPP.out.intersect, by: [0])
        .map { meta, exons, introns, gene_opp ->
            [meta.id, exons, introns, gene_opp]
        }

    AGGR_EXSHUF_GFFS(
        ch_aggr_input
    )

    AGGR_EXSHUF_GFFS.out.aggr_data_csv
        .collectFile(
            name: "candidate_events.csv",
            storeDir: "$params.outdir/summaries",
            keepHeader: true
        )


    ch_versions = ch_versions.mix(
        GAWK_GET_IDS.out.versions,
        GAWK_MATCH_KWORDS_AND_ID.out.versions,
        BEDTOOLS_SUBTRACT.out.versions,
        BEDTOOLS_INTERSECT_HITS_W_EXONS.out.versions,
        BEDTOOLS_INTERSECT_HITS_W_INTRONS.out.versions,
        BEDTOOLS_INTERSECT_HITS_W_GENE_OPP.out.versions,
        BEDTOOLS_INTERSECT_HITS_DEACTIVATED.out.versions,
        AGGR_EXSHUF_GFFS.out.versions
    )

    emit:
    introns = BEDTOOLS_SUBTRACT.out.bed
    versions = ch_versions // channel containing path to versions.yml
}
