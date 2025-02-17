process GAWK_GET_IDS {
    tag "Get IDs"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    path features
    val keyword
    val case_sensitive

    output:
    path 'ids.txt', emit: ids_file
    path 'versions.yml', emit: versions

    script:
    """
    gawk -v keyword="${keyword}" -v case_sensitive="${case_sensitive}" '
        BEGIN {
            FS = "\t"
            IGNORECASE = (case_sensitive == "false") ? 1 : 0
        }

        {
            if (\$3 == keyword) {
                split(\$9, attrs, ";")

                for (i in attrs) {
                    split(attrs[i], attr, "=")

                    if (attr[1] == "ID") {
                        print attr[2]
                        break
                    }
                }
            }
        }
    ' ${features} | sort | uniq > ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}

process GAWK_MATCH_KWORDS_AND_ID {
    tag "$id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.0' }"

    input:
    each id
    path features
    val gene_kw
    val exon_kw
    val case_sensitive

    output:
    tuple val(id), path('gene.gff'), path('exons.gff'), emit: results
    path 'versions.yml', emit: versions

    publishDir "${params.outdir}/per_gene/${id}/intermediates", mode: 'copy', overwrite: true, pattern: '*.gff'

    script:
    """
    gawk -v gene_kw="${gene_kw}" -v id="${id}" -v case_sensitive="${case_sensitive}" '
        BEGIN {
            FS = "\t"
            IGNORECASE = (case_sensitive == "false") ? 1 : 0
        }

        {
            if (\$3 == gene_kw) {
                split(\$9, attrs, ";")

                for (i in attrs) {
                    split(attrs[i], attr, "=")

                    if (attr[2] == id) {
                        print \$0
                        break
                    }
                }
            }
        }' ${features} > gene.gff

    gawk -v exon_kw="${exon_kw}" -v id="${id}" -v case_sensitive="${case_sensitive}" '
        BEGIN {
            FS = "\t"
            IGNORECASE = (case_sensitive == "false") ? 1 : 0
        }

        {
            if (\$3 == exon_kw) {
                split(\$9, attrs, ";")

                for (i in attrs) {
                    split(attrs[i], attr, "=")

                    if (attr[2] == id) {
                        print \$0
                        break
                    }
                }
            }
        }' ${features} > exons.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | awk 'NR==1{print \$3}')
    END_VERSIONS
    """
}
