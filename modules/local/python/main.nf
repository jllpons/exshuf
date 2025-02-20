process AGGR_EXSHUF_GFFS {
    tag "$id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(id), path(exon_gff), path(intron_gff), path(opp_strand_gff)

    output:
    path 'aggr_data.csv',   emit: aggr_data_csv
    path 'versions.yml',    emit: versions

    script:
    """
    aggr_exshuf_gffs.py --id $id --exon $exon_gff --intron $intron_gff --opposite $opp_strand_gff > aggr_data.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Python: \$(python -V | awk '{print \$2}')
    END_VERSIONS
    """
}

process HIT_TABLE {
    tag "HIT $kind TABLE"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    path aggr_data_csv
    val  kind

    output:
    path 'hit_*_table.csv', emit: hit_table
    path 'versions.yml',    emit: versions

    script:

    switch(kind){
        case 'type':
            """
            cat $aggr_data_csv | cut -d',' -f1,4 | mk_hit_type_table.py - > hit_type_table.csv

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                Python: \$(python -V | awk '{print \$2}')
            END_VERSIONS
            """
        case 'count':
            """
            cat $aggr_data_csv | cut -d',' -f1,5 | mk_hit_count_table.py - > hit_count_table.csv

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                Python: \$(python -V | awk '{print \$2}')
            END_VERSIONS
            """
    }
}
