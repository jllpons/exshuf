process AGGR_EXSHUF_GFFS {
    tag "$id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/python:3.11' :
        'https://depot.galaxyproject.org/singularity/python:3.11' }"

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
