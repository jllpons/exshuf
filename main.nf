import groovy.yaml.YamlBuilder

include { EXSHUF } from './workflows/exshuf'


help_message = """
exshuf.nf - Exon Shuffling Event Detection
==========================================

Detection of potential exon shuffling or events

Usage:
    nextflow run main.nf --features <features.gff> --hits <hits.gff> [options]
    nextflow run main.nf -params-file <yaml> [options]

Required Arguments:
    --features                  : Path to the GFF file containing gene models (e.g., features.gff)
    --hits                      :  Path to the GFF file with translated genome hits against (a/some) HMM (e.g., hits.gff)


Alternative to Required Arguments:
    -params-file <yaml>         : YAML/JSON file with the parameters
                                  Mutually exclusive with --features and --hits

Optional Arguments:
    --kw_type_gene <kw>         : Keyword to identify gene features in the GFF file [default: ${params.kw_type_gene}]

    --kw_type_exon <kw>         : Keyword to identify exon features in the GFF file [default: ${params.kw_type_exon}]

    --kw_type_case_sensitive    : Use case-sensitive keyword matching [default: ${params.kw_type_case_sensitive}]

    --outdir <outdir>           : Output directory [default: ${params.outdir}]

    --help                      : Print help message and exit

    --version                   : Print version and exit
"""


init_summary = """
E X S H U F   P I P E L I N E   v${params.manifest.version}
======================================
features                : ${params.features}
hits                    : ${params.hits}
kw_type_gene            : ${params.kw_type_gene}
kw_type_exon            : ${params.kw_type_exon}
kw_type_case_sensitive  : ${params.kw_type_case_sensitive}
outdir                  : ${params.outdir}

--

Run as                  : ${workflow.commandLine}
Started at              : ${workflow.start}
Config files            : ${workflow.configFiles}

--
"""
// container images : ${workflow.containerEngine}:${workflow.container}


// DESC: Validate input arguments and initialize pipeline, printing a small summary
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: `$init_summary` as log message at `INFO` level
//       `$help_message` as stdout if `--help` flag is set
//       `$version` as stdout if `--version` flag is set
//       Proper error message and exit code if required arguments are missing
// RETS: None
def validateParams() {

    // `--help` and `--version` flags
    if (params.help) {
        println help_message
        System.exit(0)
    }
    if (params.version) {
        println "${params.manifest.name} v${params.manifest.version}"
        System.exit(0)
    }

    // Check required arguments
    if (params.features == null) {
        println help_message
        log.error "Missing required argument: --features"
        System.exit(1)
    }
    if (!file(params.features).exists()) {
        log.error "File not found: ${params.features}"
        System.exit(1)
    }
    if (params.hits == null) {
        println help_message
        log.error "Missing required argument: --hits"
        System.exit(1)
    }
    if (!file(params.hits).exists()) {
        log.error "File not found: ${params.hits}"
        System.exit(1)
    }

}


// DESC: Display completion message based on workflow status
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: Completion message at `INFO` or `ERROR` level
// RETS: None
def completionMsg() {

    if (workflow.success) {
        if (workflow.stats.ignoredCount == 0) {
            log.info "Pipeline completed successfully!"
        }
        else {
            log.info "Pipeline completed successully, but with errored processes"
        }
    }
    else {
        log.error "Pipeline completed with errors"
    }

}

// DESC: Dump parameters to a YAML file
// ARGS: None, uses variables defined at the beginning of the script
// OUTS: YAML file with the parameters
// RETS: None
def dumpParametersToYaml() {
    def paramMap = [
        features: params.features,
        hits: params.hits,
        kw_type_gene: params.kw_type_gene,
        kw_type_exon: params.kw_type_exon,
        kw_type_case_sensitive: params.kw_type_case_sensitive,
        outdir: params.outdir
    ]

    def yaml = new YamlBuilder()
    yaml(paramMap)

    def outputDir = new File("${params.outdir}/pipeline_info")
    outputDir.mkdirs()
    def outputFile = new File(outputDir, "params_used.yaml")
    outputFile.text = yaml.toString()
}


// Main workflow
workflow {

    main:

    // Validate input parameters
    validateParams()
    // Initialization Summary - Everything looks good so far
    log.info init_summary

    // Dump parameters to a YAML file
    dumpParametersToYaml()

    ch_versions = Channel.empty()
    EXSHUF(
        params.features,
        params.hits,
        params.kw_type_gene,
        params.kw_type_exon,
        params.kw_type_case_sensitive,
        ch_versions,
    )
    ch_versions = ch_versions.mix(EXSHUF.out.versions)
    ch_versions.collectFile(
        storeDir: "${params.outdir}/pipeline_info/",
        name: 'versions.yml',
        sort: true,
        newLine: true
    )

    // Display any error encountered during the workflow
    workflow.onComplete {
        completionMsg()
    }

}
