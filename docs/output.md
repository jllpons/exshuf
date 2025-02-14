```
exshuf_results/
├──pipeline_info                        # Pipeline logs, error messages, and run summaries
│   ├── pipeline.log
│   ├── parameters.yaml
│   └── versions.yaml
├── per_protein/                        # Subdirectories for each protein/transcript analysis
│   ├── protein1/             
│   │   ├── intronic_hits.gff       
│   │   ├── exonic_hits.gff             
│   │   ├── opposite_strand_hits.gff    
│   │   └─ intermediates/               # Intermediate files for debugging or further analysis
│   │       ├── gene.gff
│   │       ├── exons.gff
│   │       └── introns.gff
│   ├── protein2/
│   │   ├── intronic_hits.gff
│   │   ├── exonic_hits.gff
│   │   ├── opposite_strand_hits.gff
│   │   └─ intermediates/
│   │       ├── gene.gff
│   │       ├── exons.gff
│   │       └── introns.gff
│   └── ... (one directory per protein)
├── deactivated/
│       ├── hitsDeactivated.gff
│       └─ intermediates/
│           └── combined.gff
└── summaries/                          # Final summary tables and reports
    ├── candidate_events.csv            # Aggregated table of candidate events
    ├── candidate_events.pdf            # (Optional) a visual report or summary plots
    └── overview.txt                    # A summary report with key statistics and notes
```


## Aggregated Summary Table

```tsv
Gene_ID	Transcript_ID	Protein_ID	Chromosome	Start	End	Strand	Event_Type	Overlap_Length	HMM_Hit_ID	HMM_Score
GENE1	TR1	P1	chr1	12345	12500	+	Intronic	56	hit_001	52.3
GENE1	TR1	P1	chr1	13000	13100	+	Opposite_Strand	100	hit_007	47.1
GENE2	TR2	P2	chr2	22345	22420	-	Deactivated	75	hit_015	60.0
```
