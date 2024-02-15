# The-ZmHSF20-ZmHSF4-ZmCesA2-module-regulates-heat-stress-tolerance-in-maize
## SNK for upstream transcriptome analysis

```
# Upstream transcriptome analysis using snakemake. Output: 06_transcripts_quant; 05_stringtie_merged/stringtie_merged.gtf
snakemake -j 00_Snakemake_rna-seq.smk

# 00_extract_exp.py was used to extract FPKM values for all genes of different samples
python3 00_extract_exp.py <input_path> <output_path>

# 00_ID_transition.py was used for the gene_id modification
python3 00_ID_transition.py <merge_gtf> <input_path> <output_path>

# 00_rm_MGSRT.py was used for remove unpair gene
python3 00_rm_MGSRT.py <input_paht> <output_path>
```

## DESeq2 for analyzing differential genes
## Weighted co-expression network analysis (WGCNA)
## GO enrichment analysis
