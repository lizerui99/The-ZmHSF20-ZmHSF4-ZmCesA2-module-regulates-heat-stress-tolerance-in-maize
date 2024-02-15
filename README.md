# The-ZmHSF20-ZmHSF4-ZmCesA2-module-regulates-heat-stress-tolerance-in-maize
## SNK for upstream transcriptome analysis

```
# Upstream transcriptome analysis using snakemake
snakemake -j 00_Snakemake_rna-seq.smk

# 00_extract_exp.py was used to extract FPKM values for all genes of different samples
python3 00_extract_exp.py <input>
``` 
## DESeq2 for analyzing differential genes
## Weighted co-expression network analysis (WGCNA)
## GO enrichment analysis
