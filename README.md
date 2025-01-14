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
Scripts in deseq2 folder
## Weighted co-expression network analysis (WGCNA)
Scripts in WGCNA.r
## GO and TF family enrichment analysis
### GO enrichment analysis
B73_AGPv4_GSMER.annot and go_class.txt in GO_enrichment
```
rm(list=ls())
library(clusterProfiler)
library(GOplot)
library(tidyverse)

# read annotation file
go_anno <- read.delim('B73_AGPv4_GSMER.annot', header = FALSE, stringsAsFactors = FALSE)
names(go_anno) <- c('gene_id', 'ID')
go_class <- read.delim('go_class.txt', header = FALSE, stringsAsFactors = FALSE)
names(go_class) <- c('ID', 'Description', 'Ontology')
go_anno <- merge(go_anno, go_class, by = 'ID', all.x = TRUE)

# read gene list
gene_select <- read.delim('gene_list.tab', header = T,stringsAsFactors = FALSE)
gene_cluster <- filter(gene_select)

go_rich <- enricher(gene = gene_cluster$gene,
TERM2GENE = go_anno[c('ID', 'gene_id')],
TERM2NAME = go_anno[c('ID', 'Description')],
pvalueCutoff = 0.05,
pAdjustMethod = 'BH',qvalueCutoff = 0.05, maxGSSize = 500)

write.table(go_rich, 'output.tab', sep = '\t', row.names = FALSE, quote = FALSE)

# output image file
pdf("output.pdf",h=5,w=8)
dotplot(go_rich,showCategory=15,decreasing=T,label_format =50,color = "pvalue")
dev.off()
```

### Transcription factor family enrichment analysis
TF-class.txt and TF-anno.txt in GO_enrichment
```
library(clusterProfiler)
library(GOplot)


term2gene <- read.delim('TF-class.txt',header = FALSE, sep = '\t')
term2name <- read.delim('TF-anno.txt',header = FALSE,sep = '\t')


gene <- read.delim('CytoscapeInput-nodes-brown.txt',header = 1,sep = '\t')

x <- enricher(gene$nodeName,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = 0.05,
              pAdjustMethod = 'BH',
              qvalueCutoff = 0.05, 
              maxGSSize = 500)

dotplot(x,showCategory=100,decreasing=T,label_format =50)


write.table(x, 'rich_tf.tab', sep = '\t', row.names = FALSE, quote = FALSE)
pdf("rich_TF.pdf",h=5,w=8)
dotplot(x,showCategory=15,decreasing=T,label_format =50)
dev.off()
```
