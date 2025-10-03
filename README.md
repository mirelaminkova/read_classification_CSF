# Read-level tumour classification using cfDNA from cerebrospinal fluid

Read-level classification was performed using cfDNA obtained from cerebrospinal fluid. The method was split into two parts: Use of conventional tools and applications of transformer-based model `MethylBERT`. See graphical overview: 
## Method workflow 

![Workflow](graphical_overview/method_overview-2.png)

### Methylation calling and motif enrichment analysis 

![Analysis](graphical_overview/sample_split-(2).pdf)

Samples were filtered from secondary and low-quality reads. Further, to speed up analysis and save on computational resources, samples were mixed and split into subsets of 50 thousand reads. Methylation was performed using a custom methylaiton calling script (can be found in scripts/fasta_generate.py). The output, saved in FASTA files is then applied to motif enrichment analysis, where read IDs were appended with the sample names, allowing tracing back to the sample of origin. 

### MethylBERT 
MethylBERT was applied where, DMRs were called using `modkit_pileup`, `BSSEQ` and `DSS` tools. Finetuning was performed, where sequence length was adjusted to 500 bp. To apply MethylBERT, please refer to their [Github repository](https://github.com/CompEpigen/methylbert.git)
