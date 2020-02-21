# DNA Read and Sample Embeddings

This work involves first encoding ("embedding") each sequence into a dense, 
low-dimensional, numeric vector space. We use Skip-Gram word2vec to embed 
k-mers, obtained from 16S rRNA amplicon surveys, and then leverage an 
existing sentence embedding technique to embed all sequences belonging to 
specific samples. Our work demonstrated that these representations are 
meaningful, and hence the embedding space can be exploited as a form of 
feature extraction for exploratory analysis. We showed that sequence 
embeddings preserve relevant information about the sequencing data such as 
k-mer context, sequence taxonomy, and sample class. In addition, embeddings 
are versatile features that can be used for many downstream tasks, such as 
taxonomic and sample classification. 

[Stephen Woloszynek, Zhengqiao Zhao, Jian Chen, and Gail L. Rosen. 16S rRNA 
sequence embeddings: Meaningful numeric feature representations of 
nucleotide sequences that are convenient for downstream analyses. 2019. PLOS
Computational Biology. 15(2). doi: 10.1371/journal.pcbi.1006721](https://doi.org/10.1371/journal.pcbi.1006721)

## Workflow

Code demonstrating the workflow with working examples can be found  at 
[github.com/EESI/16s\_embeddings](https://github.com/EESI/16s_embeddings).

## Datasets

Complete datasets used in the manuscript can be found at 
[gitlab.com/sw1/embeddings\_data](https://gitlab.com/sw1/embeddings_data).
