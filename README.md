# Analysis Protocol:

The top 100 codependent genes for a gene of interest (GOI) (i.e. NADK2 or MECR) was uploaded into String-db (version 12.0) and searched against either the *mus musculus* or *Homo sapiens* organism protein database (for NADK2 and MECR respectively). Subsequently, the Markov Clustering (MCL) clustering algorithm was run with default parameters on the identified proteins. A gene count threshold of 5% of input gene was applied to identified clusters. Each of the identified clusters was then separately analyzed and labeled using the Biological Process (Gene Ontology) term. Labels were assiged using a weighted combination of the gene count, strength of the functional enrichment, and false discovery rate. 

## Colored Edges (use with above)
For visualization, edges between the GOI and cluster genes have the widths and lengths determined by the correlation between the two genes. Red edges indicate a positive correlation while blue edges indicate a negative correlation. For within-cluster genes, fixed colors and widths were utilized.

## Gray Edges (use with above)
For visualization, edges between the GOI and cluster genes have the widths and lengths determined by the correlation between the two genes. For within-cluster genes, fixed widths were utilized.

# String-db Permalinks:
|   GOI     | Species | Description | Gene Count | Cluster Label | Link | 
|-----------|---------|-------------|------------|---------------|------|
|   NADK2   | *Homo Sapiens* | All Genes   | 101 | N/A | https://version-12-0.string-db.org/cgi/network?networkId=bxaddMsMZWYC |
|   NADK2   | *Homo Sapiens* | Cluster 1   | 21 | Mitochondrial Gene Expression | https://version-12-0.string-db.org/cgi/network?networkId=bK9FTMyBGoiC |
|   NADK2   | *Homo Sapiens* | Cluster 2   | 11 | Tricarboxylic Acid Cycle | https://version-12-0.string-db.org/cgi/network?networkId=ba5hjiTpsCI6 | 
|   NADK2   | *Homo Sapiens* | Cluster 3   | 6 | Iron-sulfur Cluster Assembly  | https://version-12-0.string-db.org/cgi/network?networkId=bQoOncCYobCP |
|   NADK2   | *Homo Sapiens* | Cluster 4   | 5 | Protein Lipoylation | https://version-12-0.string-db.org/cgi/network?networkId=bEN3eHzU5B4g |
|           |         |    |   | | |
|   NADK2   | *Mus Musculus* | All Genes   | 96 | N/A | https://version-12-0.string-db.org/cgi/network?networkId=bkZNmH8tDuTl |
|   NADK2   | *Mus Musculus* | Cluster 1   | 23 | Mitochondrial Gene Expression | https://version-12-0.string-db.org/cgi/network?networkId=btDCvadN5XzJ |
|   NADK2   | *Mus Musculus* | Cluster 2   | 9 | Tricarboxylic Acid Cycle | https://version-12-0.string-db.org/cgi/network?networkId=bnpEd7ZeQLUy | 
|   NADK2   | *Mus Musculus* | Cluster 3   | 6 | Iron-sulfur Cluster Assembly  | https://version-12-0.string-db.org/cgi/network?networkId=bR3UZreWGQHV |
|   NADK2   | *Mus Musculus* | Cluster 4   | 5 | Protein Lipoylation | https://version-12-0.string-db.org/cgi/network?networkId=bsTVdaUYO480 |
|   NADK2   | *Mus Musculus* | Cluster 5   | 6 | Mitochondrial Transit Peptide  | https://version-12-0.string-db.org/cgi/network?networkId=bvWFy02RJQzv |
|   NADK2   | *Mus Musculus* | Cluster 6   | 5 | Alpha-Amino Acid Metabolic Process | https://version-12-0.string-db.org/cgi/network?networkId=bzXlRGF5Zget |
|           |         |    |   | | |
|   MECR    |*Homo Sapiens* | All Genes   | 99 | N/A | https://version-12-0.string-db.org/cgi/network?networkId=bhMiyevNtv5L |
|   MECR    | *Homo Sapiens* | Cluster 1   | 40 | Mitochondrial Gene Expression | https://version-12-0.string-db.org/cgi/network?networkId=b4HDQC9hGIx1 |
|   MECR    |*Homo Sapiens* | Cluster 2   | 14 | tRNA Aminoacylation | https://version-12-0.string-db.org/cgi/network?networkId=bJiJfjUH5Pu6 | 
|   MECR    |*Homo Sapiens* | Cluster 3   | 6 | tRNA Methylation | https://version-12-0.string-db.org/cgi/network?networkId=brSHEs7Fi3px |
|   MECR    |*Homo Sapiens* | Cluster 4   | 6 | Ubiquinone Biosynthetic Process | https://version-12-0.string-db.org/cgi/network?networkId=bi4Mzo4T7Cn0  |
|   MECR    |*Homo Sapiens* | Cluster 5   | 6 | Protein Lipoylation | https://version-12-0.string-db.org/cgi/network?networkId=bUQHVhMphBwY |
|   MECR    |*Homo Sapiens* | Cluster 6   | 5 | Respiratory Chain Complex IV Assembly | https://version-12-0.string-db.org/cgi/network?networkId=bqCQaoWlgEVf |

