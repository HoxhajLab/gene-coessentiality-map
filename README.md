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

# Installation:

1) Clone the github: 
```
git clone https://github.com/HoxhajLab/gene-coessentiality-map.git
```
2) Download the  zip
- Navigate to: https://github.com/HoxhajLab/gene-coessentiality-map/tree/main and download the code as a zip file. 
- Unpack the zip into the local directory of your choice.

Installation should take no longer than 2-3 minutes on standard hardware.

# Demo:
Data and code necessary to recreate the experimental results is included in the repository. In order to reproduce the experimental results, simply run the jupyter notebook `visualization.ipynb`. Expected ouput are included in the manuscript, the supplemental text, as well as `fig` folder in the current repository.

Reproducing the results should take no longer than 5 minutes using standard hardware.

# System Requirements:
-----
* Python 3.9.18 
* Linux-3.10.0-1160.88.1.el7.x86_64-x86_64-with-glibc2.17

### Essential packages:
|Package | Version |
|--------|---------|
|matplotlib  |        3.9.4|
pandas         |     1.5.3
session_info    |    1.0.0

### Package dependencies:
|Package | Version |
|--------|---------|
PIL              |   10.0.1
asttokens        |   NA
backcall         |   0.2.0
bottleneck       |   1.3.5
cffi             |   1.15.1
cloudpickle      |   2.2.1
comm             |   0.1.4
cycler           |   0.10.0
cython_runtime   |   NA
dateutil         |   2.8.2
debugpy          |   1.6.7
decorator        |   5.1.1
defusedxml       |   0.7.1
entrypoints      |   0.4
exceptiongroup   |   1.0.4
executing        |   0.8.3
importlib_resources| NA
ipykernel          | 6.25.0
ipython_genutils   | 0.2.0
ipywidgets         | 8.1.1
jedi               | 0.18.1
kiwisolver         | 1.4.4
matplotlib_inline  | 0.1.6
mpl_toolkits       | NA
networkx           | 2.8.7
numexpr            | 2.8.7
numpy              | 1.23.5
packaging          | 23.1
parso              | 0.8.3
pexpect            | 4.8.0
pickleshare        | 0.7.5
pkg_resources      | NA
platformdirs       | 3.10.0
prompt_toolkit     | 3.0.36
psutil             | 5.9.0
ptyprocess         | 0.7.0
pure_eval          | 0.2.2
pyarrow            | 11.0.0
pydev_ipython      | NA
pydevconsole       | NA
pydevd             | 2.9.5
pydevd_file_utils  | NA
pydevd_plugins     | NA
pydevd_tracing     | NA
pygments           | 2.15.1
pyparsing          | 3.0.9
pytz               | 2023.3.post1
six                | 1.16.0
stack_data         | 0.2.0
tornado            | 6.3.3
traitlets          | 5.7.1
typing_extensions  | NA
vscode             | NA
wcwidth            | 0.2.5
zipp               | NA
zmq                | 23.2.0
zoneinfo           | NA
IPython             |8.15.0
jupyter_client      |7.4.9
jupyter_core        |5.3.0
notebook            |6.5.4


