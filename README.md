# Read depth and balance trees

> This repository holds the code for my PhD dissertation, which is a comparison of balance tree methods across different datasets and a study on how read depth affects the downstream results.

## Outline of README


## Datasets used in this project:
> All the datasets in the project are publicly available 16S datasets
### Noguera-Julian
  - Accession: [PRJNA307231](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA307231)
### Fodor 
  - Accession: [PRJNA391915](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA391915)
### Vangay
  - Accession: [PRJEB28687](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB28687)
### Jones
  - Accession: [PRJNA397450](https://www.ebi.ac.uk/ena/browser/view/PRJNA397450)

## Transformations used in this project

## For building the IQTREE trees:
  IQTREE requires aligned sequences, therefore we first must align the sequences with ClustalW.
### ClustalW
  Using input from the build fasta step from `p1_make_asv_fasta.R`
`/path/to/clustalo -in input.fasta -out output.fasta -fasta`

### Command for making IQTREE
`iqtree -s example.phy -T AUTO`
***For the Jones dataset, the GTR+F+R5 was selected because it had the lowest BIC after the modelfinder timing out at 48 hours.
 No. 	Model         	-LnL         	df  AIC          	AICc       	 BIC
8	GTR+F+R5      	1688897.219	56065 3489924.438	6290170504.438	3806380.010


   
   