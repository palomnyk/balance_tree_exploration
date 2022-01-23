#Read depth and balance trees

##This repository holds the code for my PhD dissertation, which is a comparison of balance tree methods across different datasets and a study on how read depth affects the downstream results.

##Dataset notes:
###Noguera-Julian
  Bioproject accession number: PRJNA307231, SRA accession number: SRP068240
  
###Fodor 
Data came from:
  https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=14&WebEnv=MCID_61a58f8b1fed327d1d4da156&o=acc_s%3Aa
  Bioproject accession number: PRJNA391915, SRA accession number: SAMN07277992
  
###Vangay
  Accession: PRJEB28687 https://www.ncbi.nlm.nih.gov/bioproject/PRJEB28687
  
###Jones
Project: PRJNA397450, Secondary Study Accession: SRP114985
  https://www.ebi.ac.uk/ena/browser/view/PRJNA397450?show=reads
	
##For building the IQTREE trees:
  IQTREE requires aligned sequences, therefore we first must align the sequences with ClustalW.
###ClustalW
  Using input from the build fasta step from `p1_make_asv_fasta.R`
```
   /path/to/clustalo -in input.fasta -out output.fasta -fasta
```

###Command for making IQTREE
```
iqtree -s example.phy -T AUTO
```
***For the Jones dataset, the GTR+F+R5 was selected because it had the lowest BIC after the modelfinder timing out at 48 hours.
 No. 	Model         	-LnL         	df  AIC          	AICc       	 BIC
8	GTR+F+R5      	1688897.219	56065 3489924.438	6290170504.438	3806380.010


   
   