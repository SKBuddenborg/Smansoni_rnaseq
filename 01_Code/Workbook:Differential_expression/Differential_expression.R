# Analysis of Schistosoma mansoni bulk RNAseq
```R
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install('DESeq2')
BiocManager::install('PCAtools')
BiocManager::install('pheatmap')
BiocManager::install('RColorBrewer')
BiocManager::install('gplots')
BiocManager::install('ggplot2')

library(DESeq2)
library(PCAtools)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gplots)

directory <- ('/Users/skb/Desktop/Smansoni_rnaseq/HTSeq')
sampleFiles <- grep ("htseq.o", list.files(directory),value=TRUE)
sampleFiles
sampleCondition <- gsub("_R.htseq.o","",sampleFiles)
unique(sampleCondition)
```
[1] "1d_Sporocysts"    "2d_Somules"       "32d_Sporocysts"   "5d_Sporocysts"    "Cercariae"
[6] "Eggs"             "F_26d_juveniles"  "F_2d_Somules"     "F_32d_Sporocysts" "F_Cercariae"
[11] "M_26d_juveniles"  "M_2d_Somules"     "M_32d_Sporocysts" "M_Cercariae"      "Miracidia"

```R
sampleTable <- data.frame(sampleName=sampleFiles,
                          fileName=sampleFiles,
                          condition=sampleCondition)
sampleTable
```
  sampleName                   fileName        condition
1     1d_Sporocysts_R1    1d_Sporocysts_R1htseq.o    1d_Sporocysts
2     1d_Sporocysts_R2    1d_Sporocysts_R2htseq.o    1d_Sporocysts
3     1d_Sporocysts_R3    1d_Sporocysts_R3htseq.o    1d_Sporocysts
4     1d_Sporocysts_R4    1d_Sporocysts_R4htseq.o    1d_Sporocysts
5     1d_Sporocysts_R5    1d_Sporocysts_R5htseq.o    1d_Sporocysts
6        2d_Somules_R1       2d_Somules_R1htseq.o       2d_Somules
7        2d_Somules_R2       2d_Somules_R2htseq.o       2d_Somules
8        2d_Somules_R3       2d_Somules_R3htseq.o       2d_Somules
9        2d_Somules_R4       2d_Somules_R4htseq.o       2d_Somules
10       2d_Somules_R5       2d_Somules_R5htseq.o       2d_Somules
11   32d_Sporocysts_R1   32d_Sporocysts_R1htseq.o   32d_Sporocysts
12   32d_Sporocysts_R2   32d_Sporocysts_R2htseq.o   32d_Sporocysts
13   32d_Sporocysts_R3   32d_Sporocysts_R3htseq.o   32d_Sporocysts
14   32d_Sporocysts_R4   32d_Sporocysts_R4htseq.o   32d_Sporocysts
15   32d_Sporocysts_R5   32d_Sporocysts_R5htseq.o   32d_Sporocysts
16    5d_Sporocysts_R1    5d_Sporocysts_R1htseq.o    5d_Sporocysts
17    5d_Sporocysts_R2    5d_Sporocysts_R2htseq.o    5d_Sporocysts
18    5d_Sporocysts_R3    5d_Sporocysts_R3htseq.o    5d_Sporocysts
19    5d_Sporocysts_R4    5d_Sporocysts_R4htseq.o    5d_Sporocysts
20    5d_Sporocysts_R5    5d_Sporocysts_R5htseq.o    5d_Sporocysts
21        Cercariae_R1        Cercariae_R1htseq.o        Cercariae
22        Cercariae_R2        Cercariae_R2htseq.o        Cercariae
23        Cercariae_R3        Cercariae_R3htseq.o        Cercariae
24        Cercariae_R4        Cercariae_R4htseq.o        Cercariae
25        Cercariae_R5        Cercariae_R5htseq.o        Cercariae
26             Eggs_R1             Eggs_R1htseq.o             Eggs
27             Eggs_R2             Eggs_R2htseq.o             Eggs
28             Eggs_R3             Eggs_R3htseq.o             Eggs
29             Eggs_R4             Eggs_R4htseq.o             Eggs
30             Eggs_R5             Eggs_R5htseq.o             Eggs
31  F_26d_juveniles_R1  F_26d_juveniles_R1htseq.o  F_26d_juveniles
32  F_26d_juveniles_R2  F_26d_juveniles_R2htseq.o  F_26d_juveniles
33  F_26d_juveniles_R3  F_26d_juveniles_R3htseq.o  F_26d_juveniles
34  F_26d_juveniles_R4  F_26d_juveniles_R4htseq.o  F_26d_juveniles
35  F_26d_juveniles_R5  F_26d_juveniles_R5htseq.o  F_26d_juveniles
36     F_2d_Somules_R1     F_2d_Somules_R1htseq.o     F_2d_Somules
37     F_2d_Somules_R2     F_2d_Somules_R2htseq.o     F_2d_Somules
38     F_2d_Somules_R3     F_2d_Somules_R3htseq.o     F_2d_Somules
39     F_2d_Somules_R4     F_2d_Somules_R4htseq.o     F_2d_Somules
40     F_2d_Somules_R5     F_2d_Somules_R5htseq.o     F_2d_Somules
41 F_32d_Sporocysts_R1 F_32d_Sporocysts_R1htseq.o F_32d_Sporocysts
42 F_32d_Sporocysts_R2 F_32d_Sporocysts_R2htseq.o F_32d_Sporocysts
43 F_32d_Sporocysts_R3 F_32d_Sporocysts_R3htseq.o F_32d_Sporocysts
44 F_32d_Sporocysts_R4 F_32d_Sporocysts_R4htseq.o F_32d_Sporocysts
45 F_32d_Sporocysts_R5 F_32d_Sporocysts_R5htseq.o F_32d_Sporocysts
46      F_Cercariae_R1      F_Cercariae_R1htseq.o      F_Cercariae
47      F_Cercariae_R2      F_Cercariae_R2htseq.o      F_Cercariae
48      F_Cercariae_R3      F_Cercariae_R3htseq.o      F_Cercariae
49      F_Cercariae_R4      F_Cercariae_R4htseq.o      F_Cercariae
50      F_Cercariae_R5      F_Cercariae_R5htseq.o      F_Cercariae
51  M_26d_juveniles_R1  M_26d_juveniles_R1htseq.o  M_26d_juveniles
52  M_26d_juveniles_R2  M_26d_juveniles_R2htseq.o  M_26d_juveniles
53  M_26d_juveniles_R3  M_26d_juveniles_R3htseq.o  M_26d_juveniles
54  M_26d_juveniles_R4  M_26d_juveniles_R4htseq.o  M_26d_juveniles
55  M_26d_juveniles_R5  M_26d_juveniles_R5htseq.o  M_26d_juveniles
56     M_2d_Somules_R1     M_2d_Somules_R1htseq.o     M_2d_Somules
57     M_2d_Somules_R2     M_2d_Somules_R2htseq.o     M_2d_Somules
58     M_2d_Somules_R3     M_2d_Somules_R3htseq.o     M_2d_Somules
59     M_2d_Somules_R4     M_2d_Somules_R4htseq.o     M_2d_Somules
60     M_2d_Somules_R5     M_2d_Somules_R5htseq.o     M_2d_Somules
61 M_32d_Sporocysts_R1 M_32d_Sporocysts_R1htseq.o M_32d_Sporocysts
62 M_32d_Sporocysts_R2 M_32d_Sporocysts_R2htseq.o M_32d_Sporocysts
63 M_32d_Sporocysts_R3 M_32d_Sporocysts_R3htseq.o M_32d_Sporocysts
64 M_32d_Sporocysts_R4 M_32d_Sporocysts_R4htseq.o M_32d_Sporocysts
65 M_32d_Sporocysts_R5 M_32d_Sporocysts_R5htseq.o M_32d_Sporocysts
66      M_Cercariae_R1      M_Cercariae_R1htseq.o      M_Cercariae
67      M_Cercariae_R2      M_Cercariae_R2htseq.o      M_Cercariae
68      M_Cercariae_R3      M_Cercariae_R3htseq.o      M_Cercariae
69      M_Cercariae_R4      M_Cercariae_R4htseq.o      M_Cercariae
70      M_Cercariae_R5      M_Cercariae_R5htseq.o      M_Cercariae
71        Miracidia_R1        Miracidia_R1htseq.o        Miracidia
72        Miracidia_R2        Miracidia_R2htseq.o        Miracidia
73        Miracidia_R3        Miracidia_R3htseq.o        Miracidia
74        Miracidia_R4        Miracidia_R4htseq.o        Miracidia
75        Miracidia_R5        Miracidia_R5htseq.o        Miracidia
