# ICR_CrucivirusEvolution
This project examines the use of alignment-free analyses in understanding the evolution of crucivirus.

We employed the following analyses:
1. Genome-sense prediction
We found the codon-ending preferences of each Capsid and Rep gene of our datasets. We compared the preferences of the Capsid and Rep gene of the same genome. If they have the same codon bias, the genome is unisense. Otherwise, it is ambisense. Horizontal gene transfers, which could possibily have occured between viruses and changed the composition of gene sequences of crucivirus, could have affected the sense of the crucivirus. Thus, by predicting the sense of each genome, we could observe a pattern of horizontal gene transfer.
2. Pairwise comparison
We used pairwise analysis to compare the kmer composition in between crucivurs capsids and reps. We calculated the percentage of matching kmer sequences of differing lengths of all capsid genes with all capsid genes (including to themselves) and the same for rep genes. These data are visualized using heatmaps and dendrograms. The dendrograms give predictions of the relationship and evolution of capsid genes and rep genes in crucivirus. We can analyze the pattern of evolution of capsid and rep to determine if the two genes have coexisted in the same crucivirus for a long period of time.
