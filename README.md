# ICR_CrucivirusEvolution
This project examines the use of alignment-free analyses in understanding the evolution of crucivirus.

Our programs require the following installments: matplotlib, scipy, seaborn, plotly, pylab, pandas, and numpy

How to install:
```
python3 -m pip install -U matplotlib
sudo apt-get install python3-scipy
pip install seaborn
pip install plotly
pip install pylab
pip install pandas
pip install numpy
```

We employed the following analyses:
1. Genome-sense prediction

We found the codon-ending preferences of each Capsid and Rep gene of our datasets. We compared the preferences of the Capsid and Rep gene of the same genome. If they have the same codon bias, the genome is unisense. Otherwise, it is ambisense. Horizontal gene transfers, which could possibily have occured between viruses and changed the composition of gene sequences of crucivirus, could have affected the sense of the crucivirus. Thus, by predicting the sense of each genome, we could observe a pattern of horizontal gene transfer. We showed the results on a table, composing of the name of the genome, the codon-ending bias of the CP gene, the codon bias of the Rep gene, and the sense of the genome. The accuracy of our program and the incorrect predictions allows us to find which genomes are incorrrectly predicted. Genetic strands that have coexisted in the same virus accumulate mutations at different speeds. Genetic strands that have been recently included into a virus due to genetic transfer do not have enough time to accumulate the mutations. This causes the recently acquired sequences to exhibit surprising orientations of genes. Out of all of the incorrect predictions, we identified genomes that have nnt percentages and nna percentages that are very different from each other, such as 0.05 apart, for capsid and/or rep genes. These genomes may have recent horizontal genetic transfers, which is a possible reason for why the orientation of the genes are unexpected.

Procedure:

Step 1: If you are interested in the genome-sense predictions, please use the genomeSenseTable() function from geneStrandPrediction.py and input the CP file name and Rep file name of the genomes that you want to predict the sense of. This gives you a table containing the genome names, the codon-ending preference of CP gene, the codon-ending preference of Rep gene, and the sense of the genome. 
```
genomeSenseTable("879CPs.fasta", "855Reps.fasta")
```

Step 2: Please use the genomeSenseSurprising() function from checkGenomeSensePrediction.py and input the csv. file name that contains the actual genome senses, the CP file name, and the Rep file name of the genomes that you want to predict the sense. This gives you the number of genomes that has surprising and unexpected genome sense, mainly those with percents of nnt and nna that are at least 0.05 apart for CP and/or Rep genes, and an array of the genome names.
```
countOfUnexpectedGenomes, genomesWithUnexpectedSense = genomeSenseSurprising("875crucisAnnotations.csv", "879CPs.fasta", "855Reps.fasta")
print(countOfUnexpectedGenomes, genomesWithUnexpectedSense)
```

2. Pairwise comparison

We used pairwise analysis to compare the kmer composition in between crucivurs capsids and reps. We calculated the percentage of matching kmer sequences of differing lengths of all capsid genes with all capsid genes (including to themselves) and the same for rep genes. These data are visualized using heatmaps and dendrograms. The dendrograms give predictions of the relationship and evolution of capsid genes and rep genes in crucivirus. We can analyze the pattern of evolution of capsid and rep to determine if the two genes have coexisted in the same crucivirus for a long period of time.

Procedure:

Step 1: Use the pairwiseVisual() function from VisualPairwise.py by inputting the file names of the CP and Rep genes that you want to do the pairwise comparison with and the length of kmer. Make sure that the CP and Rep genes are of the same virus type. This gives you four graphs: two heatmaps and two dendrograms for CP and Rep separately. 
```
pairwiseVisual("879CPS.fasta", "855Reps.fasta", 7)
```
3. Kmer rank for the similarity of CP and Rep genes

We proposed that the more similar two genes are, the higher the probability that they have been in the same virus for a long period of time, on the scale of evolution. In order to find the similarity between CP and Rep genes of each virus type, we compared each kmer sequence of CP gene with each kmer sequence of Rep gene with the same length. After we recheived the percentages for the similarity between CP and Rep are found for all genomes of each virus type, we made a ranking for each virus type. The ranking has the genome with the highest percentage of similarity between its CP and Rep genes on top and the genome with the lowest percentage on the bottom. Line graphs are used to visualize the rankings, with the percents of similarity between the two genes against the counts at which the percents are found in the genomes. We compiled the graphs for crucivirus, DNA, and RNA to see patterns among the three genome types.

Procedure:

Step 1: Input the CP and Rep file names of the genome data set that you want to predict the sense of and the length of kmer into the writeRank() function in CPRepKmerRank.py. This creats a txt. file containing the rank of percents of shared kmer sequences between CP and Rep from the genome with the highest percent to the genome with the lowest percent. If you want to make graphs for all three virus types, repeat this step three times for crucivirus, RNA virus, and DNA virus.
```
writeRank("879CPs.fasta", "855Reps.fasta", 7, "Cruci")
writeRank("1526RNA_CPs.fasta", "1514RdRPs.fasta", 7, "RNA")
writeRank("270CRESS_CPs.fasta", "305CRESS_Reps.fasta", 7, "DNA")
```

Step 2: If you only want to see the percents of similar kmer sequences between CP and Rep for one virus type, input the file name for the kmer rank file that you made in the previous step into the function kmerCPRepVisual() from KmerRankCPRepVisual.py. If you want to see the percents of similar kmer sequences between CP and Rep for all three virus types, crucivirus, RNA, DNA, input the file names of the crucivirus kmer rank, RNA kmer rank, and DNA kmer rank that you can make by following step 1 into kmerCPRepCompareVisual() function from KmerRankCPRepVisual.py.
```
kmerCPRepVisual(kmerRankCruci)
kmerCPRepVisual(kmerRankRNA)
kmerCPRepVisual(kmerRankDNA)
kmerCPRepCompareVisual(kmerRankCruci, kmerRankRNA, kmerRankDNA)
```

4. Comparison of kmer sequences between genomes

Kmer sequences could be used to see horizontal genetic transfer between genomes. When comparing two genomes, kmer sequences that appear frequently in one genome and less frequently in the second genome often infer to a transfer of nucleotides from the first genome to the second. We want to see if there are horizontal genetic transfers between DNA and crucivirus and where the genetic transfers are located. We first created a heatmap showing the percentage of similar kmer sequences between each DNA genome and each crucivirus. We picked the genome pairs that have the least similarity and the most similarity to create charts visualizing the counts of every shared kmer sequences. We compared the charts of the most and least similar genomes to find patterns that distinguishes the most similar from the least. This gives confidence to if the kmer sequences that appeared at a very high count for the most similar genomes possibly are due to horizontal genetic transfers. The count of kmer sequences for the two genomes separately allows us to make predictions of the direction of the horizontal genetic transfer. We may also repeat this analysis for RNA and crucivirus.

Procedure:

Step 1: Please use the virusKmerVisual() function from virusMatchingKmer.py by inputting the DNA genomes dataset file name and the crucivirus genomes dataset file name and the length of kmer. This makes a heatmap with each cell representing the percent of matching kmer sequences between a pair of DNA and crucivirus genomes. 
```
virusKmerVisual("316CRESS.fasta", "885crucis.fasta", 7)
```

Step 2: Pick the genome pairs with the highest similarity and the lowest simmilarity. The highest similarity is represented by the lightest color and the lowest similarity is the darkest color. For each genome pair, input the file name of the first genome, the file name of the second genome, the genome name of the first genome, the name of the second genome, and the length of kmer into the genomeskmerMatchVisual() function from genomesMatchingKmers.py. This gives you a line graph showing every shared kmer sequence against the count that they showed up in both genomes and an horizontal line representing the average count.
```
genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_88", "AlphaS_KT948075", 7)
genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_88", "AlphaS_KF471057", 7)
genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruCGE_296", "Bacil_MH617605", 7)
genomeskmerMatchVisual("885crucis.fasta", "316CRESS.fasta", "Cruci_CruV_87", "Bacil_AB193315", 7)
```
