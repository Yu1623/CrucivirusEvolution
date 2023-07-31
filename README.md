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
The following command line gives you a table containing the genome names, the codon-ending preference of CP gene, the codon-ending preference of Rep gene, and the sense of the genome.
It also gives you the number of unexpected genome senses and a list of the names of the genomes.
```
python3 checkGenomeSensePrediction.py
```

2. Pairwise comparison

We used pairwise analysis to compare the kmer composition in between crucivurs capsids and reps. We calculated the percentage of matching kmer sequences of differing lengths of all capsid genes with all capsid genes (including to themselves) and the same for rep genes. These data are visualized using heatmaps and dendrograms. The dendrograms give predictions of the relationship and evolution of capsid genes and rep genes in crucivirus. We can analyze the pattern of evolution of capsid and rep to determine if the two genes have coexisted in the same crucivirus for a long period of time.

Procedure:
The following command line gives you two heat maps and two dendrograms. One heat map and one dendrogram show the pairwise comparison among crucivirus CP genes. The other heat map and dendrogram show the pairwise comparison among crucivirus Rep genes.
```
python3 VisualPairwise.py
```
3. Kmer rank for the similarity of CP and Rep genes

We proposed that the more similar two genes are, the higher the probability that they have been in the same virus for a long period of time, on the scale of evolution. In order to find the similarity between CP and Rep genes of each virus type, we compared each kmer sequence of CP gene with each kmer sequence of Rep gene with the same length. After we recheived the percentages for the similarity between CP and Rep are found for all genomes of each virus type, we made a ranking for each virus type. The ranking has the genome with the highest percentage of similarity between its CP and Rep genes on top and the genome with the lowest percentage on the bottom. Line graphs are used to visualize the rankings, with the percents of similarity between the two genes against the counts at which the percents are found in the genomes. We compiled the graphs for crucivirus, DNA, and RNA to see patterns among the three genome types.

Procedure:

Step 1: This command line gives you ranks of the percentages of similar k-mers between CP and Rep genes of crucivirus, RNA, and DNA genomes.
```
python3 CPRepKmerRank.py
```

Step 2: This command line gives you line plots of the number of percentages of similar k-mers between CP and Rep and their corresponding count.
```
python3 KmerRankCPRepVisual.py
```

4. Comparison of kmer sequences between genomes

Kmer sequences could be used to see horizontal genetic transfer between genomes. When comparing two genomes, kmer sequences that appear frequently in one genome and less frequently in the second genome often infer to a transfer of nucleotides from the first genome to the second. We want to see if there are horizontal genetic transfers between DNA and crucivirus and where the genetic transfers are located. We first created a heatmap showing the percentage of similar kmer sequences between each DNA genome and each crucivirus. We picked the genome pairs that have the least similarity and the most similarity to create charts visualizing the counts of every shared kmer sequences. We compared the charts of the most and least similar genomes to find patterns that distinguishes the most similar from the least. This gives confidence to if the kmer sequences that appeared at a very high count for the most similar genomes possibly are due to horizontal genetic transfers. The count of kmer sequences for the two genomes separately allows us to make predictions of the direction of the horizontal genetic transfer. We may also repeat this analysis for RNA and crucivirus.

Procedure:

Step 1: This command line makes a heat map with each cell representing the percentage of similar k-mers between a pair of crucivirus and DNA genomes.
```
python virusMatchingKmer.py
```

Step 2: This command line gives you four line plots, each showing the shared k-mers between the two genomes of the crucivirus-DNA pair and the number of times each k-mer appeared in both genomes. The first two are for genome pairs with high similarity and the last two are for genome pairs with low similarity.
```
python3 genomesMatchingKmers.py
```

This project is licensed under GNU General Public license Version 3.
