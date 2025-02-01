#########################
#########################
Processed data ready for replicating training of SNPeBoT model:
Processed_ASB_data.csv
#########################
#########################
-First 30 columns will be the windowed E-scores for reference and alternate sequence
-Column 31 is a binary indicator of a change in motif between reference and alternate sequences
-Columns 32 to 47 (inclusive) are the non-Windowed E-scores for reference and alternate sequence
-Column 48 is the label for the data (gain, loss, or no change in binding)
-Column 49 and 50 are the p-values for the best fimo hits on the reference and alternate sequences
-Column 51 is the ASB id  

#########################
#########################
Allele Specific Binding data from chipseq experiments:
ASB_sequences.tsv
#########################
#########################
If we would like to match sequences to processed data there are only a few columns that are relevant.

example id to search for (chr10-99393234-G-A-runx3)

Corresponding line in ASB_sequences.tsv:
chr10	99393234	G	A	gm12878	runx3	25	27	False	0.4807692307692308	chr10:99393208-99393259	tccaccagcgattgccccacttgacGccgccatcctgggcgaccgcactgg	tccaccagcgattgccccacttgacAccgccatcctgggcgaccgcactgg

the following command would find the line:
grep "{id_position_1}\t{id_position_2}\t{id_position_3}\t{id_position_4}\t" ASB_sequences.tsv | grep -i {id_position_5}
Which on the example (chr10-99393234-G-A-runx3) looks like this:
grep "chr10\t99393234\tG\tA\t" ASB_sequences.tsv | grep -i runx3

Because the sequence file is not deduplicated you may receive multiple matches for the same sequence

#########################
#########################
google colab file for training of a CNN with the processed data:
CNN_training_SNPeBoT.ipynb
#########################
#########################
The colab file contains instructions on running




