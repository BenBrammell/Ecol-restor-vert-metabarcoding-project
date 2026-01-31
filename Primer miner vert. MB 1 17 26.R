# Feel free to contact Vasco per mail if you run into issues (luckylion07@googlemail.com). Enjoy!

# set the path to the PrimerMinder folder you just downloaded
setwd("~/Desktop/PrimerMiner")

# automatically downlaod & install the PrimerMiner package icl dependencies
install.packages("devtools")
library("devtools")
install_github("boldsystems-central/BOLDconnectR")
install_github("VascoElbrecht/PrimerMiner", subdir="PrimerMiner")

# load the package into R
library("PrimerMiner")

# quick guide inside R, further documentation in the Github Wiki
# https://github.com/VascoElbrecht/PrimerMiner/wiki
browseVignettes("PrimerMiner")


# Set path to sample data
setwd("Sample_Data16")

# creating configuration file and batch downloading reads
batch_config("config.txt")

# batch download and process sequence data
batch_download("taxa_small.csv", "config.txt")

#You have to generate and manually check the alignment! Only extract the region amplified by the primers (+2 * the primer length on each side). We recommend Geneious for doing this!

# You find 5 sample alignments generated with Geneious in the folder "1 COI alignments (unprocessed)". We will use these in the next steps of this tutorial.

# get alignments to process
fastafiles <- list.files("1 COI alignments (unprocessed)", full.names=T)

# define name and folder to save files in
fastafiles_export <- paste("2 COI alignments (processed)", list.files("1 COI alignments (unprocessed)"), sep="/")

# process files! This function will remove gaps from the alignments and and apply selective trimming to the primer rgions
for (i in 1:length(fastafiles)){
selectivetrim(fastafiles[i], fastafiles_export[i], trimL=25, trimR=26, gaps=0.10, minsequL=28)
}


# The processed gap free alignemts can now be plotted with PrimerMiner or processed with third party primer development software
# Make sure to chekc that all alignemts have the same length (are aligned)
alignments <- list.files("2 COI alignments (processed)", full.name=T) # find files

pdf("plot_alignments.pdf", height=4, width=100)
plot_alignments(alignments, Order_names=gsub(".*./._(.*)_.*", "\\1", alignments))
dev.off()




# in silico primer evaluation
# With PrimerMiner v0.13 the scorring tables for missmatch possition and type are integrated in PrimerMiner. You can however generate your defult tables and use them (provide them as a csv, see sample data for the default tables!)

# Batra F evaluation (1/16/26)
evaluate_primer("primer_scoring/Batra_F_first_half.fasta", "ACACCGCCCGTCACCCT", 969, 985, save="save_evaluation_table_LCO16Z.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1") 

evaluate_primer("primer_scoring/Batra_F_second_half.fasta", "ACACCGCCCGTCACCCT", 995, 1011, save="save_evaluation_table_LCO16ZA.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1") 

# Batra R evaluation (1/16/26)
evaluate_primer("primer_scoring/Batra_R_first_half.fasta", "AAGTCGTAACATGGTAAGTRTAC", 1042, 1064, save="save_evaluation_table_LCO16ZD.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1") 

evaluate_primer("primer_scoring/Batra_R_second_half.fasta", "AAGTCGTAACATGGTAAGTRTAC", 1129, 1151, save="save_evaluation_table_LCO16ZH.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1") 

# MiFish F evaluation (1/17/26)
evaluate_primer("primer_scoring/MiFish_F.fasta", "GTCGGTAAAACTCGTGCCAGC", 304, 324, save="save_evaluation_table_LCO16ZZ.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1") 

# MiFish R evaluation (1/17/26)
evaluate_primer("primer_scoring/MiFish_R_first_half_1-101.fasta", "CAAACTGGGATTAGATACCCCACTA", 528, 552, save="save_evaluation_table_118.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1")

evaluate_primer("primer_scoring/MiFish_R_second_half_102.fasta", "CAAACTGGGATTAGATACCCCACTA", 534, 558, save="save_evaluation_table_118B.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1")

#evaluate_primer("primer_scoring/15BBatra_F.fasta", ""/media/sf_Shared/Batra15.fasta"", 1016, 1038, save="save_evaluation_table_BR1.csv", mm_position ="Position_v1", adjacent=2, mm_type="Type_v1", forward=F) 

# Evaluate primer paird against each others
primer_threshold("save_evaluation_table_LCO.csv", "save_evaluation_table_BR1.csv", threshold=120)



# generate different primer versions that are contained in a degenerated primer
primerversions("GCHCCHGAYATRGCHTTYCC") # BF2

getwd()
