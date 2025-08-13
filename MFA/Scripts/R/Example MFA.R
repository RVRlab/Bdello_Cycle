source("/data/Osam/Scripts/R_functions_Ed/Osam_MFA_functions_Version_F1.R")
source("/data/Osam/Scripts/R_functions_Ed/Osam_MFA_Install_Packages.R")


#################### Read the data ##########
reference_genome_path <- "/data/Osam/Sequences/ReferenceGenomes/Bdellovibrio/HD100_ori_centered/HD100_ori_center_Rcorrected.fa"
reading_Path1 <- paste0("/data/Osam/Projects/Bdellovibrio/HiC-Third_try/Bdellovibrio_ori_center/FO52/HiC_2_MFA/Depth_FO52_HiC.txt")
savingPath <- "/data/Osam/Projects/Bdellovibrio/Analysis/HiC2MFA/"

genome_length <- get_genome_length(reference_genome_path)
bin <- 1000
s <- 0.05


#################### The MFA ##########
# If you only want the data binned
FO52_bin <- doMFABin(reading_Path1,genome_length,bin)

# If you want the Loess line
FO52_Loess <- doMFALoess(reading_Path1,genome_length,bin,s)

# If you want to get the slope of an MFA plot
window <- 5000
difference <- 3
genomeWindow <- 500

FO52_model <- doMFALinearModel(reading_Path1, genome_length, bin, window, difference, genomeWindow)

# Window, difference, and genomeWindow are variables used to smooth the MFA bins 
# (to transitorily eliminate bins that are different from what the slope trend is)
# by eliminating them, the slope reflects better the general data tendency

# Window indicates that I am doing a rolling average of 5,000 bp to start smoothing data out
# later on, I divide the genome in smaller sections:
# genomeWindow indicates that I am dividing the genome in 500 kb sections
# Within a section, everything that $difference times the standard deviations above and 
# below the mean will be consider an outlier and will be eliminated from the analysis


#### The plotting
ggplot(NULL, aes(position,average)) + 
geom_point(data=FO52_bin, alpha=3/10, color="black") +
geom_line(data=FO52_Loess, color="red") +
geom_line(data=FO52_model, color="blue") +
xlab("Genomic Position") + ylab("Value") +
scale_y_continuous(trans="log2", limits=c(5e-8, 9e-7)) +
ggtitle("T20_WT_(3rd)")
