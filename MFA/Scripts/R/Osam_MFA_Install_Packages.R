# To recognize DNA sequences
if (!require(seqinr)) install.packages('seqinr')
suppressPackageStartupMessages(library(seqinr))
# To graphs
if (!require(ggplot2)) install.packages('ggplot2')
suppressPackageStartupMessages(library(ggplot2))
# To melt a data.frame into one with categories.
if (!require(reshape2)) install.packages('reshape2')
suppressPackageStartupMessages(library(reshape2))
# To avoid the names of vectors in graphs to be too close that they cant be read. 
if (!require(ggrepel)) install.packages('ggrepel')
suppressPackageStartupMessages(library(ggrepel))