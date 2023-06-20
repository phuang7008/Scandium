# this script is used to process the data from ExCID run
# the input data should be obtained from group_bins_together.pl
# first remove everything
# rm(list=c(ls()))
#
draw_distribution <- function(file_in, lower, upper) {
  # get the base name of the file without extension
  #
  f_name <-basename(file_in)
  f_name_splitted <- unlist(strsplit(f_name, '\\.'))
  tmp_file <- paste(f_name, "_Cov_Dist.png", sep = "")
  #tmp_file <- paste(f_name_splitted[1], "_Coverage_Distribution_Smoothed.png", sep = "")
  #tmp_file <- paste(f_name_splitted[1], "_Coverage_Distribution.jpg", sep = "")

  # read everything into a table from the input file 
  #
  df <- read.table(file_in, header=FALSE)

  #png(filename=tmp_file, width=5, height=3, units="in", res=200)
  png(filename=tmp_file)

  plot(df[lower:upper, 1], df[lower:upper, 2])
  dev.off()
}

# get the command line arguments
#
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file_in <- args[1]
lower_bound <- args[2]
upper_bound <- args[3]

draw_distribution(file_in, lower_bound, upper_bound)
