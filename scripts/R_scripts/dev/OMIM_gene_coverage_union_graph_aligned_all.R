# this script is used to plot the distribution (barplot) for the coverage 
# of sequencing data that related to the OMIM genes 

coverage_plot <- function(fName, chrom_id, chrom_len, low, high, type) {
  # creating a matrix with 1 columns and chrom_len rows
  #
  mt_ <- matrix(0L, ncol=1, nrow=chrom_len)
  str(mt_)
  
  # load the individual library low coverage regions in
  #
  seq_data <- read.table(fName)

  start <- 1 

  # loop through input file
  #
  for (i in 1:nrow(seq_data)) {
    seq_row <- seq_data[i,]
    seq_start <- seq_row[1,2]
    seq_stop  <- seq_row[1,3]
    len_ <- seq_row[1,3] - seq_row[1,2]

	for (j in 1:len_) {
	  mt_[start,1] <- as.integer(seq_row[1,4] + 0.5)
	  start <- start + 1
    }
  }

  print(start)
  
  df_ <- as.data.frame(mt_)
  print("Finishing convert to data frame")
  
  library(lattice, pos=10) 
  
  # get the basename of the file without extension
  #
  f_name <-basename(fName)
  f_name_splitted <- unlist(strsplit(f_name, '\\.'))
  
  tmp_file <- paste("plot_for_", f_name_splitted[1], sep = "")
  tmp_file <- paste(tmp_file, "_union_chr", sep = "")
  tmp_file <- paste(tmp_file, chrom_id, sep = "")
  tmp_file1 <- tmp_file
  tmp_file2 <- tmp_file
  tmp_file3 <- tmp_file

  # needed for three different image files
  #
  tmp_file1 <- paste(tmp_file1, "_300x_", sep = "")
  tmp_file2 <- paste(tmp_file2, "_100x_", sep = "")
  tmp_file3 <- paste(tmp_file3, "_20x_",  sep = "")

  if (type == "1") {
    tmp_file1 <- paste(tmp_file1, "coverage.png", sep = "")
    tmp_file2 <- paste(tmp_file2, "coverage.png", sep = "")
    tmp_file3 <- paste(tmp_file3, "coverage.png", sep = "")
  } else {
    tmp_file1 <- paste(tmp_file1, "coverage_log.png", sep = "")
    tmp_file2 <- paste(tmp_file2, "coverage_log.png", sep = "")
    tmp_file3 <- paste(tmp_file3, "coverage_log.png", sep = "")
  }
  
  title_ <- paste("OMIM_Gene_Coverage For Chromosome ", chrom_id, sep="")

  if (type == "1") {
    cov_plot1 <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(low,300), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)

    cov_plot2 <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",  
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(low,100), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)

    cov_plot3 <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",  
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(low,20), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)

  } else {
    cov_plot1 <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(1,6), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)

    cov_plot2 <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",  
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(1,5), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)

    cov_plot3 <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",  
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(1,3), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_) 
  }
  
  png(filename = tmp_file1, width=15, height=4, units="in", res=200)
  plot(cov_plot1)

  png(filename = tmp_file2, width=15, height=4, units="in", res=200)                                          
  plot(cov_plot2)

  png(filename = tmp_file3, width=15, height=4, units="in", res=200)                                          
  plot(cov_plot3)

  dev.off()

  print("for loess(df_[,1] ~ df_[,2], df_, span=0.5)")
  loess(df_[,1] ~ df_[,2], df_, span=0.5)

  print("for loess(df_[,1] ~ df_[,2], df_, span=0.75)")
  loess(df_[,1] ~ df_[,2], df_, span=0.75)
  
  print("for loess(df_[,1] ~ df_[,2], df_, span=1)")
  loess(df_[,1] ~ df_[,2], df_, span=1)

  print("for loess(df_[,1] ~ df_[,2], df_) default span=")
  loess(df_[,1] ~ df_[,2], df_)

  rm(mt_)
  rm(df_)
}

# get the command line arguments
options(echo=TRUE)	# if you want to see commands in output file
#options()
args <- commandArgs(trailingOnly=TRUE)
print(args)
file_in <- args[1]
chrom_id_in <- args[2]
low <- as.numeric(args[3])
high <- as.numeric(args[4])
g_length <- as.numeric(args[5])

# type 1 is absolute value, while type 2 is log value
type <- args[6]

coverage_plot(file_in, chrom_id_in, g_length+10, low, high, type)
