# this script is ussed to plot the distribution (barplot) for genome coverage data

coverage_plot <- function(fName, chrom_id, chrom_len, low, high, type) {
  # creating a matrix with 1 columns and chrom_len rows
  mt_ <- matrix(0L, ncol=1, nrow=chrom_len)
  str(mt_)
  
  # read from file and process it accordingly to retrieve the coverage information
  con <- file(fName, "r")
  items <- c()
  start <- 1 

  # get the basename of the file without extension
  f_name <-basename(fName)
  f_name_splitted <- unlist(strsplit(f_name, '\\.'))

  tmp_file <- paste("plot_for_", f_name_splitted[1], sep = "")
  tmp_file <- paste(tmp_file, "_chr", sep = "")
  tmp_file <- paste(tmp_file, chrom_id, sep = "")
  if (type == "1") {
    tmp_file <- paste(tmp_file, "_coverage.png", sep = "")
  } else {
    tmp_file <- paste(tmp_file, "_coverage_log.png", sep = "")
  }

  title_ <- paste("OMIM_Gene_Coverage For Chromosome ", chrom_id, sep="")

  while (TRUE) {
    line = readLines(con, n=1)
    if (length(line) == 0)
      break
    
    items <- unlist(strsplit(line, '\t'))

    if (items[1] == chrom_id) {
      #print(items[1])
      #start = as.numeric(items[2])
      len <- as.numeric(items[4])
      if (type == "1") {
        cov <- as.numeric(items[5])
      } else {
        cov <- log(as.numeric(items[5]))
      }
      for (idx in c(1:len)) {
        mt_[start,1] <- cov
        #cat(f_name, "\t", start, "\t", cov, "\n")
        start <- 1 + start
      }
    } 
  }
  close(con)
  
  df_ <- as.data.frame(mt_)
  print("Finishing convert to data frame")
  
  library(lattice, pos=10) 

  if (type == "1") {
     cov_plot <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",
                        #auto.key=list(border=TRUE, cex=0.3), par.settings=simpleTheme(pch=16,cex=0.3), 
						xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(low,high), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)
  } else {
     cov_plot <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",
                        #auto.key=list(border=TRUE, cex=0.3), par.settings=simpleTheme(pch=16,cex=0.3), 
                        xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(1,5), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)
  }
  
  png(filename = tmp_file, width=15, height=4, units="in", res=200)
  plot(cov_plot)
  dev.off()

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
