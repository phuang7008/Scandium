# this script is ussed to plot the distribution (barplot) for genome coverage data

coverage_plot <- function(fName, chrom_id, chrom_len, low, high) {
  # creating a matrix with 1 columns and chrom_len rows
  mt_ <- matrix(0L, ncol=1, nrow=chrom_len)
  str(mt_)
  
  # read from file and process it accordingly to retrieve the coverage information
  con <- file(fName, "r")
  items = c()
  start = 0
  end = 0

  # get the basename of the file without extension
  f_name <-basename(fName)
  f_name_splitted <- unlist(strsplit(f_name, '\\.'))

  tmp_file <- paste("plot_for_", f_name_splitted[1], sep = "")
  tmp_file <- paste(tmp_file, "_chr", sep = "")
  tmp_file <- paste(tmp_file, chrom_id, sep = "")
  tmp_file <- paste(tmp_file, "_log_coverage.png", sep = "")
  title_ <- paste(low, "x-", high, "x_Coverage For Chromosome ", chrom_id, sep="")
  
  while (TRUE) {
    line = readLines(con, n=1)
    if (length(line) == 0)
      break
    
    items <- unlist(strsplit(line, '\t'))

    if (items[1] == chrom_id) {
      #print(items[1])
      start = as.numeric(items[2])
      len   = as.numeric(items[4])
      cov   = log(as.numeric(items[5]))
      #cov   = log10(as.numeric(items[5]))

      for (idx in c(0:len-1)) {
        pos = idx + start
        mt_[pos,1] <- cov
        #print(cat(pos, cov, "\n"))
      }
    } 
  }
  close(con)
  
  df_ <- as.data.frame(mt_)
  print("Finishing convert to data frame")
  
  library(lattice, pos=10) 

  cov_plot <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.3, fill="darkgreen",
                        #auto.key=list(border=TRUE, cex=0.3), par.settings=simpleTheme(pch=16,cex=0.3), 
						xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(1,4), na.rm = TRUE),
						#xlab = "Chromosome Position", ylab="Coverage/Position", ylim=range(c(low,as.integer(log(high))+1), na.rm = TRUE),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)
  
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
version <- args[5]

# define chromosome id name vector
id_names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
              "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")

# define lengths of each chromosome for version HG19 (v37) and v38
v <- c()

if (version == "hg38") {
	v <- c( "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555, "5"=181538259, "6"=170805979, "7"=159345973, 
			"8"=145138636, "9"=138394717, "10"=133797422, "11"=135086622, "12"=13327530, "13"=114364328, "14"=107043718, 
			"15"=101991189, "16"=90338345, "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167, "21"=46709983, 
			"22"=50818468, "X"=156040895, "Y"=57227415, "MT"=16569)
} else {
	v <- c( "1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, 
			"5"=180915260, 
			#"5"=66513080, 
			"6"=171115067, "7"=159138663,
			"8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540,
			"15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520, "21"=48129895,
			"22"=51304566, "X"=155270560, "Y"=59534049, "MT"=16569)
}

for (id in id_names) {
  if (chrom_id_in == id) {
	coverage_plot(file_in, id, v[id] + 1, low, high)
	gc()
	break;
  }
}
