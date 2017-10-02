# this script is ussed to plot the distribution (barplot) for genome coverage data
#rm(list=c(ls()))

coverage_plot <- function(fName, chrom_id, chrom_len) {
  # creating a matrix with 1 columns and chrom_len rows
  mt_ <- matrix(0L, ncol=1, nrow=chrom_len)
  str(mt_)
  
  # read from file and process it accordingly to retrieve the coverage information
  con <- file(fName, "r")
  items = c()
  start = 0
  end = 0
  tmp_file <- paste("plot_for_", chrom_id, sep = "")
  tmp_file <- paste(tmp_file, "_coverage.png", sep = "")
  title_ <- paste("Coverage For Chromosome ", chrom_id, sep="")
  chrom_id = paste(">", chrom_id, sep="")
  
  while (TRUE) {
    line = readLines(con, n=1)
    if (length(line) == 0)
      break
    
    if (grepl('>', line)) {
      x <- strsplit(line, ' ')
      print(x[[1]][1])
      if (x[[1]][1] == chrom_id) {
        start <-as.numeric( x[[1]][2])
        end <- as.numeric(x[[1]][3])
        print(start)
      } else {
		if (start > 0) {
          start = 0
          end = 0
          break
        }
		print("No match!\n")
      }
    } else {
      items <- strsplit(line, ' ')
      items <- as.numeric(unlist(items))
      print(length(items))

      if (length(items) > 0 && start > 0) {
		for (idx in c(1:length(items))) {
		  pos = idx + start
	      mt_[pos,1] <- items[idx]
		}
      }
    }
  }
  close(con)
  
  df_ <- as.data.frame(mt_)
  print("Finishing convert to data frame")
  
  library(ggplot2) 
  ggplot(df_, aes(1:nrow(df_), df_[,1])) + geom_point(color="darkgreen", size=0.1) + 
           xlab("Chromosome Position") + ylab("Coverage/Position")
  ggsave(tmp_file, width=10, height=3, dpi=150)

  # cleanup the variables
  rm(mt_)
  rm(df_)
}

# define chromosome id name vector
id_names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
              "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")

# define lengths of each chromosome for version HG19 (v37) and v38
v37 <- c("1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663,
         "8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540,
         "15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520, "21"=48129895,
         "22"=51304566, "X"=155270560, "Y"=59034049, "MT"=16569)

v38 <- c()

for (id in id_names) {
  coverage_plot('/stornext/snfs5/next-gen/scratch/phuang/dev/make_graphs/HHYFVALXX-6.hgv.bam.cov.fasta', id, v37[id])
  gc()		# garbage collection
  #coverage_plot('/Users/phuang/Downloads/cov20.fasta', id, 900000)
  #break
}


