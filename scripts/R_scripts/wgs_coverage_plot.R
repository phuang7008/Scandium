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
      if (x[[1]][1] == chrom_id) {
        print(x[[1]][1])
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
      #print(length(items))

      if (length(items) > 0) {
		for (idx in c(0:(length(items)-1))) {
		  pos = idx + start
	      mt_[pos,1] <- items[idx+1]
		}
        start = start + length(items)
      }
    }
  }
  close(con)
  
  df_ <- as.data.frame(mt_)
  print("Finishing convert to data frame")
  
  library(lattice, pos=10) 

  cov_plot <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch=16, col="darkgreen", cex=0.2, fill="darkgreen",
                        auto.key=list(border=TRUE, cex=0.2), par.settings=simpleTheme(pch=16,cex=0.2),
                        scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)

  #cov_plot <- xyplot(df_[,1] ~ 1:nrow(df_), type="p", pch='.', auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), 
  #                    scales=list(x=list(relation='same'), y=list(relation='same')), data=df_, main=title_)
  
  png(filename = tmp_file, width=12, height=4, units="in", res=200)
  plot(cov_plot)
  dev.off()

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
         "22"=51304566, "X"=155270560, "Y"=59534049, "MT"=16569)

v38 <- c()

for (id in id_names) {
  coverage_plot('/stornext/snfs5/next-gen/scratch/phuang/dev/make_graphs/HHYFVALXX-6.hgv.bam.wgs.fasta', id, v37[id] + 1)
  gc()
}
