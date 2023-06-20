# this script is ussed to plot the distribution (barplot) for genome coverage data

coverage_plot <- function(fName, chrom_id, chrom_len) {
  # create a data frame with position 'POS' initialized and everything else is 0
  #df1 <- data.frame('CID'=numeric(0), 'POS'=numeric(0), 'COV'=numeric(0))
  df2 <- data.frame(matrix(0L, ncol=3, nrow=chrom_len))
  colnames(df2) <- c('CID', 'POS', 'COV')
  df2 <- transform(df2, CID=as.numeric(CID), POS=as.numeric(POS), COV=as.numeric(COV))
  #df2[is.na(df2)] <- 0
  
  #for (pos in c(1:chrom_len))
  #  df2[pos, 'POS'] <- pos
  
  str(df2)
  
  # need to convert the coverage.fasta file format into a data frame
  con <- file(fName, "r")
  items = c()
  start = 0
  end = 0
  tmp_file <- paste("plot_for_", chrom_id, sep = "")
  tmp_file <- paste(tmp_file, "_coverage.png", sep = "")
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
        start = 0
        end = 0
        print("No match!\n")
        break
      }
    } else {
      items <- strsplit(line, ' ')
      items <- as.numeric(unlist(items))
      print(length(items))
      
      for (idx in c(1:length(items))) {
        pos = idx + start
        df2[pos, 'COV'] <- items[idx]
        #        print(df2[pos, 'COV'])
        
        #        start = start + 1
        #        if (start == end)
        #          break
      }
    }
  }
  close(con)
  
  library(lattice, pos=10) 
  #cov_plot <- xyplot(df2[,'COV'] ~ df2[,'POS'], type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), 
  #                   scales=list(x=list(relation='same'), y=list(relation='same')), data=df2, main="Coverage for Chromosome 1")
  cov_plot <- xyplot(df2[,'COV'] ~ 1:nrow(df2), type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), 
                     scales=list(x=list(relation='same'), y=list(relation='same')), data=df2, main="Coverage for Chromosome 1")
  
  png(filename = tmp_file, width=10, height=4, units="in", res=1000)
  plot(cov_plot)
  dev.off()
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
  #coverage_plot('/stornext/snfs5/next-gen/scratch/phuang/dev/make_graphs/cov2000.fasta', id, v37[id])
  #coverage_plot('/Users/phuang/Downloads/HHYFVALXX-6.hgv.bam.cov.fasta', id, v37[id])
  break
}

