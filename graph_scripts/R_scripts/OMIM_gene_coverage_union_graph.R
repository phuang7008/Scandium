# this script is used to plot the distribution (barplot) for the coverage 
# of sequencing data that related to the OMIM genes 

coverage_plot <- function(combinedFName, fName, chrom_id, chrom_len, low, high, type) {
  # creating a matrix with 1 columns and chrom_len rows
  #
  mt_ <- matrix(0L, ncol=1, nrow=chrom_len)
  str(mt_)
  
  # load the entire combined sequence in
  #
  seq_data <- read.table(combinedFName)
  #str(seq_data)

  # load the individual library low coverage regions in
  #
  updated_data <- read.table(fName)

  start <- 1 

  # loop through combined sequence data
  #
  for (i in 1:nrow(seq_data)) {
    seq_row <- seq_data[i,]
    seq_start <- seq_row[1,2]
    seq_stop  <- seq_row[1,3]
    len_ <- seq_row[1,3] - seq_row[1,2] + 1

    # loop through current library sequence data
    # if inside, set the hit to 1
    #
    hit <- 0
    prev_start <- seq_start
    prev_stop  <- seq_start
    
    # make a deep copy here so that we can still remove rows from updated_data, while keep my_data unchanged
    # until we do the re-assignment here
    #
    my_data <- data.frame(updated_data)
    
    for (j in 1:nrow(my_data)) {
      tmp_row <- my_data[j,]
      tmp_start <- tmp_row[1,2]
      tmp_stop  <- tmp_row[1,3]

      # if it is outside the range like the case below, skip it
      #                     seq_start ========================== seq_stop
      # tmp_start -------------- tmp_stop
      #
      if (seq_start > tmp_stop) {
        next
      }
 
      # if it is outside the range like the case below, break out as we don't need to loop further
      #                     seq_start ========================== seq_stop
      #                                                    tmp_start ----------------- tmp_stop
      #
      if (seq_stop < tmp_start) {
        break
      }
      
      if ((seq_start <= tmp_start && tmp_start <= seq_stop) ||
          (seq_start <= tmp_stop  && tmp_stop  <= seq_stop)) {
        hit <- 1
      
        # Need to loop through the region to see where the current library hit
        #
        idx <- 1
        pre_distance <- 0		# this is the distance between tmp_start - seq_start
        while (idx <= len_) {
          # case by case
          #
          if ((prev_stop + pre_distance) < tmp_start && tmp_stop <= seq_stop) {
            # seq_start ========================================================================= seq_stop
            #    tmp_start1 ----------------- tmp_stop1    tmp_start2 ----------------- tmp_stop2
            #
            pre_distance <- tmp_start - prev_stop
            start <- start + pre_distance
            idx <- idx + pre_distance
      
          } else if (tmp_start == (prev_stop + pre_distance) && tmp_stop <= seq_stop) {				
            # seq_start  ========================================================================= seq_stop
            # tmp_start1 ----------------- tmp_stop1       tmp_start2 ----------------- tmp_stop2
            #
            tmp_idx <- 1
            tmp_len = tmp_stop - tmp_start
            while (tmp_idx <= tmp_len) {
              mt_[start,1] <- as.integer(tmp_row[1,4] + 0.5)

              tmp_idx <- tmp_idx + 1
              start <- start + 1 
              #idx <- idx + 1 
            }
            
            # taken care of the last section between tmp_stop towards seq_stop
            #start <- start + (seq_stop - tmp_stop)
            prev_start <- tmp_start
            prev_stop  <- tmp_stop
            #tmp_start = seq_stop + tmp_len
            #tmp_stop  = seq_stop + tmp_len

            # update the updated_data by removing the rows that has been processed
            #
            updated_data <- updated_data[-1,]
            
            # set the exit condition
            idx <- len_ + 1
            break
          } else {
            #    seq_start ============================== seq_stop
            # tmp_start ----------------- tmp_stop
            # this case is impossible!
            #
            print(c(start, seq_start, seq_stop, tmp_start, tmp_stop, prev_start, prev_stop, idx))
            print ("Something is wrong: Impossible cases!")
            break;
          }
        }
      }
    }
    
    # if this region doesn't hit anything, move the start position along the X-axis
    #
    if (hit == 0) {
      start <- start + len_
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
  tmp_file <- paste(tmp_file, "_", sep = "")
  tmp_file <- paste(tmp_file, high, sep = "")
  tmp_file <- paste(tmp_file, "x", sep = "")
  if (type == "1") {
    tmp_file <- paste(tmp_file, "_coverage.png", sep = "")
  } else {
    tmp_file <- paste(tmp_file, "_coverage_log.png", sep = "")
  }
  
  title_ <- paste("OMIM_Gene_Coverage For Chromosome ", chrom_id, sep="")

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
combined_file <- args[1]
file_in <- args[2]
chrom_id_in <- args[3]
low <- as.numeric(args[4])
high <- as.numeric(args[5])
g_length <- as.numeric(args[6])

# type 1 is absolute value, while type 2 is log value
type <- args[7]

coverage_plot(combined_file, file_in, chrom_id_in, g_length+10, low, high, type)

#coverage_plot("combined_low_cov_position_for_chrom_7_sorted_merged", "OMIM_only.BED", "7", 17102336, 1, 100, 1)

