# this script is ussed to plot the distribution (barplot) for genome coverage data

fetch_whole_chromosome_data <- function(cov_file_name, output_file_name, chrom_id, chrom_len) {
  # creating a matrix with 1 columns and chrom_len rows
  mt_ <- matrix(0L, ncol=1, nrow=chrom_len)
  str(mt_)
  
  # read from file and process it accordingly to retrieve the coverage information
  con <- file(cov_file_name, "r")
  items = c()
  start = 0
  end = 0

  while (TRUE) {
    line = readLines(con, n=1)
    if (length(line) == 0)
      break
    
    items <- unlist(strsplit(line, '\t'))

    if (items[1] == chrom_id) {
      #print(items[1])
      start = as.numeric(items[2])
      len   = as.numeric(items[3]) - start		# for bed file, the end/stop position is not counted!

      if (length(items) == 5) {
        cov = as.numeric(items[5])
      } else {
        cov = as.numeric(items[4])
      }

      for (idx in c(1:len)) {
        pos = idx + start
        mt_[pos,1] <- as.integer(cov + 0.5)
        #cat(pos, cov, "\n")
      }
    } 
  }
  close(con)
  
  # now write data into a bed file
  #
  fileConn <-file(output_file_name)

  prev_val <- 0
  start_pos <- 1
  stop_pos <- 2

  for (i in 2:nrow(mt_)) {
    if (mt_[i,1] == prev_val) {
      stop <- stop+1
    } else {
      # write to the file
      writeLines(c(chrom_id, start, stop, prev_val), fileConn)

      # now update values
      start <- i
      stop <- i+1
      prev_val <- mt_[i,1]
    }
  }
  close (fileConn)

  print("Finishing dumping to data frame")
  rm(mt_)
  rm(df_)
}

# get the command line arguments
options(echo=TRUE)	# if you want to see commands in output file
#options()
args <- commandArgs(trailingOnly=TRUE)
print(args)
cov_file_in <- args[1]
output_file_in <- args[2]
chrom_id_in <- args[3]
version <- args[4]

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
	v <- c( "1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663,
			"8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540,
			"15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520, "21"=48129895,
			"22"=51304566, "X"=155270560, "Y"=59534049, "MT"=16569)
}

for (id in id_names) {
  if (chrom_id_in == id) {
	fetch_whole_chromosome_data(cov_file_in, output_file_in, id, v[id])
	gc()
	break;
  }
}
