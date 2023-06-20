# this script is used to process the data from ExCID run
# the input file should be the WGS Summary Report file
#
draw_distribution <- function(file_in, lower, upper, version, type) {
  # get the base name of the file without extension
  #
  f_name <-basename(file_in)
  f_name_splitted <- unlist(strsplit(f_name, '\\.'))
  tmp_file <- paste(f_name_splitted[1], "_", sep = "")
  tmp_file <- paste(tmp_file, lower, sep = "")
  tmp_file <- paste(tmp_file, "_", sep = "")
  tmp_file <- paste(tmp_file, upper, sep = "")
  if (type == 2) {
    tmp_file <- paste(tmp_file, "_Coverage_Distribution_Density_Curve.png", sep = "")
  } else {
    tmp_file <- paste(tmp_file, "_Coverage_Distribution_Histogram.png", sep = "")
  }

  # read from file and process it 
  #
  con <- file(file_in, "r")

  flag <- 0
  ylab <- c()
  xdat <- c()
  num_of_Ns <- 0
  if (version == 'hg38') {
    num_of_Ns <- 173893331
  } else {
    num_of_Ns <- 237019493
  }

  while (TRUE) {
    line = readLines(con, n=1)    # read one line at a time
    #if (is.na(line))
	#  next
  
    if (flag == 3)
      break
  
    items <- unlist(strsplit(line, " "))
  
    if (flag == 2) {
      xdat <-as.numeric(unlist(strsplit(line, ",")))
      xdat[1] <- xdat[1] - num_of_Ns
      flag = 3

      # if it is a density curve, we need to do extra calculation
      # convert the data to percentage/density
      #
      if (type == 2) {
        total <- sum(xdat)
        xdat  <- format(round(xdat/total, 5), nsmall=5)
      }
   }
  
    if (flag == 1) {
      ylab = unlist(strsplit(line, ","))
      flag = 2
    }
  
	#print(items[2])
	if (is.na(items[2]))
	  next

    if (items[2] == "Histogram") {
      flag <- 1
    }
  }

  close(con)

  png(filename=tmp_file)

  if (type == 2) {
    plot( ylab[c(lower:upper)], xdat[c(lower:upper)]) 
          #xaxt="n", xlim=c(lower, upper))
          #xaxt="n", xlim=c(lower, upper), yaxt="n", ylim=c(0, 1))
    #axis(1, at=lower:upper)
    #axis(2, at=0:1)
  } else {
    plot(ylab[c(lower:upper)], xdat[c(lower:upper)])
  }

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
version <- args[4]
type <- args[5]		# type is either histogram (1) or density curve (2)

draw_distribution(file_in, lower_bound, upper_bound, version, type)
