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
  ydat <- c()
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
      ydat <-as.numeric(unlist(strsplit(line, ",")))
      ydat[1] <- ydat[1] - num_of_Ns
      flag = 3

      # if it is a density curve, we need to do extra calculation
      # convert the data to percentage/density
      #
      if (type == 2) {
        total <- sum(ydat)
        ydat  <- as.numeric(format(round(ydat/total, 5), nsmall=5))
      }
   }
  
    if (flag == 1) {
      # need to convert data to numeric as we will use min() and max()
      # to get the ticker breaks
      #
      xdat = as.numeric(unlist(strsplit(line, ",")))
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

  #print(xdat)
  #print(ydat)
  print(min(ydat))
  print(max(ydat))
  library(ggplot2)
  library(dplyr)

  # here we need to convert two vectors into a data.frame
  #
  df <- data.frame(cnt=ydat, cov=xdat)
  g <- 0

  num_of_X_breaks <- (upper - lower)/10
  num_of_Y_breaks <- (max(ydat) - min(ydat))/5

  if (type == 2) {
    #g <- ggplot(df, aes(cov, cnt)) + geom_bar(stat="identity")
    g <- ggplot(df, aes(cov, cnt)) + geom_point() + xlab("Coverage") + ylab("Percentage Of Total Bases")
    num_of_Y_breaks <- round(num_of_Y_breaks * 100)/100
    if (num_of_Y_breaks == 0)
      num_of_Y_breaks <- 0.01
  } else {
    #g <- ggplot(df, aes(df$cov, df$cnt)) + geom_density(data=df, adjust = 0.5) + xlab("Coverage") + ylab("Total Base Count") 
    g <- ggplot(df, aes(cov, cnt)) + geom_point() + xlab("Coverage") + ylab("Total Base Count") 
  }

  g <- g + theme_bw(base_size = 14, base_family = "")
  g <- g + scale_x_continuous(breaks=round(seq(min(df$cov), max(df$cov), by=num_of_X_breaks),1), limits= c(lower, upper))
  g <- g + scale_y_continuous(breaks = seq(min(df$cnt), max(df$cnt), by = num_of_Y_breaks))
  #g <- g + scale_y_continuous(breaks = round(seq(min(df$cnt), max(df$cnt), by = num_of_Y_breaks),1), limits=c(max(df$cnt), min(df$cnt)))

  #print(g)
  ggsave(tmp_file, dpi=100, width=20, height=10, units="cm")

  #dev.off()
}

# get the command line arguments
#
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
file_in <- args[1]
lower_bound <- as.numeric(args[2])
upper_bound <- as.numeric(args[3])
version <- args[4]
type <- args[5]		# type is either histogram (1) or density curve (2)

draw_distribution(file_in, lower_bound, upper_bound, version, type)
