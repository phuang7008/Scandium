# this script is used to process the data from Summary File and draw histogram
# first remove everything
#
rm(list=c(ls()))
require(gridExtra)
require(ggplot2)

# passing/getting the command line arguments
#
#options(echo=TRUE)  # if you want to see commands in output file
##options()
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop("Please specify the unique pattern of the file you are processing, such as Capture.*Summary\n")
}

file_pattern <- args[1]

max_x_wgs <- 200
histogram_vals_wgs_df <- c()
PCT_of_Bases_Coverage_wgs_df <- c()
max_x_capture <- 200
histogram_vals_capture_df <- c()
PCT_of_Bases_Coverage_capture_df <- c()
binned_hist_df <- c()

process_histogram_data <- function(pattern_in, type_in) {
  tmp_file <- list.files(pattern=pattern_in, ignore.case = TRUE)
  
  # open file for reading
  #
  con <- file(tmp_file, "r")
  terminate_num <- 0
  pct_name<-c()     # Vector used to store pct name
  pct_vals<-c()     # Vector used to store pct vals
  i <- 1
  
  while (TRUE) {
    line <- readLines(con, n=1)
    if (length(line) <= 0) next
    
    # use average coverage to find the max x used on x-axis for drawing
    #
    if (grepl("average_coverage", line, ignore.case=TRUE)) {
      vals <- unlist(strsplit(line, '\t'))
      
      if (as.numeric(vals[2]) <= 30) {
        if (type_in == "wgs") max_x_wgs <<- 60
        if (type_in == "capture") max_x_capture <<- 60
      } else if (as.numeric(vals[2]) <= 55) {
        if (type_in == "wgs") max_x_wgs <<- 100
        if (type_in == "capture") max_x_capture <<- 100
      } else if (as.numeric(vals[2]) <= 80) {
        if (type_in == "wgs") max_x_wgs <<- 150
        if (type_in == "capture") max_x_capture <<- 150
      } else if (as.numeric(vals[2]) <= 120) {
        if (type_in == "wgs") max_x_wgs <<- 200
        if (type_in == "capture") max_x_capture <<- 200
      } else if (as.numeric(vals[2]) <= 160) {
        if (type_in == "wgs") max_x_wgs <<- 250
        if (type_in == "capture") max_x_capture <<- 250
      } else {
        if (type_in == "wgs") max_x_wgs <<- 300
        if (type_in == "capture") max_x_capture <<- 300
      }
    }
    
    # need to get the data frame for PCT_of_Bases_with coverage info for drawing
    #
    if (grepl("PCT_of_Bases_with", line, ignore.case=TRUE)) {
      if (grepl("6x", line, ignore.case=TRUE)) next
      if (grepl("11x", line, ignore.case=TRUE)) next
      
      pcts <- unlist(strsplit(line, '\t'))
      pcts[1] <- (gsub("\\PCT_of_Bases_with_", '', pcts[1]))
      pcts[1] <- (gsub("x\\_coverage", '', pcts[1]))
      pct_name[i] <- pcts[1]
      pct_vals[i] <- as.numeric(gsub("\\%", '', pcts[2]))
      i <- i+1
    }
    
    if (grepl(">>", line)) {
      # need to use this line to draw raw WGS histogram
      #
      tmp_data <- unlist(strsplit(line, ">>"))
      freqs <- unlist(strsplit(tmp_data[2], ','))
      Frequency<- as.numeric(freqs)
      Coverage <- seq(0, length(Frequency)-1, by=1)
      if (type_in == "wgs") histogram_vals_wgs_df <<- data.frame(Coverage, Frequency)
      if (type_in == "capture") histogram_vals_capture_df <<- data.frame(Coverage, Frequency)
      
      # now get the PCT_of_Bases_with Different number of coverage
      #
      if (type_in == "wgs") PCT_of_Bases_Coverage_wgs_df <<- data.frame(pct_name, pct_vals)
      if (type_in == "capture") {
        PCT_of_Bases_Coverage_capture_df <<- data.frame(pct_name, pct_vals)
        break
      }
    }
    
    if (grepl("==", line)) {
      # need to use this line to draw raw WGS histogram
      #
      tmp_data <- unlist(strsplit(line, "=="))
      freqs <- unlist(strsplit(tmp_data[2], ','))
      Frequency<- as.numeric(freqs)
      Coverage <- seq(0, length(Frequency)-1, by=1)
      binned_hist_df <<- data.frame(Coverage, Frequency)
      
      break
    }
  }
  
  if (type_in == "wgs") colnames(PCT_of_Bases_Coverage_wgs_df) <<- c("Coverage", "Percentage")
  if (type_in == "capture") colnames(PCT_of_Bases_Coverage_capture_df) <<- c("Coverage", "Percentage")
  
  close(con)
}

# now for Capture Annotations
#
exon_coverage <- c()
gene_coverage <- c()
transcript_coverage <- c()

process_capture_data <- function(pattern_in, type_in) {
  tmp_file <- list.files(pattern=pattern_in, ignore.case = TRUE)
  
  # open file for reading
  #
  con <- file(tmp_file, "r")
  i <- 1
  
  while (TRUE) {
    line <- readLines(con, n=1)
    if (length(line) < 1) break
    if (grepl("#", line, ignore.case = TRUE)) next
    if (grepl("_alt", line, ignore.case = TRUE)) next
    if (grepl("_random", line, ignore.case = TRUE)) next
    
    # need to get the data for gene coverage
    #
    vals <- unlist(strsplit(line, '\t'))
    if (type_in == "exon") exon_coverage[i] <<- as.numeric(gsub("\\%", '', vals[7]))
    if (type_in == "gene") gene_coverage[i] <<- as.numeric(vals[3])
    if (type_in == "transcript") transcript_coverage[i] <<- as.numeric(gsub("\\%", '',vals[6]))
    
    i <- i+1
  }
  
  close(con)
}


#########################################################################################
# Invoke the calls
#########################################################################################

# First for raw coverage count distribution
#
process_histogram_data(paste(file_pattern, ".*wgs.*summary", sep=""), "wgs")
process_histogram_data(paste(file_pattern, ".*capture.*summary", sep=""), "capture")

# I got the following error message if I don't do as the web suggested
# https://stackoverflow.com/questions/29278153/plotting-with-ggplot2-error-discrete-value-supplied-to-continuous-scale-on-c
# Error: Discrete value supplied to continuous scale
#
#PCT_of_Bases_Coverage_wgs_df$Coverage=as.numeric(levels(PCT_of_Bases_Coverage_wgs_df$Coverage))[PCT_of_Bases_Coverage_wgs_df$Coverage]
#PCT_of_Bases_Coverage_capture_df$Coverage=as.numeric(levels(PCT_of_Bases_Coverage_capture_df$Coverage))[PCT_of_Bases_Coverage_capture_df$Coverage]

png(filename = paste(file_pattern, "_coverage_graphs.png", sep=""), width=6, height=5, units="in", res=200)

wgs_hist <- ggplot(histogram_vals_wgs_df, aes(x=Coverage, y=Frequency)) +
  geom_point(color="green2", size=0.3) + xlim(0, max_x_wgs) +labs(title="Raw WGS Coverage Frequency Distribution") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# need to sort the factor level by the numeric number
#
PCT_of_Bases_Coverage_wgs_df$Coverage <- factor(PCT_of_Bases_Coverage_wgs_df$Coverage, 
                                                levels=sort(as.numeric(levels(PCT_of_Bases_Coverage_wgs_df$Coverage)))) 

wgs_pct <- ggplot(PCT_of_Bases_Coverage_wgs_df, aes(x=Coverage, y=Percentage)) + 
  geom_col(color="lightsalmon1", fill="lightsalmon1") + labs(title="WGS Percentage Base Coverage") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size=4),
    axis.text.y = element_text(size=4)
  )

capture_hist <- ggplot(histogram_vals_capture_df, aes(x=Coverage, y=Frequency)) + 
  geom_point(color="tomato1", size=0.3) + xlim(0,max_x_capture) + labs(title="Raw Capture Coverage Frequency Distribution") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

PCT_of_Bases_Coverage_capture_df$Coverage <- factor(PCT_of_Bases_Coverage_capture_df$Coverage, 
                                                    levels=sort(as.numeric(levels(PCT_of_Bases_Coverage_capture_df$Coverage))))

capture_pct <- ggplot(PCT_of_Bases_Coverage_capture_df, aes(x=Coverage, y=Percentage)) + 
  geom_col(color="slategray4", fill="slategray4") + labs(title="Capture Percentage Base Coverage") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size=4),
    axis.text.y = element_text(size=4)
  )

binned_graph <- ggplot(binned_hist_df, aes(x=Coverage, y=Frequency)) + 
  geom_point(color="green4", size=0.3) + xlim(0,max_x_capture) + labs(title="Smoothed WGS Coverage Frequency Distribution") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# Next for Capture Annotation Graphing
#
process_capture_data(paste(file_pattern, ".*Capture_Exon_pct.txt", sep=""), "exon")
process_capture_data(paste(file_pattern, ".*Capture_Gene_pct", sep=""), "gene")
process_capture_data(paste(file_pattern, ".*Capture_Transcript_pct",sep=""), "transcript")

# For Exon
#
#tmp_exon_hist <- hist(exon_coverage, plot=FALSE)
#highestExonCount <- max(tmp_exon_hist$counts)
exon_coverage_df <- data.frame(exon_coverage)

exon_hist <- ggplot(exon_coverage_df, aes(x=exon_coverage)) + geom_histogram(color="cyan4", fill="cyan4") +
  labs(title="Exon Coverage Histogram", x="Exon Coverage", y="Frequency") +
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# For Transcript
#
#tmp_trabscruot_hist <- hist(transcript_coverage, plot=FALSE)
#highestTranscriptCount <- max(tmp_transcript_hist$counts)
transcript_coverage_df <- data.frame(transcript_coverage)

transcript_hist <- ggplot(transcript_coverage_df, aes(x=transcript_coverage)) + 
  geom_histogram(color="cornsilk4", fill="cornsilk4") +
  labs(title="Transcript Coverage Histogram", x="Transcript Coverage", y="Frequency") +
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# For Gene
#
gene_coverage_df <- data.frame(gene_coverage)

gene_hist <- ggplot(gene_coverage_df, aes(x=gene_coverage)) + 
  geom_histogram(color="goldenrod2", fill="goldenrod2") + 
  labs(title="Gene Coverage Histogram", x="Gene Coverage", y="Frequency") +
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )
###############################
# Combine everything together
#
grid.arrange(wgs_pct, capture_pct, wgs_hist, capture_hist, binned_graph, exon_hist, gene_hist, transcript_hist, ncol=2)

dev.off()

