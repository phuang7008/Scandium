# this script is used to process the data from Summary File and draw histogram
# first remove everything
#
rm(list=c(ls()))
require(gridExtra)
require(ggplot2)
require(scales)

# passing/getting the command line arguments
#
#options(echo=TRUE)  # if you want to see commands in output file
##options()
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("\nPlease specify the prefix of the file you are processing, such as H72CKCCXX-1\n\n")
  stop("\nPlease specify the capture file prefix you are processing, such as Parkinson_Genes_hg38\n\n")
  stop("\nPlease specify the minimal base threshold, such as 20\n\n")
}

file_pattern <- args[1]
capture_pattern <- args[2]
min_base_cutoff <- args[3]

max_x_wgs <- 200
histogram_vals_wgs_df <- c()
PCT_of_Bases_Coverage_wgs_df <- c()
max_x_capture <- 200
histogram_vals_capture_df <- c()
PCT_of_Bases_Coverage_capture_df <- c()
binned_hist_df <- c()
categories <- c("<90%", "90-97%", "97-99%", "100%")

process_histogram_data <- function(pattern_in, type_in) {
  tmp_file <- list.files(pattern=glob2rx(pattern_in), ignore.case = TRUE)
  print(tmp_file)
  
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
      Frequency <- Frequency * 100 / sum(Frequency)
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
      Frequency <- Frequency * 100 / sum(Frequency)
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
exon_coverage <- c(0,0,0,0)
gene_coverage <- c(0,0,0,0)
transcript_coverage <- c(0,0,0,0)

process_capture_data <- function(pattern_in, type_in) {
  #print(pattern_in)
  #pattern_tmp = paste(pattern_in, pattern_in_cap, sep='*')
  tmp_file <- list.files(pattern=glob2rx(pattern_in), ignore.case = TRUE)
  print(tmp_file)
  
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
    if (type_in == "exon") 
      grouping_data(type_in, as.numeric(gsub("\\%", '', vals[7])))

    if (type_in == "gene") 
      grouping_data(type_in, as.numeric(vals[3]))

    if (type_in == "transcript")
      grouping_data(type_in, as.numeric(gsub("\\%", '',vals[6])))

  }
  
  close(con)
}

grouping_data <- function(type_in, val_in) {
  if (val_in == 100) {
    if (type_in == "exon") exon_coverage[4] <<- exon_coverage[4] + 1
    if (type_in == "gene") gene_coverage[4] <<- gene_coverage[4] + 1
    if (type_in == "transcript") transcript_coverage[4] <<- transcript_coverage[4] + 1
  } else if (val_in >= 97) {
    if (type_in == "exon") exon_coverage[3] <<- exon_coverage[3] + 1
    if (type_in == "gene") gene_coverage[3] <<- gene_coverage[3] + 1
    if (type_in == "transcript") transcript_coverage[3] <<- transcript_coverage[3] + 1
  } else if (val_in >= 90) {
    if (type_in == "exon") exon_coverage[2] <<- exon_coverage[2] + 1
    if (type_in == "gene") gene_coverage[2] <<- gene_coverage[2] + 1
    if (type_in == "transcript") transcript_coverage[2] <<- transcript_coverage[2] + 1
  } else {
    if (type_in == "exon") exon_coverage[1] <<- exon_coverage[1] + 1
    if (type_in == "gene") gene_coverage[1] <<- gene_coverage[1] + 1
    if (type_in == "transcript") transcript_coverage[1] <<- transcript_coverage[1] + 1
  }
}

#########################################################################################
# Invoke the calls
#########################################################################################

# First for raw coverage count distribution
#
process_histogram_data(paste(file_pattern, "*WGS*Summary*", sep=""), "wgs")
process_histogram_data(paste(paste(file_pattern, capture_pattern, sep="*"), "*Capture*Summary*", sep=""), "capture")

tmp_pattern = paste(file_pattern, capture_pattern, sep="_")
png(filename = paste(tmp_pattern, "_coverage_graphs.png", sep=""), width=6, height=5, units="in", res=200)
#pdf(paste(file_pattern, "_coverage_graphs.pdf", sep=""), width=7, height=5)

wgs_hist <- ggplot(histogram_vals_wgs_df, aes(x=Coverage, y=Frequency)) +
  geom_point(color="green2", size=0.3) + xlim(0, max_x_wgs) +labs(title="C: Raw WGS Coverage Frequency Distribution") + ylab("Frequency (%)") +
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
  geom_col(color="lightsalmon1", fill="lightsalmon1") + labs(title="A: WGS Percentage Base Coverage") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size=4),
    axis.text.y = element_text(size=4)
  )

capture_hist <- ggplot(histogram_vals_capture_df, aes(x=Coverage, y=Frequency)) + 
  geom_point(color="tomato1", size=0.3) + xlim(0, max_x_capture) + labs(title="D: Raw Capture Coverage Frequency Distribution") + ylab("Frequency (%)") +
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
  geom_col(color="slategray4", fill="slategray4") + labs(title="B: Capture Percentage Base Coverage") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size=4),
    axis.text.y = element_text(size=4)
  )

binned_graph <- ggplot(binned_hist_df, aes(x=Coverage, y=Frequency)) + 
  geom_point(color="green4", size=0.3) + xlim(0, max_x_capture) + labs(title="E: Smoothed WGS Coverage Frequency Distribution") + ylab("Frequency (%)") + 
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# Next for Capture Annotation Graphing
#
tmp_pattern = paste(file_pattern, capture_pattern, sep="*")
process_capture_data(paste(tmp_pattern, "*Capture_CDS_pct*", sep=""), "exon")
process_capture_data(paste(tmp_pattern, "*Capture_Gene_pct*", sep=""), "gene")
process_capture_data(paste(tmp_pattern, "*Capture_Transcript_pct*",sep=""), "transcript")

# For Exon
#
num_of_exons <- sum(exon_coverage)
ec <- exon_coverage
exon_coverage <- round(exon_coverage * 100 / sum(exon_coverage), 2)
exon_coverage_df <- data.frame(exon_coverage)
exon_coverage_df <- cbind(categories, exon_coverage_df)
ec_text <- c()
for (id in 1:4) {
	ec_text[id] = paste(exon_coverage[id], "%\n", sep="")
	ec_text[id] = paste(ec_text[id], "(", sep="")
	ec_text[id] = paste(ec_text[id], ec[id], sep="")
	ec_text[id] = paste(ec_text[id], ")", sep="")
}

# lock in factor level order
exon_coverage_df$categories <- factor(exon_coverage_df$categories, levels = exon_coverage_df$categories)

e_title = paste("F: Exon Coverage Stats (total # of exons: ", num_of_exons)
e_title = paste(e_title, ")")

print(exon_coverage[4])
y_limit = 125
if (exon_coverage[4] > 80.0) y_limit = 135
if (exon_coverage[4] > 90.0) y_limit = 150
if (exon_coverage[4] < 60.0) y_limit = 110
#print(y_limit)

x_label <- paste("Base Coverage Percentages of Exons (min base coverage >= ", min_base_cutoff, sep="")
x_label <- paste(x_label, ")", sep="")

exon_hist <- ggplot(exon_coverage_df, aes(y=exon_coverage, x=categories)) + 
  geom_col(color="cyan4", fill="cyan4", width=0.5) + geom_text(aes(label=ec_text), size=1.6, vjust=-0.5) +
  labs(title=e_title, x=x_label, y="Coverage Stats") +  ylim(0, y_limit) +
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# For Transcript
#
num_of_transcripts <- sum(transcript_coverage)
tc <- transcript_coverage
transcript_coverage <- round(transcript_coverage * 100 / sum(transcript_coverage), 2)
transcript_coverage_df <- data.frame(transcript_coverage)
transcript_coverage_df <- cbind(categories, transcript_coverage_df)
tc_text <- c()
for (id in 1:4) {
	tc_text[id] = paste(transcript_coverage[id], "%\n", sep="")
    tc_text[id] = paste(tc_text[id], "(", sep="")
    tc_text[id] = paste(tc_text[id], tc[id], sep="")
    tc_text[id] = paste(tc_text[id], ")", sep="")
}

# lock in factor level order
transcript_coverage_df$categories <- factor(transcript_coverage_df$categories, levels = transcript_coverage_df$categories)

t_title = paste("G: Transcript Coverage Stats (total # of transcript: ", num_of_transcripts)
t_title = paste(t_title, ")")

y_limit = 120
if (transcript_coverage[4] > 80.0) y_limit = 135
if (transcript_coverage[4] > 90.0) y_limit = 150
if (transcript_coverage[4] < 60.0) y_limit = 110

x_label <- paste("Base Coverage Percentage of Transcripts (min base coverage >= ", min_base_cutoff, sep="")
x_label <- paste(x_label, ")", sep="")

transcript_hist <- ggplot(transcript_coverage_df, aes(y=transcript_coverage, x=categories)) + 
  geom_col(color="cornsilk4", fill="cornsilk4", width=0.5) + geom_text(aes(label=tc_text), size=1.6, vjust=-0.5) +
  labs(title=t_title, x=x_label, y="Coverage Stats") +  ylim(0, y_limit) +
  theme(
    plot.title = element_text(color="blue4", size=6, face="bold.italic"),
    axis.title.x = element_text(color="black", size=5, face="bold"),
    axis.title.y = element_text(color="black", size=5, face="bold"),
    axis.text.x = element_text(size=4),
    axis.text.y = element_text(size=4)
  )

# For Gene
#
num_of_genes <- sum(gene_coverage)
gc <- gene_coverage
gene_coverage <- round(gene_coverage * 100 / sum(gene_coverage), 2)
gene_coverage_df <- data.frame(gene_coverage)
gene_coverage_df <- cbind(categories, gene_coverage_df)
gc_text <- c()
for (id in 1:4) {                                                                                          
    gc_text[id] = paste(gene_coverage[id], "%\n", sep="")
    gc_text[id] = paste(gc_text[id], "(", sep="")
    gc_text[id] = paste(gc_text[id], gc[id], sep="")
    gc_text[id] = paste(gc_text[id], ")", sep="")
}

# lock in factor level order
gene_coverage_df$categories <- factor(gene_coverage_df$categories, levels = gene_coverage_df$categories)

g_title = paste("H: Gene Coverage Stats (total # of genes: ", num_of_genes)
g_title = paste(g_title, ")")

y_limit = 120
if (gene_coverage[4] > 80.0) y_limit = 135
if (gene_coverage[4] > 90.0) y_limit = 150
if (gene_coverage[4] < 60.0) y_limit = 110

x_label <- paste("Base Coverage Percentages of Genes (min base coverage >= ", min_base_cutoff, sep="")
x_label <- paste(x_label, ")", sep="")
gene_hist <- ggplot(gene_coverage_df, aes(y=gene_coverage, x=categories)) + 
  geom_col(color="goldenrod2", fill="goldenrod2", width=0.5) + geom_text(aes(label=gc_text), size=1.6, vjust=-0.5) +
  labs(title=g_title, x=x_label, y="Coverage Stats") +  ylim(0, y_limit) +
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
grid.arrange(wgs_pct, capture_pct, wgs_hist, capture_hist, binned_graph, exon_hist, transcript_hist, gene_hist, ncol=2)

dev.off()

