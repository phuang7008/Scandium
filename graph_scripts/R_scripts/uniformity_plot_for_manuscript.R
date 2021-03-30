# this script is ussed to plot the distribution (barplot) for genome coverage data
# this is a simplied graph as it will draw a single point for each block no matter how big the block would be
#
require(ggplot2)

# creating two coverage vectors
df1 <- c()
df2 <- c()

coverage_plot <- function(fName, chrom_id, chrom_len, type) {
  Coverage <- c(integer(chrom_len))

  # read from file and process it accordingly to retrieve the coverage information
  con <- file(fName, "r")
  items = c()
  start = 0
  end = 0

  while (TRUE) {
    line = readLines(con, n=1)
    if (length(line) == 0)
      break
    
    items <- unlist(strsplit(line, '\t'))
	items[1] = gsub("^chr","", items[1])
	items[1] = gsub("^Chr","", items[1])
	items[1] = gsub("^CHR","", items[1])

    if (items[1] == chrom_id) {
      #print(items[1])
      start = as.numeric(items[2])
      len   = as.numeric(items[3]) - start + 1

      if (length(items) > 4) {
        cov = as.numeric(items[5])
      } else {
        cov = as.numeric(items[4])
      }

      Coverage[start] <- as.integer(cov + 0.5)
    } 
  }

  Position <- seq(1, chrom_len, by=1)
  if (type == 1) df1 <<- data.frame(Position, Coverage)
  if (type == 2) df2 <<- data.frame(Position, Coverage)

  close(con)
  
}

# get the command line arguments
options(echo=TRUE)	# if you want to see commands in output file
#options()
args <- commandArgs(trailingOnly=TRUE)
print(args)
file_in1 <- args[1]
file_in2 <- args[2]
chrom_id_in <- args[3]
version <- args[4]

# define chromosome id name vector
id_names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", 
              "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT")

# define lengths of each chromosome for version HG19 (v37) and v38
v <- c()

if (version == "hg38") {
    v <- c( "1"=248956422, "2"=242193529, "3"=198295559, "4"=190214555, "5"=181538259, "6"=170805979, "7"=159345973,
            "8"=145138636, "9"=138394717, "10"=133797422, "11"=135086622, "12"=133275309, "13"=114364328, "14"=107043718,
            "15"=101991189, "16"=90338345, "17"=83257441, "18"=80373285, "19"=58617616, "20"=64444167, "21"=46709983,
            "22"=50818468, "X"=156040895, "Y"=57227415, "MT"=16569)
} else {
    v <- c( "1"=249250621, "2"=243199373, "3"=198022430, "4"=191154276, "5"=180915260, "6"=171115067, "7"=159138663,
            "8"=146364022, "9"=141213431, "10"=135534747, "11"=135006516, "12"=133851895, "13"=115169878, "14"=107349540,
            "15"=102531392, "16"=90354753, "17"=81195210, "18"=78077248, "19"=59128983, "20"=63025520, "21"=48129895,
            "22"=51304566, "X"=155270560, "Y"=59373566, "MT"=16569)
}

for (id in id_names) {
  if (chrom_id_in == id) {
    #coverage_plot(file_in1, 7, 159138663, 1)
    coverage_plot(file_in1, id, v[id], 1)
    #coverage_plot(file_in2, 7, 159138663, 2)
    coverage_plot(file_in2, id, v[id], 2)
    gc()
    break;
  }
}

print("Finishing convert to data frame")
  
cov_plot <- ggplot(NULL, aes(x=Position, y=Coverage)) + 
		geom_point(data=df1, aes(color="orangered3"), size=0.05) +
		geom_point(data=df2, aes(color="green4"),  size=0.05) +
		scale_color_manual(labels = c("Sample A Uniformity Metric: 0.958", "Sample B Uniformity Metric: 0.269"), values=c("green4","orangered3"), name="Legend") +
		guides(colour = guide_legend(title.hjust = 0, keywidth=0.1, keyheight=1.2, override.aes = list(size=2.0))) +
		labs(title="Coverage Uniformity Comparison (Human Chromosome 7)") + ylim(1, 120) + theme_bw() +
		theme(
			plot.title = element_text(color="black", size=17, face="bold.italic"),
			axis.title.x = element_text(color="black", size=15, face="bold"),
			axis.title.y = element_text(color="black", size=15, face="bold"),
			axis.text.x = element_text(size=14,face="bold"),
			axis.text.y = element_text(size=14,face="bold"),
			legend.title=element_text(size=14,face="bold"),
			legend.text=element_text(size=14,face="bold"),
			legend.justification=c(0.9,0.9),
			legend.position=c(0.9, 0.9)
	)
  
png(filename = "comparison7.png", width=19, height=3.5, units="in", res=200)
plot(cov_plot)
dev.off()

rm(Coverage1)
rm(df1)
