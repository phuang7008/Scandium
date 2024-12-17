<img src="images/BCM-HGSC-Logo.png" width=250>

---

### Table of Contents

- [Description](#description)
- [Build Scandium](#build-scandium)
- [How To Use](#how-to-use)
- [Contact Info](#Contact-Info)

---

## Description

Scandium is an open-source software providing comprehensive sequence quality control and coverage metrics for short-read NGS data from BAM/CRAM files. These include a novel whole-genome coverage uniformity metric with user-friendly coverage summary and uniformity plots for easy visualization. Scandium provides annotations and calculates coverage percentages of user-defined inputs at various granularities from SNPs to exon to gene based on a user-defined minimal base coverage threshold to aid clinical studies. Because Scandium requires minimal computational resources, it is suitable for large sequencing centers or core laboratories. It can be run on desktop terminals, local in-house clusters or cloud-computing platforms such as AWS.

Note: Scandium has only been tested on Linux-like systems. Thanks!

[Back To The Top](#Table-of-Contents)

## Build Scandium

See [INSTALL](INSTALL) for complete details. Please download the [release tarballs](https://github.com/phuang7008/Scandium/releases) of your choice and build Scandium using the following steps:

    wget https://github.com/phuang7008/Scandium/releases/download/Scandium_vXXX/Scandium-vXXX.tar.gz
    tar zxvf Scandium-vXXX.tar.gz 

build htslib (see htslib install instruction at the htslib website)
    
    autoreconf -i
    ./configure --includedir="path-to-mysql-header" LDFLAGS=$LDFLAGS -L"path-to-libmysqlclient"
    make
    make install or make prefix=<your dir choice> install'

[Back To The Top](#Table-of-Contents)

## How To Use

To run Scandium, here is an example of the run command: 

    scandium -i input_bam -o output_dir -R reference -r list_of_chromosomes_to_be_processed -T num_of_threads -b minimal_base_qual -m minimal_mapping_qual -D reference_version -w WGS_coverage_summary -t target_bedfile -f annotation_bedfile

Note: Starting from version 3+, the -f option is no longer available. Instead, use the -t option with a file that contains both the target BED file and annotation BED file pairs, separated by a tab, one pair per line.

Note: The resources/ folder contains pre-built bedfiles as examples to run Scandium on human genomes. Files are available for both hg37 and hg38 genome builds. Be sure to select the correct version to match your dataset.

- The human genome references can be downloaded from the following links:

        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa .
        wget https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/references/GRCh37/hs37d5.fa.gz .

- For the annotations, you can choose the Matched Annotation from NCBI and EMBL-EBI (MANE) transcripts or any other user-defined annotations.

[Back To The Top](#Table-of-Contents)

## Contact Info

If users encounter any issues, please contact pemhuang@gmail.com

[Back To The Top](#Table-of-Contents)
