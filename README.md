<img src="BCM-HGSC-logo.png" width=250>

---

### Table of Contents

- [Description](#description)
- [Build Scandium](#build-scandium)
- [How To Use](#how-to-use)
- [Manuscript Abstract](#Manuscript-Abstract)
- [Contact Info](#Contact-Info)

---

## Description

Scandium is a software program that was developed at Baylor College of Medicine Human Genome Sequencing Center. It produces a set of quantitative and quality metrics, especially the uniformity metric to evaluate Next-Generation sequencing runs for quality control purposes.

[Back To The Top](#Table-of-Contents)

## Build Scandium

See [INSTALL](INSTALL) for complete details. Please download the [release tarballs](https://github.com/phuang7008/kaggle/releases) of your choice and build Scandium using the following steps:

    wget https://github.com/phuang7008/Scandium/releases/download/Scandium_vXXX/Scandium-vXXX.tar.gz
    tar zxvf Scandium-vXXX.tar.gz 

    build htslib (see htslib install instruction at the htslib website)
    
    autoreconf -i
    ./configure --includedir="path-to-mysql-header" LDFLAGS=$LDFLAGS -L"path-to-libmysqlclient"
    make
    make install or make prefix=<your dir choice> install'

[Back To The Top](#Table-of-Contents)

## How To Use

For Scandium v1 series, here is an example of the run command: 

    scandium-v1+ -i input_bam -o output_dir -R reference -t target_bedfile -f annotation_file -r list_of_chromosomes_to_be_processed -T num_of_threads -b minimal_base_qual -m minimal_mapping_qual -D reference_version -w WGS_coverage_summary -N Ns_regions_file

For Scandium v3 series, here is an example of the run command: 

    scandium-v3+ -i input_bam -o output_dir -R reference -r list_of_chromosomes_to_be_processed -T num_of_threads -b minimal_base_qual -m minimal_mapping_qual -D reference_version -w WGS_coverage_summary -t pair_of_target_bedfile_vs_annotation_file_separated_by_tab_one_per_line

[Back To The Top](#Table-of-Contents)

## Manuscript Abstract

<img src="Scandium_Abstract_new.png" width=1250>

[Back To The Top](#Table-of-Contents)

## Contact Info

If users encounter any issues, please contact pemhuang@gmail.com

[Back To The Top](#Table-of-Contents)
