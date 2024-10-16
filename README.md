<img src="images/BCM-HGSC-Logo.png" width=250>

### Table of Contents

- [Description](#description)
- [Installation](#installation)
- [How To Use](#how-to-use)
- [Manuscript](#manuscript)
- [Contact Info](#Contact-Info)

---

## Description

Scandium is a software program that was developed at Baylor College of Medicine Human Genome Sequencing Center. It produces a set of quantitative and quality metrics to evaluate Next-Generation sequencing runs for quality control purposes.

[Back To The Top](README.md)

---

## Installation

#### Requirements

Building Scandium requires a few programs and libraries to be present.

Check to ensure the system provides the followings

    GNU make
    Both C and C++ compiler (e.g. gcc/g++ or related modules)

The build requires autotools scripts like,

    autoheader
    autoconf
    autoreconf

The following external libraries are needed for the Scandium installation

    libbz2
    libcrypto
    libcurl
    libdl
    liblzma
    libm
    libmysqlclient
    libpthread
    librt
    libssl
    libz

#### Building Configure

Before the installation of Scandium, users need to first install 
    [htslib](https://github.com/samtools/htslib)
package

Scandium has been tested against htslib from version 1.10 to 1.19.1, but Do Not use v1.12 as there is a bug in v1.12.

Important Note:
htslib must be installed inside scandium directory. Please create a 'htslib' 
sub-directory like scandium/htslib if the 'htslib' sub-directory doesn't exist.
The 'htslib' sub-directory name should not include any version numbers.

Please follow the htslib 'INSTALL' instructions all the way to the final 
step 'make install' such as

    make prefix=<path to scandium/htslib directory> install

Once htslib is successfully installed, users can proceeds to build Scandium.
Please download any release package you prefer and decompress the package

    tar zxvf Scandium-release-version.tar.gz
    cd Scandium-release-version

Run the following to create a brand new configure for your system

    autoreconf -i

Once users have generated the configure file, then do the following

    ./configure --includedir="path-to-mysql-header" LDFLAGS=$LDFLAGS -L"path-to-libmysqlclient"
    make
    make install

The './configure' generates Makefiles for the compilation.

The 'make' compiles source code in both src/ and uniformity_graph/
and generates executables.

The 'make install' command installs the compiled executables into the /usr/local
directory. Users can change the installation location by adding --prefix=DIR
option to the ./configure run or via 'make prefix=DIR install'.

#### NOTES

In order to generate uniformity graph for each chromosome of your choice, 
please ensure your system has 'gnuplot' installed. 

#### PERL

To generate gene annotation from the scripts provided, users need to have
MariaDB (mysql) database and DBD/mysql.pm module installed.

#### PYTHON

To generate gene annotation from the scripts provided, users need to have
python's mysql and mysql-connector installed.

#### BEDOPS

The scripts used to generate gene annotation also need to use bedops. Please have
it installed.

[Back To The Top](README.md)

---

## How To Use

For Scandium v1 series, here is an example of the run command:
    scandium-v1+ -i input_bam -o output_dir -R reference -t target_bedfile -f annotation_file -r list_of_chromosomes_to_be_processed -T num_of_threads -b minimal_base_qual -m minimal_mapping_qual -D reference_version -w WGS_coverage_summary -N Ns_regions_file

For Scandium v3 series, here is an example of the run command:
    scandium-v3+ -i input_bam -o output_dir -R reference -r list_of_chromosomes_to_be_processed -T num_of_threads -b minimal_base_qual -m minimal_mapping_qual -D reference_version -w WGS_coverage_summary -t pair_of_target_bedfile_vs_annotation_file_separated_by_tab_one_per_line

[Back To The Top](README.md)

---

## 2019 ASHG Poster

![Scandium 2019 Poster](images/2019_ASHG_poster_PeterHuang_v12_PC.png)

[Back To The Top](README.md)

---

## Contact Info

If users encounter any issues, please contact pemhuang@gmail.com

[Back To The Top](README.md)
