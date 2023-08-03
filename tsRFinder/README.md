# tsRFinder

# Requirements:

- Linux system, enough disk space and Ram depending on the size of RNA deep sequencing data. (Tested system: ubuntu 12.04 LTS, ubuntu 16.04 LTS)

# Installation  

- Download tsRFinder pipeline package from https://github.com/zhlingl/tsRFun  

  ```shell
  wget https://github.com/zhlingl/tsRFun/archive/refs/heads/main.zip
  ```

- Download necessary software, packages and reference databases as listed below:

  - Perl 5 (https://www.perl.org) (Tested version: v5.14.2, v5.22.1); Perl 5 might be already installed in the linux system.
  - Bowtie [1] (http://bowtie-bio.sourceforge.net/index.shtml) (Tested version: 1.1.1, 1.2.1.1)
  - Bedtools(https://bedtools.readthedocs.io/en/latest/) (Tested version: 2.25.0, 2.29.2)
  - Samtools(http://www.htslib.org/)( Tested version:1.9,1.13)
  - Reference database (See lists and download link of all pre-compiled species’ databases in Pre-compiled Databases Instruction)

- Download tsRFinder pipeline package 

  - Install tsRFinder from https://github.com/zhlingl/tsRFun

    - Unpack tsRFinder package  

      ```shell
      unzip tsRFun-master.zip
      ```

    - Attach the tsRFinder directory to your PATH:  

      ```shell
      echo 'export PATH=$PATH:your_path_to_tsRFun-master/tsRFinder' >> ~/.bashrc
      chmod 755 your_path_to_tsRFun-master/tsRFinder/tsRFinder.pl
      ```

  - Install Bowtie  

    - Unpack bowtie-1.x.x-linux-x86_64.zip.  

      ```shell
      unzip bowtie-1.x.x-linux-x86_64.zip
      ```

    - Attach the bowtie directory to your PATH:  

      ```shell
      echo 'export PATH=$PATH:your_path_to_bowtie' >> ~/.bashrc 
      ```

    - Start a new shell session to apply changes to environment variables:  

      ```shell
      source ~/.bashrc
      ```

  - Install bedtools  

    - Unpack bedtools-2.25.0.tar.gz.  

      ```shell
      wget https://github.com/arq5x/bedtools2/archive/v2.25.0.tar.gz
      tar xzvf v2.25.0
      cd bedtools2-2.25.0/
      make
      cd bin/
      ```

    - Attach the bedtools directory to your PATH:  

      ```shell
      echo 'export PATH=$PATH:your_path_to_bedtools' >> ~/.bashrc
      ```

    - Start a new shell session to apply changes to environment variables:  

      ```shell
      source ~/.bashrc
      ```

  - Install samtools  

    - Unpack bedtools-2.25.0.tar.gz.  

      ```shell
      wget -c https://github.com/samtools/samtools/releases/download/1.9/samtools-
      1.9.tar.bz2
      tar jxvf samtools-1.9.tar.bz2
      cd samtools-1.9/
      ./configure --prefix= your_path_of_samtools(example:
      /home/vip47/biosoft/samtools-1.9)
      make
      make install
      ```

    - Attach the samtools directory to your PATH:  

      ```shell
      echo 'export PATH=$PATH:your_path_to_ samtools' >> ~/.bashrc
      ```

    - Start a new shell session to apply changes to environment variables:  

      ```shell
      source ~/.bashrc
      ```

  - Install perl modules  

    ```
    cpan -i Net::Server;
    cpan -i Getopt::Std;
    cpan -i File::Find;
    cpan -i File::Basename;
    cpan -i Cwd;
    cpan -i Math::CDF;
    cpan -i FileHandle;
    ```

- Download tsRFinder pipeline package  

  - Test if everything is installed properly:  

    ```shell
    perl -v
    tsRFinder.pl -h
    bowtie
    samtools
    bedtools
    ```

- If you get any error messages you should install the software or perl modules once again.  

# Script description

- The input files of tsRFinder.pl are:  

  ```shell
  Options:
      -i <file> Input could be: .fastq/.fq or .fasta/.fa file.
      -o<file> Output address of annotation results
      -t <int> Number of threads to launch (default = 4)
      -x <str> Address of bowtie index tRNA information
  Alignment:
      -l <int> The minimal length of the output sequences (default = 15)
      -L <int> The maximal length of the output sequences (default = 45)
      -M <int> The total number of mismatches in the entire alignment (default = 		  0)
      -p <float> The p-value threshold to determine whether the fragment is 		   tsRNA.
  Others:
      -v Print version information
      -h Print this usage message
  ```

  

# Example of use:  

```
tsRFinder.pl -i PATH_of_example/fasta -o PATH_of_example/ -x hg38_index/
```

# Output file:  

The information columns of the output file are as follows:

1. Type of tsRNA, such as tRF-3, tRF-5, and so on.
2. RPM (Reads Per Million) value of the tsRNA.
3. AlignScore value.
4. Length of the tsRNA.
5. tRNA from which the tsRNA originates.
6. Position of the tsRNA on its source tRNA.
7. Sequence of the tsRNA.
