# tsRTarget

# Requirements:

- Linux system, enough disk space and Ram depending on the size of RNA deep sequencing data. (Tested system: ubuntu 12.04 LTS, ubuntu 16.04 LTS)

# Installation

- Download tsRTarget pipeline package from https://github.com/zhlingl/tsRFun

  ```shell
  wget https://github.com/zhlingl/tsRFun/archive/refs/heads/main.zip
  ```

- Download necessary software, packages and reference databases as listed below:

  - The CLIP Tool Kit (CTK) (https://zhanglab.c2b2.columbia.edu/index.php/CTK_DocumentationTested version:1.1.3)
  - RNAhybrid (https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid/)
  - BLAST(https://blast.ncbi.nlm.nih.gov/Blast.cgi)( Tested version:2.10.1)
  - Reference database (See lists and download link of all pre-compiled species’ databases in Pre-compiled Databases Instruction)

- Download tsRTarget pipeline package

  - Install tsRTarget from https://github.com/zhlingl/tsRFun

    - Unpack tsRTarget package

      ```shell
      unzip main.zip
      ```

    - Attach the tsRTarget directory to your PATH:

      ```shell
      echo 'export PATH=$PATH:your_path_to_tsRFun-master/tsRTarget' >> ~/.bashrc.
      chmod 755 your_path_to_tsRTarget-master/tsRTarget/tsRTarget.sh
      ```

  - Install CTK

    - Unpack czplib-1.0.x.zip

      ```shell
      unzip czplib-1.0.x.zip
      mv czplib-1.0.x /usr/local/lib/czplib
      ```

    - Add the library path to the environment variable, so perl can find it.

      ```shell
      export PERL5LIB=your_path_to_czplib
      ```

    - Download CTK code and likewise decompress and move to whatever directory you like (as an example, we use /usr/local/)

      ```shell
      unzip ctk-1.0.x.zip
      mv ctk-1.0.x /usr/local/CTK
      ```

    - Attach the CTK directory to your PATH:

      ```shell
      echo 'export PATH=$PATH:your_path_to_CTK'
      ```

  - Install RNAhybrid

    - Unpack RNAhybrid-2.1.2.tar.gz.

      ```shell
      wget https://bibiserv.cebitec.unibielefeld.de/applications/rnahybrid/resources/downloads/RNAhybrid-2.1.2.tar.gz
      tar -xzvf RNAhybrid-2.1.2.tar.gz
      cd /path/to/ RNAhybrid-2.1.2
      ./configure
      make
      make install
      cd bin/
      ```

    - Attach the bedtools directory to your PATH:  

      ```shell
      echo 'export PATH=$PATH:your_path_to_RNAhybrid' >> ~/.bashrc
      ```

    - Start a new shell session to apply changes to environment variables:

      ```shell
      source ~/.bashrc
      ```

  - Install BLAST

    - Unpack ncbi-blast-2.10.1+-x64-linux.tar.gz from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/.

      ```shell
      tar -zxvf ncbi-blast-2.10.1+-x64-linux.tar.gz
      ```

    - Attach the samtools directory to your PATH:

      ```shell
      echo 'export PATH=$PATH:your_path_to_BLAST' >> ~/.bashrc
      ```

    - Start a new shell session to apply changes to environment variables:

      ```shell
      source ~/.bashrc
      ```

  - Download tsRFun pipeline package

    - Test if everything is installed properly:

      ```shell
      perl -v
      bash tsRTargt.sh -h
      ctk -h
      RNAhybrid -h
      blastn -h
      ```

    - If you get any error messages you should install the software or perl modules once again.

# Script description

- The input files of tsRTarget.sh are:

  ```shell
  Options:
      -i <inputfile> 		Input could be: a .fastq/.fq or .fasta/.fa file.
      -t data type clear/clash or CLIP (iclip/eclip/par-clip/hits-clip)
      -o outputdir address of annotation results
      -a index_adress
      -m mismatch default=1
      -w min tsRNA_length word_size default=14
      -n max tsRNA_length word_size default=40
      -g min target_length default=10
      -h max -g min target_length default(70 for clear clash data; 140 for clip        data)
      -s matched-length default=6
      -e gapopen default=1
      -P Collapse PCR duplicates default=T
      -S seed T/N default=T
      -R seed start default=2
      -E seed end default=8
      -M Min Free Energy default=-10
      -c match_clip the word_size of RNA-RNA target; only used for CLIP seq 		   data,
  Others:
  	-h print this usage message
  ```

- Example of use:

  ```shell
  bash tsRTarget.sh -i example.fa -t clear -o example_result -a hg38_index
  ```

  