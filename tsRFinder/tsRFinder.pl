#!/usr/bin/perl


##################################################
# tsRFinder.pl find tsRNAs from deep-sequencing data
# Lingling Zheng(zhengll33@mail.sysu.edu.cn)
# Data: 2021-9-2
##################################################

use warnings;
use strict;
use Getopt::Std;
use Cwd 'abs_path';
use File::Find;
use File::Basename;
use Cwd;
use Math::CDF;
use FileHandle;

print  "

###################################
#                                 #
# tsRNA_exp_in_small-RNA data     #
#                                 #
#    last change:2/9/2021         #
#                                 #
################################### 

";

use vars qw($opt_i $opt_l $opt_L $opt_M $opt_t $opt_x $opt_p $opt_o $opt_v $opt_h);
getopts('i:l:L:M:t:x:p:o:v:h');
my $inputfile		= $opt_i;
my $MintRF		= $opt_l ? $opt_l : 15;
my $tRFLarge		= $opt_L ? $opt_L : 45;
my $mismatchNumber		= $opt_M ? $opt_M : 0;
my $thread		= $opt_t ? $opt_t : 4;
my $index_address	= $opt_x;
my $pMin = $opt_p ? $opt_p : 0.05;
my $outputdir	= $opt_o;
my $version		= $opt_v ? 1 : 0;
my $help		= $opt_h ? 1 : 0;

my $usage="
Description:	Perl script used to find tsRNA from small RNA sequencing data.
		Query input files support sequence formats: .fastq/.fq, .fasta/.fa.

Usage Examples: tsRFinder.pl -i fasta/fastq -o outputfile

The input files are:

Options:
  -i <inputfile>	Input could be: 
		  a .fastq/.fq or .fasta/.fa file. 
  -o output address of annotation results
  -t <int>	number of threads to launch (default = 4)
  -x <str>	address of bowtie index tRNA information 

Alignment:
  -l <int>	the minimal length of the output sequences (default = 15)
  -L <int>	the maximal length of the output sequences (default = 45)
  -M <int>	the total number of mismatches in the entire alignment (default = 0)
  -p <float> the p-value threshold to determine whether the fragment is tsRNA.
Others:
  -v		print version information
  -h		print this usage message

Example of use: 
tsRFinder.pl -i PATH_of_example/fastq -o PATH_of_example/ -x PATH_of_example/index/

";

my $version_info = "1.1.1";

if ($version) {
	print "\ntsRFun version: $version_info\n";
	exit;
}elsif($help) {
	print $usage;
	exit;
}

##determine if input file and genome file are defined

unless (defined $inputfile && 	defined $index_address){
	print "\nInput file and index address should be specified!\n\n";
	print $usage;
	exit;
}

unless (-e $inputfile){
	print "\nInput file is not exist!\n\n";
	print $usage;
	exit;
}

my $input_address = abs_path($inputfile);

if (-f $input_address){
	my @input = split(/\//, $input_address);
	pop (@input);
	$input_address = join('/', @input) . '/';
}

unless (defined $outputdir){
	$outputdir = $input_address;
}

print "\ntsRFun version: $version_info\n\n";
print "Please cite: .\n\n";
print "Input file address: $input_address\n";
print "output file address: $outputdir\n";
print "Reference address: $index_address\n";

my @array = split(/\//,$inputfile); 
my @array1 = split(/\./,pop(@array)); 
my $barcode = $array1[0];

my $tRNAMatureindex = $index_address."bowtie/mature";
my $tRNAMatureFasta = $index_address."hg38-tRNA.fa";
my $tRNAPRIindex = $index_address."bowtie/pri";
my $tRNAPRIFasta = $index_address."hg38-pri50tRNA.fa";
my $tRNAMatureChromeSize = $index_address."hg38-tRNA.sizes";
my $tRNAPRIChromeSize = $index_address."hg38-pri50tRNA.sizes";
my $tRNAposition = $index_address."tRNA_position.txt";
my $tsRNA_TCGA = $index_address."tsRNA_TCGA.txt";

my $dir = $outputdir."Result";

unless (-e $dir) {
    system ("mkdir $dir");
}

###############################################################################
##
## Don't change anything below this line
##
###############################################################################

## measuring times
my ($sTime,$eTime,$stime,$etime,$sTimeG, $eTimeG, $stimeG, $etimeG);
my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings);

($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
$second = "0$second" if($second =~ /^\d$/);
$sTimeG = "$hour:$minute:$second";
$stimeG = time;

print  "handle_raw_reads started at $sTimeG\n\n\n";

#############################################################Get Options##############################

my $unknownReads=$dir."/unknown_reads.list";
my $tRNAMaturebedGraph=$dir."/tRNAMature.bg";
my $tRNAPRIbedGraph=$dir."/tRNAPri.bg";
my $tRNAMatureBam=$dir."/".$barcode."_tRNAMature_bam";
my $tRNAPRIBam=$dir."/".$barcode."_tRNAPri_bam";
my $tRNAMaturebed=$dir."/tRNAMature.bed";
my $tRNAPRIbed=$dir."/tRNAPri.bed";

my $unmap2MatureFile   = $dir."/unmap2Mature.fa"; 
my $unmap2PriFile      = $dir."/unmap2Pri.fa";
my $reads2tRNAMature   =$dir."/reads2tRNAMature.sam";
my $reads2tRNAPri      =$dir."/reads2tRNAPri.sam";
my $tsRFfinderOutFile   =$dir."/".$barcode."_tsRfinderResult_pMin_".$pMin.".txt";
my $logFile            = $dir . "/tsRfinder.log";
my $logFileHandle      = FileHandle->new(">$logFile");


my $ErrorInfo;

my %hash=();
my %hashMod=();
my $InputFile2 = $dir."/input2.fa";


####################################### 1.change fastq/fastq file to collapsed fasta format (awk)  ###########################################

open (IN,$inputfile) || die "Can not open intput file $inputfile";
  my $line=<IN>;
    chomp $line;

  if($line =~ /^@/){
    my $cmd = " awk \'BEGIN{flag=0;}{if(\$0 ~ /\^@/) {flag=1;k=1;} else {if(flag) {hash[\$0]++; flag=0;}} } END{for (i in hash)  {print \">seq_\"k\"_x\"hash[i]\"\\n\"i; k++;}}\' ".$inputfile." > ". $InputFile2;
    print " 3.change fastq to collaspsed fasta\n",$cmd,"\n";
    system $cmd;
  }
  elsif($line =~ /^>seq_/){
    my $cmd = "cp $inputfile  $InputFile2";
    system $cmd;
  }
  elsif($line =~ /^>/){
    my $cmd = " awk \'BEGIN{flag=0;}{if(\$0 ~ /\^>/) {flag=1;k=1;} else {if(flag) {hash[\$0]++; flag=0;}} } END{for (i in hash)  {print \">seq_\"k\"_x\"hash[i]\"\\n\"i; k++;}}\' ".$inputfile." > ". $InputFile2;
    print " 3.change fasta to collasped fasta\n",$cmd,"\n";
    system $cmd;

  }
  else{
    print " please print fastq or fasta file\n";
  }
close IN;

$ErrorInfo = CheckInputFile($InputFile2);

my $InputFileLength = $InputFile2."_length";
LengthLimit($InputFile2,$InputFileLength,$MintRF,$tRFLarge);

my %tRNAHash=();
my %tRNA;
my %tRNAMatureLen;
my %tRNAPriLen;
my %tRNAmatureSeq;
my %tRNAPriSeq;
my %tRNAposition;

if ($ErrorInfo){ 
  #If there is error information, the specific content will be output according to the error information
  &HTMLPrint($ErrorInfo);
}
else{
  main();
}

sub main{
  ####################################### 2.mapping reads to tRNA reference  ###########################################
  GettRNANames($tRNAMatureFasta);
  gettRNALen($tRNAMatureFasta,$tRNAPRIFasta);
  gettRNAposition($tRNAposition);

  my ($MatureMapped,$PriMapped)=0;
  $MatureMapped=ExeBowtie($InputFile2,$unmap2MatureFile,$tRNAMatureindex,$reads2tRNAMature,$mismatchNumber,$logFile,"normal");
  $PriMapped=ExeBowtie($unmap2MatureFile,$unmap2PriFile,$tRNAPRIindex,$reads2tRNAPri,$mismatchNumber,$logFile,"normal");

  chomp($MatureMapped);
  chomp ($PriMapped);
  print $logFileHandle "There are $MatureMapped reads mapped to mature tRNA\n";
  print $logFileHandle "There are $PriMapped reads mapped to pri tRNA\n";

  SamtoBed($reads2tRNAMature,$tRNAMaturebed,$logFileHandle);
  BedtoBedGraph($tRNAMaturebed,$tRNAMatureChromeSize,$tRNAMaturebedGraph,$logFileHandle);
  SamtoBam($reads2tRNAMature,$tRNAMatureBam,$logFileHandle);

  SamtoBed($reads2tRNAPri,$tRNAPRIbed,$logFileHandle);
  BedtoBedGraph($tRNAPRIbed,$tRNAPRIChromeSize,$tRNAPRIbedGraph,$logFileHandle);
  SamtoBam($reads2tRNAPri,$tRNAPRIBam,$logFileHandle);

  ####################################### 3.classify small RNAs fragment  ###########################################
  ####################################### 4.quantify tsRNAs  ###########################################
  Quantify_tRNAExp($tRNAMaturebedGraph,$tRNAPRIbedGraph,$pMin,$MatureMapped,$PriMapped,$reads2tRNAMature,$reads2tRNAPri);

  ####################################### 5.write all the information to file  ###########################################
  printHTML($tsRFfinderOutFile);
  my $cmd = "cp $tsRFfinderOutFile $outputdir";
  system ($cmd);
  system ("rm -rf $dir");

}


sub CheckInputFile{
  my ($file) = @_;
  my $result=0;

  if (-z $InputFile2) { 
    #If the file is empty, the upload fails and error 1 is returned
    $result =1;
    return $result;
  }

  open (IN,$file) || die "Can not open intput file $file";
  my $line=<IN>;
  chomp $line;

  if($line !~ /^>\S+/){
    $result = 3;
    return $result;
  }

  if($line =~ /\s/){
    $result = 4;
    return $result;

  }

  if(<IN> !~ /^[ACGTUNacgtun]*$/){
    $result = 5;
    return $result;
  }

  if($line !~ /^>\S+\_\S+\_x\d+/){
    $result = 6;
    return $result;
  }
}

sub ExeBowtie{
  my ($fastaFile_collapse,$unmapFile,$index,$outFile,$mismatchNumber,$logFileHandle,$type)=@_;
  my $cmd;
  if ($type eq "large") {
    
    $cmd = "/public/home/wangzh/data/bowtie-1.2.3-linux-x86_64/bowtie --large-index -p ".$thread." ". $index." -f ".$fastaFile_collapse." -v ".$mismatchNumber." -a --best --strata -S --norc  > ".$outFile  ." --un ". $unmapFile;
  }
  elsif ($type eq "normal"){
     $cmd = "/public/home/wangzh/data/bowtie-1.2.3-linux-x86_64/bowtie -p ".$thread." ". $index." -f ".$fastaFile_collapse." -v ".$mismatchNumber." -a --best --strata -S --norc  > ".$outFile  ." --un ". $unmapFile;
  }
  else {die "do not know the type $type, this program can only take normal and large index type\n";}
 
  system $cmd;
  return StatMappedNumber($outFile,$logFileHandle);
}

sub StatMappedNumber{
  my ($file) = @_;
  my $cmd = "samtools view -@ 96 -bS ". $file ."| samtools view -@ 96 -F 4 - | awk \'BEGIN{FS=\"x\"} { if(\$0 ~ /^\\\[/) {next;} else { total +\=\$2;}} END{print total;}'";

  my $number = `$cmd`;
  return $number;
}

sub SamtoBed{
  my ($input,$output,$logFileHandle)=@_;
  my $command = 'samtools view -@ 96 -Sb '.$input.' | bamToBed -i stdin > '.$output;
    print $logFileHandle $command,"\n";
    system($command);

}

sub GettRNANames{
  my ($file)=@_;
  my @tRNA=("tRF-5","tRF-3","tRF-1","tRF-i");
  open (IN,$file) || die "can not open input file $file\n";
  while (<IN>) {
    chomp $_;
    if ($_ =~ /^>(\S+)$/) {
      my $name=$1;
      foreach my $type (@tRNA) {
        $tRNAHash{$name."_".$type}{0}=$barcode;
        $tRNAHash{$name."_".$type}{1}=0; #Store tsRNA expression
        $tRNAHash{$name."_".$type}{2}="-";  #Store tsRNA position
      }
    }
    
  }
  close IN;
}

sub gettRNALen{
  my ($tRNAMatureFasta,$tRNAPRIFasta) = @_;
  open (IN,$tRNAMatureFasta) || die "can not open $tRNAMatureFasta";
  my $name;
  while (<IN>) {
    chomp $_;
    if ($_ =~ /^>(\S+)$/) {
      $name=$1;
    }
    else {
      $tRNAMatureLen{$name}=length($_);
      $tRNAmatureSeq{$name}=$_;
    }
  }
  close IN;

  open (IN,$tRNAPRIFasta) || die "can not open $tRNAPRIFasta";
  
  while (<IN>) {
    chomp $_;
    if ($_ =~ /^>(\S+)$/) {
      $name=$1;
    }
    else {
      $tRNAPriLen{$name}=length($_);
      $tRNAPriSeq{$name}=$_;
    }
  }
  close IN;
}

sub gettRNAposition{
  my ($tRNAposition) = @_;

  open (IN,$tRNAposition) || die "can not open $tRNAposition";

  while (<IN>) {
    chomp $_;
    my ($k,$v)=split (/\t/,$_);
    #Use "=" to split each row. Perl is stored in the default after each row is in read $_ 
    
    $tRNAposition{$k}=$v;
    #Save to hash
  }
  close IN;

}

sub LengthLimit{
  my ($InputFile,$InputFileLength,$MintRF,$tRFLarge) = @_;
  open(IN, $InputFile) || die "can not open $InputFile";
  open(OUT, ">$InputFileLength") || die "can not open $InputFileLength"; 
  my ($name,$seq);
  while (<IN>) {
    chomp $_;
    if ($_=~ /^>/) {
      $name=$_;
    }
    else {$seq=$_;
      if ((length($seq) < ($tRFLarge +1)) && (length($seq) > ($MintRF -1) )) {
        print OUT $name,"\n",$seq,"\n";
      }
    }
  }
  close IN; 
  close OUT;
}


sub SamtoBam{
  my($inputFile,$outputFile,$logFileHandle) =@_;
  my $OriginBam = $inputFile."_bam";
  my $command = " samtools view -@ 96 -bS ".$inputFile." > ".$OriginBam;
  print $logFileHandle $command,"\n";
  system($command);

  $command = " samtools sort ".$OriginBam." -o ".$outputFile;
  print $logFileHandle $command,"\n";
  system($command);

  $command ="samtools index ".$outputFile; 
  #Samtools sort will automatically add. BAM to the end of the file

  print $logFileHandle $command,"\n";
   system($command);
}

sub Quantify_tRNAExp{
  my ($tRNAMaturebedGraph,$tRNAPRIbedGraph,$pMin,$MatureMapped,$PriMapped,$reads2tRNAMature,$reads2tRNAPri)=@_;

  handleMature($pMin,$MatureMapped);
  handlePri($pMin,$PriMapped);
}

sub handleMature{
  my ($pMin,$MatureMapped)=@_;

  ##Read the position of tRNA anti codon ring and judge the type of subsequent tiRNA
  my %tRNAposition;

  open (IN,$tRNAposition) || die "can not open $tRNAposition";

  while (<IN>) {
    chomp $_;
    my ($k,$v)=split (/\t/,$_);
    $tRNAposition{$k}=$v;#save to hash
  }
  close IN;

  open (IN,$tRNAMaturebedGraph);
  while (<IN>) {
    chomp $_;
    my @cols = split(/\s+/,$_);
    my $Tname=$cols[0];
    my $Tlength=$tRNAMatureLen{$Tname};

    my $start = $cols[1];
    my $end = $cols[2];
    my $Snum=$cols[3];

    $hash{$Tname}{"Maturetotal"}+=$Snum;    #The total expression of reads on tRNA mate was recorded
    $hash{$Tname}{"MatureLength"}=$Tlength;

    for(my $i=$start; $i<$end; $i++){
      $hash{$Tname}{"Mature"}{$i}+=$Snum;
    }

  }
  close IN;

  my %cluster=();
  my $tsRNA_sum;
  # hash  GettRNANames; TRNA precursor and mature information read in the module, key is tRNA
  foreach my $key (keys %hash){
    my $Tlength = $hash{$key}{"MatureLength"};
    if($Tlength==0){
      next;
    }

    my $Total = $hash{$key}{"Maturetotal"};
    if($Total==0){
      next;
    }
    %cluster=();
    my $largestNum=0; 
    ## Determine whether the tsrna fragment conforms to the binomial distribution and obtain the start and end position of tsRNA
    for (my $i=0; $i< $Tlength; $i++){
      if (defined $hash{$key}{"Mature"}{$i}) {
        my $Snum=$hash{$key}{"Mature"}{$i};
        my $pvalue=Math::CDF::pbinom($Total-$Snum,$Total,1- 1/$Tlength);

        if ($pvalue <$pMin) {
          
          unless (defined $cluster{"start"}) {
            $cluster{"start"}=$i;
            $cluster{"end"}=$i;
            next;
          }

          if ($i-$cluster{"end"} == "1") {
            $cluster{"end"}=$i;
          }
        }
        else {
          unless (defined $cluster{"start"}) {
            %cluster = ();
            next;
          }

            my $region = $key.":".$cluster{"start"}."-".$cluster{"end"};
            my $regionNumber = 0;

            for (my $i = $cluster{"start"}; $i <= $cluster{"end"}; $i++) {
              $regionNumber+=$hash{$key}{"Mature"}{$i};
            }

            if ($regionNumber > $largestNum) {
              $largestNum = $regionNumber;
            }
            my $len=$hash{$key}{"MatureLength"};
            my $type="";


            my $SigCov =getCoverage($tRNAMatureBam,$key,$region,$regionNumber); 

            my ($tRNA_temp,$position_temp)=split (/:/,$SigCov);
            my ($position_start,$position_end)=split (/-/,$position_temp);
            my ($steamstart,$steamend)=split (/-/,$tRNAposition{$key});
            my $regionScore = 0;

            if (($position_start<2) || ($len-$position_end<2)) {

              if((($steamend - $position_end +3>0) && ($steamstart - $position_end <3))&& ($SigCov)){
                $type="tiRNA-5";
                $region=$SigCov;
                $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$MatureMapped;
                $hash{$key}{$type}[1]=$region;
                $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
              }
              elsif((($steamend- $position_start +3>0) && ($steamstart - $position_start <3))&& ($SigCov)){
                $type="tiRNA-3";
                $region=$SigCov;
                $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$MatureMapped;
                $hash{$key}{$type}[1]=$region;
                $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
              }
              elsif (($position_start<2) && $SigCov){
                $type="tRF-5";
                $region=$SigCov;
                $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$MatureMapped;
                $hash{$key}{$type}[1]=$region;
                $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
              }
              elsif( ($len-$position_end<2) && $SigCov){
                $type="tRF-3";
                $region=$SigCov;
                $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$MatureMapped;
                $hash{$key}{$type}[1]=$region;
                $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
              } 
            }
            #If it is neither tRF-5 nor tRF-3, but the expression is the highest of all tRFs in the tRNA, it is named tRF-i
            else {
              $type ="tRF-i";
              $region=$SigCov;
              $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$MatureMapped;
              $hash{$key}{$type}[1]=$region;
              $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
            }
            %cluster=();
          }
      }
    }

    if (defined $cluster{"start"}) {

      my $region = $key.":".$cluster{"start"}."-".$cluster{"end"};
      $tsRNA_sum=$hash{$key}{"Mature"}{$cluster{"end"}};

      my $regionNumber = 0;
      my $regionScore = 0;
      for (my $i = $cluster{"start"}; $i <= $cluster{"end"}; $i++) {
        $regionNumber+=$hash{$key}{"Mature"}{$i};
      }


      if ($regionNumber > $largestNum) {
        $largestNum = $regionNumber;
      }

      my $len=$hash{$key}{"MatureLength"};
      my $type="";
  
      my $SigCov =getCoverage($tRNAMatureBam,$key,$region,$regionNumber); 

      my ($tRNA_temp,$position_temp)=split (/:/,$SigCov);
      my ($position_start,$position_end)=split (/-/,$position_temp);

      ### The start and end positions of tRNA anti codon were obtained
      my ($steamstart,$steamend)=split (/-/,$tRNAposition{$key});

      if (($position_start<2) || ($len-$position_end<2)) {

        if(($steamend+3-$position_end>0) && ($steamstart-$position_end<3)){
          $type="tiRNA-5";
          $region=$SigCov;
          $hash{$key}{$type}[0]=sprintf "%0.2f",$tsRNA_sum*1000000/$MatureMapped;
          $hash{$key}{$type}[1]=$region;
          $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
        }
        elsif(($steamend+3>$position_start) && ($steamstart-3<$position_start)){
          $type="tiRNA-3";
          $region=$SigCov;
          $hash{$key}{$type}[0]=sprintf "%0.2f",$tsRNA_sum*1000000/$MatureMapped;
          $hash{$key}{$type}[1]=$region;
          $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
        }
        elsif($position_start<2){
          $type="tRF-5";
          $region=$SigCov;
          $hash{$key}{$type}[0]=sprintf "%0.2f",$tsRNA_sum*1000000/$MatureMapped;
          $hash{$key}{$type}[1]=$region;
          $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
        }
        elsif($len-$position_end<2) {
          $type="tRF-3";
          $region=$SigCov;
          $hash{$key}{$type}[0]=sprintf "%0.2f",$tsRNA_sum*1000000/$MatureMapped;
          $hash{$key}{$type}[1]=$region;
          $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
        }
      }
      else{
        $type ="tRF-i";
        $region=$SigCov;
        $hash{$key}{$type}[0]=sprintf "%0.2f",$tsRNA_sum*1000000/$MatureMapped;
        $hash{$key}{$type}[1]=$region;
        $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
      }
    }
  }
}

sub handlePri{
  my ($pMin,$PriMapped)=@_;
  my %cluster=();

  open (IN,$tRNAPRIbedGraph);
  while (<IN>) {
    chomp $_;
    my @cols = split(/\s+/,$_);
    my $Tname=$cols[0];
    my $Tlength=$tRNAPriLen{$Tname};

    my $start = $cols[1];
    my $end = $cols[2];
    my $Snum=$cols[3];

    $hash{$Tname}{"Pritotal"}+=$Snum;    #The total expression of reads on Pri-tRNA  was recorded
    for(my $i=$start; $i<$end; $i++){
      $hash{$Tname}{"Pri"}{$i}+=$Snum;
    }
  }
  close IN;

  foreach my $key (keys %hash){
    my $Tlength = $tRNAPriLen{$key};
    
    my $Total = $hash{$key}{"Pritotal"};
    %cluster=();
    for (my $i=0; $i< $Tlength; $i++){
      if (defined $hash{$key}{"Pri"}{$i}) {
        my $Snum=$hash{$key}{"Pri"}{$i};
        my $pvalue=Math::CDF::pbinom($Total-$Snum,$Total,1- 1/$Tlength);

        if ($pvalue <$pMin) {
          unless (defined $cluster{"start"}) {
            $cluster{"start"}=$i;
            $cluster{"end"}=$i;
            next;
          }
          if ($i-$cluster{"end"} == "1") {
            $cluster{"end"}=$i;
          }
        }
        else {
          unless (defined $cluster{"start"}) {
            %cluster=();
            next;
          }
          my $len=$tRNAPriLen{$key};
          my $type="";

          my $region = $key.":".$cluster{"start"}."-".$cluster{"end"};
          #print $key,"\t",$len,"\t",$cluster{"start"},"\t",$cluster{"end"},"\n";
          #my $cmd = "  samtools view -@ 96 ".$tRNAPRIBam.".bam ".$region."| awk \'BEGIN{FS=\"x\";} {total+=\$2; } END {print total;}\'";
          
          my $regionNumber = 0;
          my $regionScore = 0;
          for (my $i = $cluster{"start"}; $i <= $cluster{"end"}; $i++) {
            $regionNumber+=$hash{$key}{"Pri"}{$i};
            # $regionScore+=$hash{$key}{"PriScore"}{$i};
          }
              
          my $SigCov=getCoverage($tRNAPRIBam,$key,$region,$regionNumber);
          if ($SigCov) {
            $region=$SigCov;
            $region =~ /.+\:(.+)\-.+/;
            my $RegionStart =$1;

            if ( $RegionStart>= ($len-50)  ){
                $type="tRF-1";
            }

            $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$PriMapped;
            $hash{$key}{$type}[1]=$region;
            $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
          }
          %cluster=();
        }
      }
    }

    if (defined $cluster{"start"}) {
      my $len=$tRNAPriLen{$key};
      my $type="";
      my $region = $key.":".$cluster{"start"}."-".$cluster{"end"};
      #my $cmd = "  samtools view -@ 96 ".$tRNAPRIBam.".bam ".$region."| awk \'BEGIN{FS=\"x\";} {total+=\$2; } END {print total;}\'";
      my $regionNumber = 0;
      my $regionScore = 0;

      for (my $i = $cluster{"start"}; $i <= $cluster{"end"}; $i++) {
        $regionNumber+=$hash{$key}{"Pri"}{$i};
      }

      my $SigCov=getCoverage($tRNAPRIBam,$key,$region,$regionNumber); 
      if ($SigCov) {
        $region=$SigCov;
        $region =~ /.+\:(.+)\-.+/;
        my $RegionStart =$1;

        if ( $RegionStart>= ($len-50)  ){
          $type="tRF-1";
        }
        $hash{$key}{$type}[0]=sprintf "%0.2f",$regionNumber*1000000/$PriMapped;
        $hash{$key}{$type}[1]=$region;
        $hash{$key}{$type}[2]=sprintf "%0.2f",$regionScore;
      }
    }
  }
}

sub BedtoBedGraph{
  my ($inputFile,$ChromeSizeFile,$outputFile,$logFileHandle) =@_;
  my $command = "sortBed -i ".$inputFile." | awk \'BEGIN{OFS=\"\\t\"}{split(\$4,a,\"x\"); while (a[2]>0) {print \$0,a[2]--} }' | genomeCoverageBed -bg -split -i stdin -g  ".$ChromeSizeFile." > ".$outputFile;

  print $logFileHandle $command,"\n";
  system($command);
}

sub getCoverage{
  my ($bamFile,$key,$region,$regionNumber) =@_;
  my $cmd="samtools view ".$bamFile." ".$region." |  awk 'BEGIN{max=0; start=0; end=0;}{split(\$1,a,\"x\"); len=length(\$10);   if(a[2]*len > max) {max=a[2]*len; start=\$4-1; end=\$4-1+len-1;} } END {print start\"-\"end\"_\"max;}\'";

  my $MaxResult=`$cmd`;
  chomp $MaxResult;

  my ($Maxregion,$MaxReadsNum)=split(/\_/,$MaxResult);
  my $coverage=$MaxReadsNum/($regionNumber+1);

  if ($coverage >= 0) {
    return $key.":".$Maxregion; 
  }
  else {return 0;}
}

sub printHTML{
  my ($tsRFfinderOutFile) = @_;

  open (OUT,">$tsRFfinderOutFile") || die "can not open output file $tsRFfinderOutFile";

  my $htmlstrmiddle = "";
  my $position; my $length;
  
  my %tsRNA_TCGA;
  open (IN,$tsRNA_TCGA) || die "can not open $tsRNA_TCGA";

  while (<IN>) {
    chomp $_;
    my ($k,$v)=split (/\t/,$_);
    $tsRNA_TCGA{$k}=$v;#save to hash
  }
  close IN;
  my $Score = 0;

  foreach my $key (keys %hash){
    if (defined $hash{$key}{"tRF-5"}[0]) {
      $position=$hash{$key}{"tRF-5"}[1];
      my ($tRNA,$sites)=split(/:/,$position);
      my ($start,$end)=split(/-/,$sites);
      $length=$end-$start+1;
      my $tsRNA_name=$tsRNA_TCGA{$key};

      if(!$tsRNA_name) {
        $tsRNA_name = "NA";
      }

      my $RPM = $hash{$key}{"tRF-5"}[0];

      my $seq=uc(substr($tRNAmatureSeq{$key},$start,$length));

      print OUT "tRF-5\t".$RPM."\t".$Score."\t".$length."\t$tRNA\t$sites\t$seq\n";
      
    }

    if (defined $hash{$key}{"tRF-3"}[0]) {
      $position=$hash{$key}{"tRF-3"}[1];
      my ($tRNA,$sites)=split(/:/,$position);
      my ($start,$end)=split(/-/,$sites);
      $length=$end-$start+1;
      
      my $tsRNA_name=$tsRNA_TCGA{$key};

      if(!$tsRNA_name) {
        $tsRNA_name = "NA";
      }

      my $RPM = $hash{$key}{"tRF-3"}[0];
      my $seq=uc(substr($tRNAmatureSeq{$key},$start,$length));
      print OUT "tRF-3\t".$RPM."\t".$Score."\t".$length."\t$tRNA\t$sites\t$seq\n";
    }

    if (defined $hash{$key}{"tiRNA-5"}[0]) {
      $position=$hash{$key}{"tiRNA-5"}[1];
      my ($tRNA,$sites)=split(/:/,$position);
      my ($start,$end)=split(/-/,$sites);
      $length=$end-$start+1;

      my $tsRNA_name=$tsRNA_TCGA{$key};

      if(!$tsRNA_name) {
        $tsRNA_name = "NA";
      }

      my $RPM = $hash{$key}{"tiRNA-5"}[0];

      my $seq=uc(substr($tRNAmatureSeq{$key},$start,$length));
      $htmlstrmiddle.=" <tr>\n<td>tiRNA-5</td>\n<td>$Score</td>\n<td>$length</td>\n<td>$tRNA</td>\n<td>$sites</td>\n<td>$RPM</td>\n<td><button>View</button></td>\n<td><button>View</button></td>\n<td>$tRNA"."_$sites</td>\n<td>$seq</td>\n</tr>\n";    
      print OUT "tiRNA-5\t".$hash{$key}{"tiRNA-5"}[0]."\t".$Score."\t".$length."\t$tRNA\t$sites\t$seq\n";
    }

    if (defined $hash{$key}{"tiRNA-3"}[0]) {
      $position=$hash{$key}{"tiRNA-3"}[1];
      my ($tRNA,$sites)=split(/:/,$position);
      my ($start,$end)=split(/-/,$sites);
      $length=$end-$start+1;

      my $Score=sprintf("%.2f",$hash{$key}{"tiRNA-3"}[2]/$length);

      my $tsRNA_name=$tsRNA_TCGA{$key};

      if(!$tsRNA_name) {
        $tsRNA_name = "NA";
      }

      my $RPM = $hash{$key}{"tiRNA-3"}[0];

      my $seq=uc(substr($tRNAmatureSeq{$key},$start,$length));
    
      print OUT "tiRNA-3\t".$hash{$key}{"tiRNA-3"}[0]."\t".$Score."\t".$length."\t$tRNA\t$sites\t$seq\n";
    }

    if ( defined $hash{$key}{"tRF-1"}[0] ) {
      $position=$hash{$key}{"tRF-1"}[1];
      my ($tRNA,$sites)=split(/:/,$position);
      my ($start,$end)=split(/-/,$sites);
      $length=$end-$start+1;

      my $tsRNA_name=$tsRNA_TCGA{$key};

      if(!$tsRNA_name) {
        $tsRNA_name = "NA";
      }

      my $RPM = $hash{$key}{"tRF-1"}[0];

      my $LastNuc= substr($tRNAPriSeq{$key},$end,1);  #取得tRF-1的结尾字符，使其保证为U
      if ($LastNuc eq "T") {
        my $seq=uc(substr($tRNAPriSeq{$key},$start,$length));
        print OUT "tRF-1\t".$hash{$key}{"tRF-1"}[0]."\t".$Score."\t".$length."\t$tRNA\t$sites\t$seq\n";
      }    
    }

    if ( defined $hash{$key}{"tRF-i"}[0] ) {
      $position=$hash{$key}{"tRF-i"}[1];
      my ($tRNA,$sites)=split(/:/,$position);
      my ($start,$end)=split(/-/,$sites);
      $length=$end-$start+1;
      
      my $tsRNA_name=$tsRNA_TCGA{$key};

      if(!$tsRNA_name) {
        $tsRNA_name = "NA";
      }
      my $RPM = $hash{$key}{"tRF-i"}[0];
      my $seq=uc(substr($tRNAmatureSeq{$key},$start,$length));
      print OUT "tRF-i\t".$hash{$key}{"tRF-i"}[0]."\t".$Score."\t".$length."\t$tRNA\t$sites\t$seq\n";    
    }
  }
  close OUT;
}
