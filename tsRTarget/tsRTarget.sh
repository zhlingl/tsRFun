#!/bin/bash 

##################################################
# tsRTarget.sh find tsRNA-target interaction from 
# CLIP sequencing data
# Lingling Zheng(zhengll33@mail.sysu.edu.cn)
# Data: 2021-9-13
##################################################

###################################
#                                 #
#   tsRNA_target_in_CLIP data     #
#                                 #
#     last change:9/9/2021        #
#                                 #
################################### 

usage="Description:	bash script used to find tsRNA-target interaction from CLIP sequencing data.       
		Query input files support sequence formats: .fastq/.fq, .fasta/.fa.                            
Usage Examples: bash tsRTarget.sh -i fasta/fastq -t clear -o outputdir                                 
Options:                                                                                               
    -i <inputfile>	Input could be:                                                                    
		  a .fastq/.fq or .fasta/.fa file.                                                             
    -t data type clear/clash or CLIP (iclip/eclip/par-clip/hits-clip)                                  
    -o outputdir  address of annotation results                                                        
    -a index_adress                                                                                    
Alignment:                                                                                             
    -m mismatch  default=1                                                                             
    -w min tsRNA_length word_size default=14                                                           
    -n max tsRNA_length word_size default=40                                                           
    -g min target_length  default=10                                                                   
    -h max -g min target_length default(70 for clear clash data; 140 for clip data)                    
    -s matched-length  default=6                                                                       
    -e gapopen	default=1                                                                              
    -P Collapse PCR duplicates default=T                                                               
    -S seed T/N  default=T                                                                             
    -R seed start default=2                                                                            
    -E seed end default=8                                                                              
    -M Min Free Energy  default=-10                                                                    
    -c match_clip the word_size of RNA-RNA target; only used for CLIP seq data,                        
Others:                                                                                                
    -h		print this usage message                                                                   
Example of use:                                                                                        
    bash tsRTarget.sh -i fa -t clear -o outputdir -a index_dir                                         
    bash tsRTarget.sh -i example.fa -t clear -o example_result -a hg38_index                           
"

## Set default parameters
target_length=20
target_maxlength=70
mismatch=1
tsRNA_length=14
tsRNA_maxlength=40
match_l=6
gapopen=0
seed="T"
seeds=2
seede=7
MFE=-10
pcr="T"
match_clip=6
help=""

## Extract the input parameters from the input of users
while getopts i:t:o:a:m:w:n:g:G:s:e:P:S:R:E:M:c:? opts; do
    case $opts in

        i) input=$OPTARG ;;
        t) type=$OPTARG ;;
        o) output=$OPTARG ;;
        a) index_adress=$OPTARG ;;
        m) mismatch=$OPTARG ;;
        w) tsRNA_length=$OPTARG ;;
        n) tsRNA_maxlength=$OPTARG ;;
        g) target_length=$OPTARG ;;
		G) target_maxlength=$OPTARG ;;
        s) match_l=$OPTARG ;;
        e) gapopen=$OPTARG ;;
		P) pcr=$OPTARG ;;
		S) seed=$OPTARG ;;
		R) seedstart=$OPTARG ;;
		E) seedend=$OPTARG ;;
		M) MFE=$OPTARG ;;
		c) match_clip=$OPTARG ;;
        h) help=$OPTARG ;;

        ?) ;;
    esac
done

sample=$(basename $input)
finame=${sample%%.*}

# ctk="/public/home/wangzh/bin/ctk"
# RNAhybrid="/public/home/wangzh/soft/RNAhybrid/src/RNAhybrid"

##  set the detail index
if [[ ! $input ]];
then	
    echo "$usage"
    echo "Please enter the right inputFile"
    exit
fi

if [[ ! $type ]];
then	
    echo "$usage"
    echo "Please enter the right CLIP type"
    exit
fi

if [[ ! $output ]];
then	
    echo "$usage"
    echo "Please enter the right outputDir"
    exit
fi

if [[ ! $index_adress ]];
then	
    echo "$usage"
    echo "Please enter the right index_adress"
    exit
fi

tRNA_bed=$index_adress"/tRNA/tRNA_exons.bed"
tsRNA_db=$index_adress"/tRNA/blast/tRNA"
tsRNA_index=$index_adress"/tRNA/bowtie/mature"
miRNA_bed=$index_adress"/miRNA/miRNA.bed"
miRNA_db=$index_adress"/miRNA/blast/miRNA_blast"
miRNA_index=$index_adress"/miRNA/bowtie/mature"
ENST_bed=$index_adress"/ENST_bed.bed"
name=$index_adress"/tsRNA_TCGA_for_NAR.txt"

genome_db="/public/home/genome/hg38/ucsc/bowtie_index/hg38"
genome_fa="/public/home/genome/hg38/ucsc/hg38.fa"
ENST="/public/home/genome/hg38/TCGA_tRF2Cancer/enst_blast/enst"


if [[ $help != ""* ]];
then	
    echo $usage
    exit
fi


type(){

    file=$input
    line1=$(head -n 1 $file)
    file2="$output/collapsed.fa"

    if [[ $line1 == "@"* ]] || [[ $line1 == ">"* ]];
    then	
        if [[ $line1 == ">seq"* ]];
        then
            cp $file $file2;
        else
        cmd=$(awk 'BEGIN{flag=0;}{if($0 ~ /^@/) {flag=1;k=1;} else {if(flag) {hash[$0]++; flag=0;}} } END{for (i in hash)  {print ">seq_"k"_x"hash[i]"\n"i; k++;}}' $file >  $file2)
        fi
    fi

    echo "    inputFile: $input;
    outputDir: $output;
    dataType: $type"

    if [[ $type == 'clear' ]] || [[ $type == 'clash' ]];

    then

        #####         0、process the input              #######"

        cmd=$(awk '{if($0 ~ /^>/) {printf($0"\t")} else {print $0} }' $file2 | awk -F '[_\t]' '{if(length($4)-"'$tsRNA_length'">-1){$4=toupper($4);print $1"_"$2"_"$4"_"$3"\n"$4}}' >$output/trimed.fa)
	                                                                                  
        miRNA_target(){

            if [[ ! -d $output/miRNA/ ]]; then
                mkdir $output/miRNA/
            fi

            echo "################        seach for miRNA-target interactions          ##############"

            ############  Map reads miRNA ref, and the map length is more than 14nt 

            cmd=$(blastn -query $output/trimed.fa  -out $output/miRNA/$finame.1_miRNAmapped.blast -db  $miRNA_db  \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 16 -word_size $tsRNA_length)
            
            ######### mapped part >= 14 & unmapped part >10 &  mapped part locate at first or last reads  && 16<= reads <= 50 
            cmd=$(awk '{if( ($12 < ("'$tsRNA_maxlength'"+1)) &&  ( $12 > ("'$tsRNA_length'" -1) )  && ( $5 < 5 || ($2 - $6 ) <5 ) && ( ($5 - "'$target_length'" > 0)  || ($2 - "'$target_length'" > $6) ) && ( $14 < ("'$mismatch'"+1) ) && ( $15 < ("'$gapopen'"+1) ) ) \
                print $1 }' $output/miRNA/${finame}.1_miRNAmapped.blast |uniq |  awk -F '_' '{print ">"$0"\n"$3}' >$output/miRNA/${finame}.2_miRNA.for_tsRNA.blast)


            ##### revome the reads mapped to tsRNA ref
            ######        3、revome the reads mapped to tsRNA ref       ######" 
            cmd=$(blastn -query $output/miRNA/${finame}.2_miRNA.for_tsRNA.blast  -out $output/miRNA/${finame}.3_miRNA.tsRNA_candidates.blast -db  $tsRNA_db  \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 16 -word_size $tsRNA_length)

            ######  去除 miRNA的潜在序列 获得只能map到tRNA上的reads  ####
            cmd=$(awk -F ' ' '{print $1}'  $output/miRNA/${finame}.3_miRNA.tsRNA_candidates.blast |  awk NF | sort | uniq >$output/miRNA/${finame}.test1.txt;)
            cmd=$(awk -F '>' '{print $2}' $output/miRNA/${finame}.2_miRNA.for_tsRNA.blast | awk NF|sort| uniq > $output/miRNA/${finame}.test2.txt;)

            if	[[ -s $output/miRNA/${finame}.test1.txt ]] && [[ -s $output/miRNA/${finame}.test2.txt ]]
            then
                cmd=$(awk 'NR==FNR{a[$0]}NR>FNR{ if(($0 in a)) print $0}' $output/miRNA/${finame}.test1.txt $output/miRNA/${finame}.test2.txt | awk -F '_' '{print ">"$0"\n"$3}'>$output/miRNA/${finame}.3_miRNA.Reads_unalignedtsRNA.fa)

            else
                grep -vwf $output/miRNA/${finame}.test1.txt $output/miRNA/${finame}.test2.txt | awk -F '_' '{print ">"$0"\n"$3}'>$output/miRNA/${finame}.3_miRNA.Reads_unalignedtsRNA.fa
            fi


            ######    4、revome the reads (fake reads) mapped to genome ref   #######" 
            cmd=$(bowtie -f -a -m 20 -v $mismatch --un $output/miRNA/${finame}.4_miRNA.Genome_unaligned.fa $genome_db  \
                    $output/miRNA/${finame}.3_miRNA.Reads_unalignedtsRNA.fa  $output/miRNA/${finame}.4_miRNA.Genome.bwt )

            cmd=$(blastn -query $output/miRNA/${finame}.4_miRNA.Genome_unaligned.fa -out $output/miRNA/${finame}.4_miRNA.miRNA_mapped.blast \
                -db $miRNA_db -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 16 -word_size $tsRNA_length)
            ######    5、select the mti_candidates     #######" 

            cmd=$(awk '{ if( $12 - "'$tsRNA_maxlength'" < 1 &&  ( $12 - "'$tsRNA_length'" +1 > 0 ) && ( $5 < 5 || ($2 - $6 ) <5 ) && (( $5 - "'$target_length'" > 0) || ( $2-$6 - "'$target_length'" > 0 )) && ($14 < ("'$mismatch'"+1)) && ($15 < ("'$gapopen'"+1))  )  print }'  $output/miRNA/${finame}.4_miRNA.miRNA_mapped.blast 	>$output/miRNA/${finame}.5_mti_candidates.blast)

            ########## control 8nt unmapped   mismatch<2  gapopen=0		
            ######    6、select the miRNA and targets of mti_candidates      ######"
            
            cmd=$(awk -F '[ _\t]' '{if($8>4){print $1"_"$2"#"$6"#"$12"_"$10"-"$11"_"$8"_"$9,$12,$1"_"$2"_"$4"_"substr($3,1,$8-1)"_1_"$8-1,substr($3,1,$8-1)} \
                else { print $1"_"$2"#"$6"#"$12"_"$10"-"$11"_"$8"_"$9,$12,$1"_"$2"_"$4"_"substr($3,$9+1,$2)"_"$9+1"_"$5,substr($3,$9+1,$2)}}' \
                $output/miRNA/${finame}.5_mti_candidates.blast | awk  '{gsub("-","",$2);print}' | sort  >$output/miRNA/${finame}.6_mti_candidates2.fa)
                    
            ######### 将gapopen的 “-” 替换为 空  并提取tsRNA和targets序列

            cmd=$(awk -F '#' 'BEGIN{tRNA="";} {if($3!=tRNA){hash[$3]=$1"_"$2; tRNA=$3}  else {hash[$3]=hash[$3]"#"$2 ;}  } \
                END{for (i in hash) {print hash[i]"_"i}}' $output/miRNA/${finame}.6_mti_candidates2.fa | awk '{if(length($4)<"'$target_maxlength'") print }' >$output/miRNA/${finame}.6_mti_brife.fa)

            cmd=$(awk -F ' ' '{print ">"$3,$4}' $output/miRNA/${finame}.6_mti_brife.fa  | awk -F ' ' '{print $1"_\n"$2}' \
                > $output/miRNA/${finame}.6_targets_candidates.fa)

            cmd=$(awk -F ' ' '{print ">"$1"\n"$2}' $output/miRNA/${finame}.6_mti_brife.fa >$output/miRNA/${finame}.6_miRNA_candidates.fa)
            
            ####  RNAhybrid search for tsRNA and target

            if [[ $seed == "T" ]] || [[ $seed == "t" ]];
            then
                
                if [[  -f $output/miRNA/${finame}.6.RNAhybrid_brief_final.txt ]]
                then
                    cmd=$( rm $output/miRNA/${finame}.6.RNAhybrid_brief_final.txt )
                fi
                
                cmd=$(cat $output/miRNA/${finame}.6_mti_brife.fa | while read reads; do  echo $reads > $output/miRNA/temp.read.fa ; awk -F' ' \
                    '{print ">"$1"_\n"$2}' $output/miRNA/temp.read.fa >$output/miRNA/temp.miRNA.fa ; awk -F' ' '{print ">"$3"_\n"$4}' $output/miRNA/temp.read.fa >$output/miRNA/temp.target.fa ;\
                    RNAhybrid -c -b 1 -u 2 -v 2 -f $seeds,$seede -n $tsRNA_maxlength -e $MFE -m $target_maxlength -s 3utr_human -t $output/miRNA/temp.target.fa -q $output/miRNA/temp.miRNA.fa >>$output/miRNA/${finame}.6.RNAhybrid_brief_final.txt ; \
                    done)

            fi	

            ########### blast search for tsRNA and target
            
            ######     7、select the miRNA and targets that can partly match        #######
            cmd=$(makeblastdb -in $output/miRNA/${finame}.6_targets_candidates.fa  -dbtype nucl -out $output/miRNA/blast_miRNAindex/miRNA)
            
            cmd=$(blastn -query $output/miRNA/${finame}.6_miRNA_candidates.fa  -out $output/miRNA/${finame}.7_miRNA.target_map.blast \
                -db $output/miRNA/blast_miRNAindex/miRNA \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand minus -num_threads 16 -word_size $match_l)

            cmd=$(awk -F '[\t_]' '{if( $2==$10 && (($14-$6 < 5 && $14-$6>0)||($5-$15 < 5 && $5-$15>0))) {print }}' \
                $output/miRNA/${finame}.7_miRNA.target_map.blast >$output/miRNA/${finame}.7_mti_candidates_final.txt)

            ######     8、select the miRNA and targets that can partly match        #######" 

            cmd=$(awk '{print $3}' $output/miRNA/${finame}.7_mti_candidates_final.txt | sort | uniq | awk -F '[_]' '{print ">unseed:"$1"_"$2"_"$3"\n"$4}'>$output/miRNA/${finame}.8.target.fa)
            
            cmd=$(awk -F ':' '{print $1}' $output/miRNA/${finame}.6.RNAhybrid_brief_final.txt | sort | uniq | awk -F '_' '{print ">seed:"$1"_"$2"_"$3"\n"$4}' >>$output/miRNA/${finame}.8.target.fa)

            ## ENST
            cmd=$(blastn -query $output/miRNA/${finame}.8.target.fa -out $output/miRNA/${finame}.8.target_ENSG.txt -db $ENST \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 26 -word_size $target_length )

            cmd=$(awk -F'[:\t]' '{if($3-$16 < 5 )print $1,$2,$4,$6,$7,$9,$10,$11,$12}' $output/miRNA/${finame}.8.target_ENSG.txt |  sed 's/(.*)//'| awk -F '[ -]' '{print $1,$2,$3,$4,$5+$9-1,$5+$10,$7,$8}' | sort -k 2 |uniq >$output/miRNA/${finame}.8.target_ENSG_final.txt)
                
            cmd=$(awk -F '[_:]' '{print $1"_"$2"_"$3,":target "$4,$5,$6,":tsRNA "$11"_"$13"|"$11"_"$13"|"$12,$14"-"$15,":MFE "$18,":p-value "$19,":match at target "$20,":target_boundary:"$21":"$22":"$23":"$24}' \
                $output/miRNA/${finame}.6.RNAhybrid_brief_final.txt | sort |	awk -F '|' 'BEGIN{while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$2);$2=map[$2]; if($2~/^$/) {$2="NA";print}else{print} }' | \
                awk -F '[:]' '{printf $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9":"$10":"$11" ";gsub(" ","",$8);gsub(" ","",$9);printf ":match_target length "length($8)+length($9)"\n"}' >$output/miRNA/${finame}.8.seed_final.txt)
            

            cmd=$(awk -F '[_\t]' '{rev="rev|sed 's/T/U/g'";sed="sed 's/T/U/g'"; printf $9"_"$10"_"$11" :target "$12" "$13" "$14" :miRNA "$3"_"$5" "$3"_"$5" "$4" "$6"-"$7" :evalue "$23":mismatch "$26":gapopen "$27; \
                $12=substr($12,$20,$19-$20+1);printf $12|&rev;close(rev,"to");rev|& getline out;close(rev); printf substr($4,$17,$18-$17+1)|& sed;close(sed,"to"); sed|& getline out1;close(sed); \
                printf " :target_match "$17"\t"$18":\t\t"out1":miRNA_match "$19"\t\t"$20":"out"\n" }' $output/miRNA/${finame}.7_mti_candidates_final.txt | awk -F ':' '{print $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":\t\t"$10":"$9}' | \
                sort |	awk 'BEGIN{FS=OFS=" ";while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$8);$8= map[$8]; if($8~/^$/) {$8="NA";print}else{print} }' >$output/miRNA/${finame}.8.unseed_final.txt)
        
            cmd=$(join -t $' t' -12 -21  $output/miRNA/${finame}.8.target_ENSG_final.txt  $output/miRNA/${finame}.8.seed_final.txt | sort |uniq | awk '{if($2=="seed") print}'> $output/miRNA/${finame}.miRNA_match_final.txt)

            cmd=$(join -t $' t' -12 -21 $output/miRNA/${finame}.8.target_ENSG_final.txt  $output/miRNA/${finame}.8.unseed_final.txt | sort | uniq | awk '{if($2=="unseed") print}'  >> $output/miRNA/${finame}.miRNA_match_final.txt)
            
            cmd=$(sort -t $' '  -k 1 -k 2 $output/miRNA/${finame}.miRNA_match_final.txt > $output/miRNA_${finame}.match_final.txt)

        }
        

        tsRNA_target(){

            if [[ ! -d $output/tsRNA/ ]]; then
                mkdir $output/tsRNA/
            fi

            echo "################        seach for miRNA-target interactions          ##############"

            ######将count的reads map to tsRNA_db ref，要求14nt的比对范围
            cmd=$(blastn -query $output/trimed.fa  -out $output/tsRNA/$finame.1_tsRNAmapped.blast -db  $tsRNA_db  \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 26 -word_size $tsRNA_length)
                
            cmd=$(awk '{if( $12 < "'$tsRNA_maxlength'"+1 &&  ( $12 > ("'$tsRNA_length'" -1) )  && ( $5 < 5 || ($2 - $6 ) <5 ) && ( ($5 - "'$target_length'" > 0)  || ($2 - "'$target_length'" > $6) ) && ( $14 < ("'$mismatch'"+1) ) && ( $15 < ("'$gapopen'"+1) ) ) \
                print $1 }' $output/tsRNA/${finame}.1_tsRNAmapped.blast  | uniq | awk -F '_' '{print ">"$0"\n"$3}' >$output/tsRNA/${finame}.2_tsRNA.for_miRNA.blast )

            #### revome the reads mapped to miRNA ref
            cmd=$(blastn -query $output/tsRNA/${finame}.2_tsRNA.for_miRNA.blast  -out $output/tsRNA/${finame}.3_tsRNA.miRNA_candidates.blast -db  $miRNA_db  \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 26 -word_size $tsRNA_length)

            # 将count的reads map to miRNA ref，

            ##  去除 miRNA的潜在序列 获得只能map到tRNA上的reads
            ######  去除 tsRNA的潜在序列 获得只能map到tRNA上的reads  ####
            cmd=$(awk -F ' ' '{print $1} '  $output/tsRNA/${finame}.3_tsRNA.miRNA_candidates.blast |  awk NF | sort | uniq >$output/tsRNA/${finame}.test1.txt;)
            cmd=$(awk -F '>' '{print $2} ' $output/tsRNA/${finame}.2_tsRNA.for_miRNA.blast | awk NF|sort| uniq > $output/tsRNA/${finame}.test2.txt;)

            if	[[ -s $output/tsRNA/${finame}.test1.txt ]] && [[ -s $output/tsRNA/${finame}.test2.txt ]]
            then
                cmd=$(awk  'NR==FNR{a[$0]}NR>FNR{ if(($0 in a)) print $0}' $output/tsRNA/${finame}.test1.txt $output/tsRNA/${finame}.test2.txt | awk -F '_' '{print ">"$0"\n"$3}'>$output/tsRNA/${finame}.3_tsRNA.Reads_unalignedtsRNA.fa)

            else

                cmd=$(grep -vwf $output/tsRNA/${finame}.test1.txt $output/tsRNA/${finame}.test2.txt | awk -F '_' '{print ">"$0"\n"$3}'>$output/tsRNA/${finame}.3_tsRNA.Reads_unalignedtsRNA.fa)

            fi

            ########  将reads进一步map到genome  去除能匹配到genome的 fake reads	

            cmd=$(bowtie -f -a -m 20 -v $mismatch --un $output/tsRNA/${finame}.4_tsRNA.Genome_unaligned.fa $genome_db  \
                $output/tsRNA/${finame}.3_tsRNA.Reads_unalignedtsRNA.fa  $output/tsRNA/${finame}.4_tsRNA.Genome.bwt )

            ######将候选的tsRNA序列再次比对到tsRNA ref 获得嵌合体组合

            cmd=$(blastn -query $output/tsRNA/${finame}.4_tsRNA.Genome_unaligned.fa -out $output/tsRNA/${finame}.5_tsRNA.tsRNA_mapped.blast \
                -db $tsRNA_db -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 26 -word_size $tsRNA_length)

            cmd=$(awk '{ if( $12 - "'$tsRNA_maxlength'" < 1 &&  ( $12 - "'$tsRNA_length'" +1 > 0 ) && ( $5 < 5 || ($2 - $6 ) <5 ) && (( $5 - "'$target_length'" > 0) || ( $2-$6 - "'$target_length'" > 0 )) && ($14 < ("'$mismatch'"+1)) && ($15 < ("'$gapopen'"+1))  )  print }' $output/tsRNA/${finame}.5_tsRNA.tsRNA_mapped.blast 	>$output/tsRNA/${finame}.5_tti_candidates.blast)
            
            ##### control 8nt unmapped 片段在tRNA起始或终止区域  mismatch<2  gapopen<2
            ######    6、select the tsRNA and targets of tti_candidates      ######"

            cmd=$(awk -F '[ _\t]' '{if($8>4){print $1"_"$2"#"$6"#"$12"_"$10"-"$11"_"$8"_"$9,$12,$1"_"$2"_"$4"_"substr($3,1,$8-1)"_1_"$8-1,substr($3,1,$8-1)} \
                else { print $1"_"$2"#"$6"#"$12"_"$10"-"$11"_"$8"_"$9,$12,$1"_"$2"_"$4"_"substr($3,$9+1,$2)"_"$9+1"_"$5,substr($3,$9+1,$2)}}' \
                $output/tsRNA/${finame}.5_tti_candidates.blast | awk  '{gsub("-","",$2);print}' | sort  >$output/tsRNA/${finame}.6_tti_candidates2.fa)
            ########将gapopen的 “-” 替换为 空  并提取tsRNA和targets序列

            cmd=$(awk -F '#' 'BEGIN{tRNA="";} {if($3!=tRNA){hash[$3]=$1"_"$2; tRNA=$3}  else {hash[$3]=hash[$3]"#"$2 ;}  } \
            END{for (i in hash) {print hash[i]"_"i}}' $output/tsRNA/${finame}.6_tti_candidates2.fa | awk '{if(length($4)<"'$target_maxlength'") print }' >$output/tsRNA/${finame}.6_tti_brife.fa)

            cmd=$(awk -F ' ' '{print ">"$3,$4}' $output/tsRNA/${finame}.6_tti_brife.fa  | awk -F ' ' '{print $1"_\n"$2}' \
            > $output/tsRNA/${finame}.6_targets_candidates.fa)

            cmd=$(awk -F ' ' '{print ">"$1"\n"$2}' $output/tsRNA/${finame}.6_tti_brife.fa >$output/tsRNA/${finame}.6_tsRNA_candidates.fa)
        
            ############# RNAhybrid 种子区比对tsRNA and target
            if [[ $seed == "T" ]] || [[ $seed == "t" ]];
            then

                if [[  -f $output/tsRNA/${finame}.6.RNAhybrid_brief_final.txt ]]; then
                    cmd=$(rm $output/tsRNA/${finame}.6.RNAhybrid_brief_final.txt)
                fi

                cmd=$(cat $output/tsRNA/${finame}.6_tti_brife.fa | while read reads; do  echo $reads > $output/tsRNA/temp.read.fa ; awk -F' ' \
                    '{print ">"$1"_\n"$2}' $output/tsRNA/temp.read.fa >$output/tsRNA/temp.tsRNA.fa ; awk -F' ' '{print ">"$3"_\n"$4}' $output/tsRNA/temp.read.fa >$output/tsRNA/temp.target.fa ;\
                    RNAhybrid -c -b 1 -u 2 -v 2 -f $seeds,$seede -n $tsRNA_maxlength -e $MFE -m $target_maxlength -s 3utr_human -t $output/tsRNA/temp.target.fa -q $output/tsRNA/temp.tsRNA.fa >>$output/tsRNA/${finame}.6.RNAhybrid_brief_final.txt ; \
                    done)

            fi	


            ######     7、select the tsRNA and targets that can partly match        #######

            cmd=$(makeblastdb -in $output/tsRNA/${finame}.6_targets_candidates.fa  -dbtype nucl -out $output/tsRNA/blast_tsRNAindex/tsRNA)
                
            cmd=$(blastn -query $output/tsRNA/${finame}.6_tsRNA_candidates.fa  -out $output/tsRNA/${finame}.7_tsRNA.target_map.blast \
                -db $output/tsRNA/blast_tsRNAindex/tsRNA \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand minus -num_threads 26 -word_size $match_l)

            cmd=$(awk -F '[\t_]' '{if( $2==$10 && (($14-$6 < 5 && $14-$6>0)||($5-$15 < 5 && $5-$15>0))) {print }}' $output/tsRNA/${finame}.7_tsRNA.target_map.blast \
                >$output/tsRNA/${finame}.7_tti_candidates_final.txt)

            ######     8、select the tsRNA and targets that can partly match        #######

            cmd=$(awk '{print $3}' $output/tsRNA/${finame}.7_tti_candidates_final.txt | sort | uniq | awk -F '[_]' '{print ">unseed:"$1"_"$2"_"$3"\n"$4}'>$output/tsRNA/${finame}.8.target.fa)
            
            cmd=$(awk -F ':' '{print $1}' $output/tsRNA/${finame}.6.RNAhybrid_brief_final.txt | sort | uniq | awk -F '_' '{print ">seed:"$1"_"$2"_"$3"\n"$4}' >>$output/tsRNA/${finame}.8.target.fa)
            
            ## ENST
            cmd=$(blastn -query $output/tsRNA/${finame}.8.target.fa -out $output/tsRNA/${finame}.8.target_ENSG.txt -db $ENST \
                -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
                -strand plus -num_threads 26 -word_size $target_length)

            cmd=$(awk -F'[:\t]' '{if($3-$16 < 5 )print $1,$2,$4,$6,$7,$9,$10,$11,$12}' $output/tsRNA/${finame}.8.target_ENSG.txt | sed 's/(.*)//'| awk -F '[ -]' '{print $1,$2,$3,$4,$5+$9-1,$5+$10,$7,$8}' | sort -k 2 |uniq >$output/tsRNA/${finame}.8.target_ENSG_final.txt)

            cmd=$(awk -F '[_:]' '{print $1"_"$2"_"$3,":target "$4,$5,$6,":tsRNA "$11"_"$13"|"$11"_"$13"|"$12,$14"-"$15,":MFE "$18,":p-value "$19,":match at target "$20,":target_boundary:"$21":"$22":"$23":"$24}' \
                $output/tsRNA/${finame}.6.RNAhybrid_brief_final.txt | sort |	awk -F '|' 'BEGIN{while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$2);$2=map[$2]; if($2~/^$/) {$2="NA";print}else{print} }' | \
                awk -F '[:]' '{printf $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9":"$10":"$11" ";gsub(" ","",$8);gsub(" ","",$9);printf ":match_target length "length($8)+length($9)"\n"}' >$output/tsRNA/${finame}.8.seed_final.txt)

            cmd=$(awk -F '[_\t]' '{rev="rev|sed 's/T/U/g'";sed="sed 's/T/U/g'"; printf $9"_"$10"_"$11" :target "$12" "$13" "$14" :tsRNA "$3"_"$5" "$3"_"$5" "$4" "$6"-"$7" :evalue "$23":mismatch "$26":gapopen "$27; \
                $12=substr($12,$20,$19-$20+1);printf $12|&rev;close(rev,"to");rev|& getline out;close(rev); printf substr($4,$17,$18-$17+1)|& sed;close(sed,"to"); sed|& getline out1;close(sed); \
                printf " :target_match "$17"\t"$18":\t\t"out1":tsRNA_match "$19"\t\t"$20":"out"\n" }' $output/tsRNA/${finame}.7_tti_candidates_final.txt | awk -F ':' '{print $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":\t\t"$10":"$9}' | \
                sort |	awk 'BEGIN{FS=OFS=" ";while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$8);$8= map[$8]; if($8~/^$/) {$8="NA";print}else{print} }' >$output/tsRNA/${finame}.8.unseed_final.txt)
        
            cmd=$(join -t $' t' -12 -21  $output/tsRNA/${finame}.8.target_ENSG_final.txt  $output/tsRNA/${finame}.8.seed_final.txt | sort |uniq | awk '{if($2=="seed") print}'> $output/tsRNA/${finame}.tsRNA_match_final.txt)

            cmd=$(join -t $' t' -12 -21 $output/tsRNA/${finame}.8.target_ENSG_final.txt  $output/tsRNA/${finame}.8.unseed_final.txt | sort | uniq | awk '{if($2=="unseed") print}'  >> $output/tsRNA/${finame}.tsRNA_match_final.txt)
            
            cmd=$(sort -t $' '  -k 1 -k 2 $output/tsRNA/${finame}.tsRNA_match_final.txt > $output/tsRNA_${finame}.match_final.txt)

            date
        
        }

        miRNA_target

        tsRNA_target

        target_maxlength=150

        if [  -f $output/miRNA_${finame}.match_final.txt ] && [  -f $output/tsRNA_${finame}.match_final.txt ];
        then

            if [[ ! -d $output/tsRNA/cenet/ ]]; then
                cmd=$(mkdir $output/tsRNA/cenet/)
            fi
                cmd=$(awk '{print $3,$14}' $output/miRNA_${finame}.match_final.txt |sort > $output/tsRNA/cenet/${finame}.1.miRNA_target.txt)
                cmd=$(awk '{print $3,$14}' $output/tsRNA_${finame}.match_final.txt|sort > $output/tsRNA/cenet/${finame}.1.tsRNA_target.txt)
                cmd=$(join -t $' t' -11 -21 $output/tsRNA/cenet/${finame}.1.miRNA_target.txt $output/tsRNA/cenet/${finame}.1.tsRNA_target.txt > $output/tsRNA/cenet/${finame}.2.cenet.txt)

            if [[ -s $output/tsRNA/cenet/${finame}.2.cenet.txt ]];
            then 

                cmd=$(awk '{print $3"_"$14}' $output/miRNA_${finame}.match_final.txt > $output/tsRNA/cenet/${finame}.2.miRNA.txt)
                
                cmd=$(paste -d'\t' $output/tsRNA/cenet/${finame}.2.miRNA.txt $output/miRNA_${finame}.match_final.txt |  sed  's/seq_.*_x[0-9]* /NA /'| sort -k 1|uniq > $output/tsRNA/cenet/${finame}.3.miRNA.txt)

                cmd=$(awk '{print $3"_"$14}' $output/tsRNA_${finame}.match_final.txt > $output/tsRNA/cenet/${finame}.2.tsRNA.txt)
                
                cmd=$(paste -d'\t' $output/tsRNA/cenet/${finame}.2.tsRNA.txt $output/tsRNA_${finame}.match_final.txt | sed  's/seq_.*_x[0-9]* /NA /'| sort -k 1|uniq  > $output/tsRNA/cenet/${finame}.3.tsRNA.txt)

                cmd=$(awk '{$2=$1"_"$2;$3=$1"_"$3;print $2"\t"$3"\t"}' $output/tsRNA/cenet/${finame}.2.cenet.txt | sort -k 1 >$output/tsRNA/cenet/${finame}.3.tr_cenet.txt)

                cmd=$(join   -t $'\t' -11 -21 $output/tsRNA/cenet/${finame}.3.tr_cenet.txt  $output/tsRNA/cenet/${finame}.3.miRNA.txt |  sed 's/$/&|/g'| sed  's/[^\t]*\t/NA\t/' |sort -k 2 | uniq >$output/tsRNA/cenet/${finame}.4.miRNA_cenet.txt)
                
                cmd=$(join  -t $'\t' -12 -21 $output/tsRNA/cenet/${finame}.4.miRNA_cenet.txt  $output/tsRNA/cenet/${finame}.3.tsRNA.txt    |  sed -e 's/[^ ]* /NA /' |sort | uniq  >$output/tsRNA/cenet/${finame}.5.cenet.txt)
                
                cmd=$(cp $output/tsRNA/cenet/${finame}.5.cenet.txt $output/cenetwork_${finame}.match_final.txt)
            else
                echo "No common targets of tsRNA and miRNA was found"
                cmd=$(touch $output/cenetwork_${finame}.match_final.txt)
            fi
            
        else
            echo "No common targets of tsRNA and miRNA was found"
            touch $output/cenetwork_${finame}.match_final.txt

        fi
        
        ###### $process ;

    elif [[ $type == 'iclip' ]] || [[ $type == 'eclip' ]] ||  [[ $type == 'par-clip' ]] || [[ $type == 'hits-clip' ]];

    then

        if [[ ! -d $output/tsRNA/ ]]; then
            mkdir $output/tsRNA/
        fi

        date

        cmd=$(awk '{if($0 ~ /^>/) {printf($0"\t")} else {print $0} }' $file2 | awk -F '[_\t]' '{if(length($4)-"'$tsRNA_length'">-1){$4=toupper($4);print $1"_"$2"_"$4"_"$3"\n"$4}}' >$output/trimed.fa)

        #####         1、map to genome              #######"

        # echo $cmd
        cmd=$(bowtie -p 16  $genome_db -f  $output/trimed.fa -v $mismatch --best -S --un $output/tsRNA/$finame.1.genome_unmap.fa --al $output/tsRNA/$finame.1.genome_map.fa $output/tsRNA/$finame.1.genome_map.sam )


        # ########  Collapse PCR duplicates 
        if [[ $pcr == "T" ]] || [[ $pcr == "t" ]]; 
        then
            
            #####         2、revome pcr Duplication              #######
            cmd=$(samtools view -b -S $output/tsRNA/$finame.1.genome_map.sam >$output/tsRNA/$finame.2.genome_map.bam)
            cmd=$(samtools rmdup -s  $output/tsRNA/$finame.2.genome_map.bam  $output/tsRNA/$finame.2.genome_map.rmdup.bam)
            cmd=$(bedtools bamtobed -i $output/tsRNA/$finame.2.genome_map.rmdup.bam | sort -k 4 >$output/tsRNA/$finame.2.genome_map.rmdup.bed)

        fi

        ######	Peak calling

        cmd=$(perl ctk/tag2peak.pl -big -ss -v --valley-seeking --valley-depth 0.9  $output/tsRNA/$finame.2.genome_map.rmdup.bed  \
            $output/tsRNA/$finame.3.genome_map.rmdup.peak.bed  --out-boundary $output/tsRNA/$finame.3.genome_map.rmdup.peak.boundary.bed  >call_peak.log)

        ## ENST
        cmd=$(bedtools intersect -a $output/tsRNA/$finame.3.genome_map.rmdup.peak.bed -b $ENST_bed -s -wb |awk '{if( ($3- "'$target_length'")>$2) \
            print $1"\t"$2"\t"$3"\t"$4"_"$10}'| sort | uniq > $output/tsRNA/$finame.4.genome_ensg.txt)

        cmd=$(bedtools intersect -a $output/tsRNA/$finame.3.genome_map.rmdup.peak.bed -b $tRNA_bed -s  -wb |awk '{if( ($3- "'$tsRNA_length'")>$2) \
            print $1"\t"$2"\t"$3"\t"$4"_"$10"_"$2-$8"_"$3-$8}'| sort | uniq > $output/tsRNA/$finame.4.genome_tRNA.txt)

        cmd=$(bedtools intersect -a $output/tsRNA/$finame.3.genome_map.rmdup.peak.bed -b $miRNA_bed -s  -wb |awk '{if( ($3- "'$tsRNA_length'")>$2) \
            print $1"\t"$2"\t"$3"\t"$4"_"$10"_"$2-$8"_"$3-$8}'| sort | uniq > $output/tsRNA/$finame.4.genome_miRNA.txt)

        ####  get target sequence
        cmd=$(bedtools getfasta -fi $genome_fa -bed $output/tsRNA/$finame.4.genome_ensg.txt -s -name |   awk -F ':' '{gsub("_","#",$3);gsub("chr","chr.",$3);print }'| awk '{if($0 ~ /^>/) {printf $0"\t"} else {print $0} }' | \
            awk -F'[_\t ]' '{if(length($7)<"'$tsRNA_maxlength'"+1){$7=toupper($7);print ">"$3"_"$7"_"$5"_"$6"\n"$7}}' | sed 's/:*//g' | sed 's/\.//g' | sed 's/(.*)//g' |sed 's/-/_/g' >$output/tsRNA/$finame.4.genome_target.fa )

        ####  get tsRNA sequence
        cmd=$(bedtools getfasta -fi $genome_fa -bed $output/tsRNA/$finame.4.genome_tRNA.txt -s -name | awk -F ':' '{gsub("_","#",$3);gsub("chr","chr.",$3);print }'| awk '{if($0 ~ /^>/) {printf $0"\t"} else {print $0} }' | \
            awk -F'[_\t ]' '{if(length($9)<"'$target_maxlength'"+1){$9=toupper($9);print ">"$3"_"$9"_"$4"_"$5"\n"$9}}' >$output/tsRNA/$finame.4.genome_tRNA.fa )

        ####  get miRNA sequence
        cmd=$(bedtools getfasta -fi $genome_fa -bed $output/tsRNA/$finame.4.genome_miRNA.txt -s -name | awk -F ':' '{gsub("_","#",$3);gsub("chr","chr.",$3);print }'| awk '{if($0 ~ /^>/) {printf $0"\t"} else {print $0} }' | \
            awk -F'[_\t ]' '{if(length($9)<"'$target_maxlength'"+1){$9=toupper($9);print ">"$3"_"$9"_"$4"_"$5"\n"$9}}' >$output/tsRNA/$finame.4.genome_miRNA.fa )


        ##### 2.所有序列直接匹配获得tRNA miRNA片段 target序列仍从peak中获得


        # echo "#####         5、种子区域比对结果              #######
        cmd=$(RNAhybrid -c -b 1 -u 1 -v 1 -f $seeds,$seede -n $tsRNA_maxlength -e $MFE -m $target_maxlength -s 3utr_human -t $output/tsRNA/$finame.4.genome_target.fa \
            -q $output/tsRNA/$finame.4.genome_tRNA.fa >$output/tsRNA/${finame}.5.tRNAhybrid_brief_final.txt)

        #####         5、种子区域比对结果              #######
        cmd=$(RNAhybrid -c -b 1 -u 1 -v 1 -f $seeds,$seede -n $tsRNA_maxlength -e $MFE -m $target_maxlength -s 3utr_human -t $output/tsRNA/$finame.4.genome_target.fa \
            -q $output/tsRNA/$finame.4.genome_miRNA.fa >$output/tsRNA/${finame}.5.miRNAhybrid_brief_final.txt)


        ######     6、非种子区域比对结果        #######"
        cmd=$(makeblastdb -in $output/tsRNA/$finame.4.genome_tRNA.fa  -dbtype nucl -out $output/tsRNA/blast_tsRNAindex/tsRNA)
                
        cmd=$(makeblastdb -in $output/tsRNA/$finame.4.genome_miRNA.fa  -dbtype nucl -out $output/tsRNA/blast_tsRNAindex/miRNA)
        

        cmd=$(blastn -query $output/tsRNA/$finame.4.genome_target.fa  -out $output/tsRNA/${finame}.6.tsRNA.target_map.blast \
            -db $output/tsRNA/blast_tsRNAindex/tsRNA \
            -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
            -strand minus -num_threads 26 -word_size $match_clip)

        cmd=$(blastn -query $output/tsRNA/$finame.4.genome_target.fa  -out $output/tsRNA/${finame}.6.miRNA.target_map.blast \
            -db $output/tsRNA/blast_tsRNAindex/miRNA \
            -outfmt '6 qseqid  qlen  sseqid  slen  qstart  qend  sstart  send  qseq   sseq   evalue   length  pident  mismatch  gapopen' \
            -strand minus -num_threads 26 -word_size $match_clip)

        cmd=$(awk -F '[_:]' '{print $1"_"$7"_"$9"-"$10" seed",$1,$3,$4,$5,"NA NA :target",$2,"NA :tsRNA "$7"_"$9"-"$10"|"$7"_"$9"-"$10"|"$8,"NA :MFE "$12,":p-value "$13,":match at target "$14,":target and tsRNA boundary:"$15":"$16":"$17":"$18}' \
            $output/tsRNA/${finame}.5.tRNAhybrid_brief_final.txt | sort |	awk -F '|' 'BEGIN{while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$2);$2= map[$2]; if($2~/^$/) {$2="NA";print}else{print} }'| \
            awk -F '[:]' '{printf $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9":"$10":"$11" ";gsub(" ","",$8);gsub(" ","",$9);printf ":match_target length "length($8)+length($9)"\n"}' >$output/tsRNA/${finame}.6.tsRNA.seed_final.txt)

        cmd=$(awk -F '[_:]' '{print $1"_"$7"_"$9"-"$10" seed",$1,$3,$4,$5,"NA NA :target",$2,"NA :miRNA "$7"_"$9"-"$10"|"$7"_"$9"-"$10"|"$8,"NA :MFE "$12,":p-value "$13,":match at target "$14,":target and miRNA boundary:"$15":"$16":"$17":"$18}' \
            $output/tsRNA/${finame}.5.miRNAhybrid_brief_final.txt | sort |	awk -F '|' 'BEGIN{while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$2);$2= map[$2]; if($2~/^$/) {$2="NA";print}else{print} }'| \
            awk -F '[:]' '{printf $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":"$9":"$10":"$11" ";gsub(" ","",$8);gsub(" ","",$9);printf ":match_target length "length($8)+length($9)"\n"}' >$output/tsRNA/${finame}.6.miRNA.seed_final.txt)

        cmd=$(awk -F '[_\t]' '{rev="rev|sed 's/T/U/g'";sed="sed 's/T/U/g'"; printf $1"_"$7"_"$9"-"$10" unseed "$1" "$3" "$4" "$5" NA NA :target "$2" NA :tsRNA "$7"_"$9"-"$10" "$7"_"$9"-"$10" "$8" NA :evalue "$18":mismatch "$21":gapopen "$22;\
            $8=substr($8,$15,$14-$15+1);printf $8|&rev;close(rev,"to");rev|& getline out;close(rev);printf substr($2,$12,$13-$12+1)|& sed;close(sed,"to"); sed|& getline out1;close(sed);\
            printf " :tsRNA target_match "$12"\t"$13":\t\t"out1":tsRNA_match   "$14"\t\t"$15":"out"\n" }' $output/tsRNA/${finame}.6.tsRNA.target_map.blast | awk -F ':' '{print $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":\t\t"$10":"$9}' | \
            sort |	awk 'BEGIN{while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$14);$14= map[$14]; if($14~/^$/) {$14="NA";print}else{print} }' >$output/tsRNA/${finame}.6.tsRNA.unseed_final.txt)
            
        cmd=$(awk -F '[_\t]' '{rev="rev|sed 's/T/U/g'";sed="sed 's/T/U/g'"; printf $1"_"$7"_"$9"-"$10" unseed "$1" "$3" "$4" "$5" NA NA :target "$2" NA :miRNA "$7"_"$9"-"$10" "$7"_"$9"-"$10" "$8" NA :evalue "$18":mismatch "$21":gapopen "$22;\
            $8=substr($8,$15,$14-$15+1);printf $8|&rev;close(rev,"to");rev|& getline out;close(rev);printf substr($2,$12,$13-$12+1)|& sed;close(sed,"to"); sed|& getline out1;close(sed);\
            printf " :miRNA target_match "$12"\t"$13":\t\t"out1":miRNA_match   "$14"\t\t"$15":"out"\n" }' $output/tsRNA/${finame}.6.miRNA.target_map.blast | awk -F ':' '{print $1":"$2":"$3":"$4":"$5":"$6":"$7":"$8":\t\t"$10":"$9}' | \
            sort |	awk 'BEGIN{while( (getline line < "'$name'")>0) {split(line,f);map[f[1]]=f[2];}}{gsub(".*#","",$14);$14= map[$14]; if($14~/^$/) {$14="NA";print}else{print} }' >$output/tsRNA/${finame}.6.miRNA.unseed_final.txt)

        cmd=$(cat $output/tsRNA/${finame}.6.tsRNA.seed_final.txt > $output/tsRNA_${finame}.match_final.txt)
        cmd=$(cat $output/tsRNA/${finame}.6.tsRNA.unseed_final.txt >> $output/tsRNA_${finame}.match_final.txt)

        cmd=$(cat $output/tsRNA/${finame}.6.miRNA.seed_final.txt > $output/miRNA_${finame}.match_final.txt)
        cmd=$(cat $output/tsRNA/${finame}.6.miRNA.unseed_final.txt >> $output/miRNA_${finame}.match_final.txt)

        if [  -f $output/miRNA_${finame}.match_final.txt ] && [  -f $output/tsRNA_${finame}.match_final.txt ];
        then

            if [[ ! -d $output/tsRNA/cenet/ ]]; then
                cmd=$(mkdir $output/tsRNA/cenet/)
            fi

            cmd=$(awk '{print $3,$13}' $output/miRNA_${finame}.match_final.txt |sort | uniq  > $output/tsRNA/cenet/${finame}.1.miRNA_target.txt)
            cmd=$(awk '{print $3,$13}' $output/tsRNA_${finame}.match_final.txt|sort | uniq > $output/tsRNA/cenet/${finame}.1.tsRNA_target.txt)
            cmd=$(join -t $' t' -11 -21 $output/tsRNA/cenet/${finame}.1.miRNA_target.txt $output/tsRNA/cenet/${finame}.1.tsRNA_target.txt > $output/tsRNA/cenet/${finame}.2.cenet.txt)

            if [[ -s $output/tsRNA/cenet/${finame}.2.cenet.txt ]];
            then 

                cmd=$(sort -k 1 $output/tsRNA_${finame}.match_final.txt >$output/tsRNA/cenet/${finame}.2.tsRNA_match_final.txt)
                cmd=$(sort -k 1 $output/miRNA_${finame}.match_final.txt >$output/tsRNA/cenet/${finame}.2.miRNA_match_final.txt)

                cmd=$(awk '{$2=$1"_"$2;$3=$1"_"$3;print $2,$3}' $output/tsRNA/cenet/${finame}.2.cenet.txt >$output/tsRNA/cenet/${finame}.3.tr_cenet.txt)

                cmd=$(join   -t $' t' -11 -21 $output/tsRNA/cenet/${finame}.3.tr_cenet.txt  $output/tsRNA/cenet/${finame}.2.miRNA_match_final.txt | sed  's/[^ ]* /NA /' |  sed 's/$/&|NA /g' |sort -k 2 | uniq   >$output/tsRNA/cenet/${finame}.4.miRNA_cenet.txt)
                cmd=$(join  -t $' t' -12 -21 $output/tsRNA/cenet/${finame}.4.miRNA_cenet.txt  $output/tsRNA/cenet/${finame}.2.tsRNA_match_final.txt  |  sed -e 's/[^ ]* //' |sort | uniq  >$output/tsRNA/cenet/${finame}.5.cenet.txt)
                cmd=$(cp $output/tsRNA/cenet/${finame}.5.cenet.txt $output/cenetwork_${finame}.match_final.txt)

            else

                echo "No common targets of tsRNA and miRNA was found"
                cmd=$(touch $output/cenetwork_${finame}.match_final.txt)

            fi
        else

            echo "No common targets of tsRNA and miRNA was found"
            cmd=$(touch $output/cenetwork_${finame}.match_final.txt)

        fi
	else
		echo "please select the correct CLIP type";
	fi
}

type $input

if [[  -d $output/tsRNA/ ]]; then
	rm -r $output/tsRNA/
fi

if [[  -d $output/miRNA/ ]]; then
	rm -r $output/miRNA/
fi

rm $output/collapsed.fa
rm $output/trimed.fa