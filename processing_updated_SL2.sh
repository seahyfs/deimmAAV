#!/bin/bash
# Stop on errors.
set -e
set +x

PROGRAM="NGS | PHIP-seq Analysis"

### CHANGE HERE
THREADS=16
REFERENCE="small_library_trim.fasta"
### END

BASEDIR=`pwd` #$(dirname "$0")

function check_file {
  if [ ! -f $1 ]; then
    echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} $1 Does not Exist"
    exit_abnormal
  fi
}

function create_folder {
  if [ -d $1 ]; then
    rm -r $1
  fi
  
  mkdir $1
  echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} $1 Output directory: $OUTPUT/$1"
}

function check_tool {
    if [ -z `which $1` ]; then
        echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} Please ensure $1 is installed properly, Suggestion: see install.md"
        exit_abnormal
    fi
}

function SORT_SAM {
  samtools sort -O BAM
}

function DEPTH_SAM {
  samtools depth -aa -m 100000000 - 
}

function AWK_PROCESS {
  awk -F "\t" 'BEGIN {OFS = FS} {counts[$1] = ($3 < counts[$1]) ? counts[$1] : $3} END {for (c in counts) {print c, counts[c]}}'
}

function bowtie_process {
  #remove .fastq
  #FASTQ_1=(${$1//./ })
  #FASTQ_2=(${$2//./ })

  echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} PROCESSING $1 $2"
  IFS='.' read -ra FASTQ_1 <<< "$1"
  IFS='.' read -ra FASTQ_2 <<< "$2"
  
  gzip -d $1
  fastx_trimmer -f 54 -i ${FASTQ_1[0]}.fastq -o ${FASTQ_1[0]}_trim.fastq
  
  gzip -d $2
  fastx_trimmer -f 33 -i ${FASTQ_2[0]}.fastq -o ${FASTQ_2[0]}_trim.fastq

  IFS='_' read -ra ADDR <<< $1
  TITLE="${ADDR[0]}"
  FILE="${BASEDIR}/sample_counts/${TITLE}.tsv"
  echo -e "id\t${TITLE}" > $FILE
  bowtie2 -p 16 -x oligosref2 -1 ${FASTQ_1[0]}_trim.fastq -2 ${FASTQ_2[0]}_trim.fastq --threads $THREADS | SORT_SAM | DEPTH_SAM | AWK_PROCESS | sort -k 1 >> $FILE
      #| samtools sort -O BAM \\
      #| samtools depth -aa -m 100000000 - \\
      #| awk 'BEGIN {OFS=\"\\t\"} {counts[$1] = ($3 < counts[$1]) ? counts[$1] : $3} END {for (c in counts) {print c, counts[c]}}' \\
      #| sort -k 1 \\
      #>> FILE
  gzip ${FASTQ_1[0]}.fastq 
  gzip ${FASTQ_2[0]}.fastq
}

function phip_seq_process {
  phip merge-columns -m outer -i sample_counts -o counts.tsv
}

check_tool phip
check_tool bowtie2
check_tool gzip
check_tool fastx_trimmer

# use nullglob in case there are no matching files
shopt -s nullglob

echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} GETTING FILES"
# create an array with all the filer/dir inside ~/myDir
arr=($( ls * ))
declare -a FASTQ

for ((i=0; i<${#arr[@]}; i++)); do
  if [[ ${arr[$i]} == *".fastq.gz" ]] || [[ ${arr[$i]} == *".fq.gz" ]]  ; then
    FASTQ+=(${arr[$i]})
  fi
done

for ((i=0; i<${#FASTQ[@]}; i++)); do
  #do something to each element of array
  echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} 1 ${FASTQ[$i]}"
done

IFS=$'\n' FASTQ_SORTED=($(sort <<<"${FASTQ[*]}")); unset IFS

FASTQ_LENGTH=${#FASTQ[@]}
if (( FASTQ_LENGTH % 2 != 0 )); then
  for (( i=0; i<${#FASTQ_SORTED[@]} ; i+=2 )) ; do
    #do something to each element of array
    echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} ${FASTQ_SORTED[i]} ${FASTQ_SORTED[i+1]}"
  done
  
  echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} Error! Number of Fastq Files are not even!"
  exit 125
fi


create_folder sample_counts
bowtie2-build $REFERENCE oligosref2
echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} BOWTIE PROCESSING"
for (( i=0; i<${#FASTQ_SORTED[@]} ; i+=2 )) ; do
  bowtie_process ${FASTQ_SORTED[i]} ${FASTQ_SORTED[i+1]}
done
create_folder sample_counts_inputs
mv sample_counts/Input-1.tsv sample_counts_inputs/.
mv sample_counts/Input-2.tsv sample_counts_inputs/.

echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} BOWTIE PROCESSING COMPLETED"
phip_seq_process
echo $(date +"%Y-%m-%d %H:%M:%S") "- ${PROGRAM} Completed!"