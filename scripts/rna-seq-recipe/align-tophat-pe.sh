#!/bin/bash
# align-tophat-pe.sh

script_name="align-tophat-pe.sh"
script_ver="1.0.0"


#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --Path to read1 fastq files a comma-separated list"
  echo "-s  --Path to read1 fastq files a comma-separated list"
  echo "-i  --Path to BOWTIE2 index of genome"
  echo "-a  --Path to the annotation file"
  echo "-t  --Path to the index of annotation file"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.fastq.gz' -s 'foo_2.fastq.gz' -i '/path/to/Bowtie2Index/genome' -a '/path/to/gtf/gencode.v19.annotation.gtf' -t '/path/to/gtf/index/gencode.v19.annotation' [-o '/path/to/output/dir/']"
  echo "Example: $script_name -f 'foo_1.fastq.gz,fastq2_1.fastq.gz' -s 'foo_2.fastq.gz,fastq2_2.fastq.gz' -i '/path/to/Bowtie2Index/genome' -a '/path/to/gtf/gencode.v19.annotation.gtf' -t '/path/to/gtf/index/gencode.v19.annotation' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load perl/5.18.2
    module load bowtie2/gcc/2.2.3
    module load tophat/2.0.12
    module load iGenomes/2013-03-25
    module load python/2.7.x-anaconda
    module load samtools/intel/0.1.19

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:i:a:t:o:vh opt
        do
            case $opt in
                f) raw_file_1=$OPTARG;;
                s) raw_file_2=$OPTARG;;
                i) annotation=$OPTARG;;
                a) transcriptome=$OPTARG;;
                t) transcriptome_index=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $raw_file_1 ]] || [[ -z $raw_file_2 ]] || [[ -z $annotation ]] || [[ -z $transcriptome ]] || [[ -z $transcriptome_index ]]; then
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    raw_file_array_1=(${raw_file_1//,/ })
    raw_file_array_2=(${raw_file_2//,/ })
    if [ -z $out ]; then
        out_dir=$(dirname "${raw_file_array[0]}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [[(${#raw_file_array_1[@]} -ne ${#raw_file_array_2[@]}) ]]; then
        echo "The number of arguments are not equal."
        exit 1
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    echo "* Value of reads: '$raw_file_1'"
    echo "* Value of reads: '$raw_file_2'"
    echo "* Value of genome index : '$annotation'"
    echo "* Value of annotation: '$transcriptome'"

    # Align if file doesn't exist
    aln_fn=$(basename "${raw_file_array_1[0]}")
    if [ ! -e $out_dir/$aln_fn.accepted_hits.bam ]; then
        echo "* Map reads..."
        tophat -p 8 --keep-fasta-order --no-coverage-search --library-type fr-firststrand -G $transcriptome \
            --transcriptome-index $transcriptome_index -g 10 --output-dir $out_dir $annotation $raw_file_1 $raw_file_2

        # Rename and move files
        echo "* Moving files..."
        rename_files=("accepted_hits.bam" "align_summary.txt"  "deletions.bed"  "insertions.bed"  "junctions.bed"  "prep_reads.info" "unmapped.bam")
        for file in ${rename_files[@]}
        do
            mv $out_dir\/$file $out_dir\/$aln_fn.$file
        done

        # Get input and output files and then print out metadata.json file
        input_files=("$transcriptome" "$transcriptome_index" "$annotation" "${raw_file_array_1[@]}", "${raw_file_array_2[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$aln_fn*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$aln_fn.metadata.json
        echo "* Finished."
    else
        echo "* $aln_fn has already been aligned"
    fi
}

main "$@"
