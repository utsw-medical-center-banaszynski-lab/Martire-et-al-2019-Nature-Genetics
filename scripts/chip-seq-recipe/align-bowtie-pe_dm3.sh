#!/bin/bash
# align-bowtie-pe.sh

script_name="align-bowtie-pe.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --Path to read1 fastq files a comma-separated list"
  echo "-s  --Path to read2 fastq files a comma-separated list"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.fastq.gz' -s 'foo_2.fastq.gz' -r 'hg19' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load bowtie2/2.2.8-intel
    module load samtools/0.1.19
    module load bedtools/2.17.0
    module load iGenomes/2013-03-25

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:r:o:vh opt
        do
            case $opt in
                f) raw_file_1=$OPTARG;;
                s) raw_file_2=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $raw_file_1 ]] || [[ -z $raw_file_2 ]] || [[ -z $ucsc_reference ]]; then
        usage
    fi

    # Define the index to use
    # Establish the genome file version to use and BOWTIEINDEX
    if [ $ucsc_reference = 'hg19' ]; then
        index=$iGENOMES_DB_DIR\/Homo_sapiens/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
    elif [ $ucsc_reference = 'mm10' ]; then
        index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
    elif [ $ucsc_reference = 'mm9' ]; then
        index=$iGENOMES_DB_DIR\/Mus_musculus/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
    elif [ $ucsc_reference = 'dm3' ]; then
        index=$iGENOMES_DB_DIR\/Drosophila_melanogaster/UCSC/$ucsc_reference/Sequence/Bowtie2Index/genome
    else
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        out_dir=$(dirname "${raw_file_1}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Define the output file name, based on the first file
    raw_fn_1=$(basename "${raw_file_1}")
    raw_fn_2=$(basename "${raw_file_2}")
    output_fp=${raw_fn_1%_*}

    # Align if file doesn't exist
    if [ ! -e $out_dir/$output_fp.sorted.bam ]; then

        #unizip files
        file_un_zip_1=${raw_fn_1%.gz}
        file_un_zip_2=${raw_fn_2%.gz}
        gunzip -c $raw_file_1 > $out_dir/$file_un_zip_1
        gunzip -c $raw_file_2 > $out_dir/$file_un_zip_2

        bowtie2 -t -x $index -1 $out_dir/$file_un_zip_1 -2 $out_dir/$file_un_zip_2 -S $out_dir/$output_fp.sam

        # Convert sam to bam
        samtools view -bh -S $out_dir/$output_fp.sam > $out_dir/$output_fp.unsorted.bam

        # Sort bam
    	samtools sort $out_dir/$output_fp.unsorted.bam $out_dir/$output_fp.sorted

    	rm $out_dir/$output_fp.unsorted.bam
        rm $out_dir/$output_fp.sam
        rm $out_dir/$file_un_zip_1 $out_dir/$file_un_zip_2

        # Get input and output files and then print out metadata.json file
        input_files=("$index" "$raw_file_1" "$raw_file_2")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$output_fp*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$output_fp.metadata.json

    else
        echo "* $output_fp has already been aligned"
    fi
}

main "$@"
