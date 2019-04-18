#!/bin/bash
# quality-metric-fastqc.sh

script_name="quality-metric-fastqc.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --Path to fastq file"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.fastq.gz' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load fastqc/0.11.2

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:o:vh opt
        do
            case $opt in
                f) raw_file=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $raw_file ]]; then
        usage
    fi

    #Run FastQC if file doesn't exist
    raw_fn=$(basename "${raw_file}")
    raw_prefix=${raw_fn%.fastq.gz}
    if [[ -z $out ]]; then
        out_dir=$(dirname "${raw_file}")\/$script_name-$script_ver
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    # Make Out directory if doesn't exist
    if [ ! -e $out_dir ]; then
        mkdir $out_dir
    fi

    # Run FastQC
    if [ ! -e $out_dir\/$raw_prefix\_fastqc.zip ]; then
        echo "* Value of annotations: '$raw_fn'"
        echo "* Generating Report... "
        fastqc $raw_file -o $out_dir

        # Unzip key summary.txt and fastqc_data.txt to be viewed
        unzip -p $out_dir\/$raw_prefix\_fastqc.zip $raw_prefix\_fastqc\/summary.txt > $out_dir\/$raw_prefix.summary.txt
        unzip -p $out_dir\/$raw_prefix\_fastqc.zip $raw_prefix\_fastqc\/fastqc_data.txt > $out_dir\/$raw_prefix.fastqc_data.txt

        # Get input and output files and then print out metadata.json file
        input_files=("$raw_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$raw_prefix*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$raw_fn.metadata.json
        echo "* Finished."
    else
        echo "* QC metrics has been calculated $raw_fn"
    fi
}

main "$@"
