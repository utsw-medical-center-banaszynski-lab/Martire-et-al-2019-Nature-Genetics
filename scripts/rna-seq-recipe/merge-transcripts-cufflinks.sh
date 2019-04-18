#!/bin/bash
# merge-transcripts-cufflinks.sh

script_name="merge-transcripts-cufflinks.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-o  --Path to the output directory"
  echo "-a  --File path to reference annotation"
  echo "-r  --File path to reference genome"
  echo "-g  --File path to transcripts"
  echo "-v  --Version of script"
  echo "Example: $script_name -o 'experiment_152_merge' -r '/path/to/reference/hg19.fa' -a '/path/to/gtf/gencode.v19.annotation.gtf' -g 'foo_rep1.gtf foo_rep2.gtf man_rep1.gtf man_rep2.gtf'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load cufflinks/2.2.1
    module load python/2.7.x-anaconda

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :o:a:r:g:vh opt
        do
            case $opt in
                o) out=$OPTARG;;
                a) annotation=$OPTARG;;
                r) reference_genome=$OPTARG;;
                g) transcripts=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $out ]] || [[ -z $annotation ]] || [[ -z $reference_genome ]] || [[ -z $transcripts ]]; then
        usage
    fi

    # Make Out directory if doesn't exist
    if [ ! -e $out ]; then
        mkdir $out
    fi

    # Compare expression if not already done so
    out_dir=$out\/$script_name-$script_ver
    if [ ! -d $out_dir ]; then
        mkdir $out_dir

        # Make a GTF_list file
        array_transcripts=($transcripts)
        printf '%s\n' "${array_transcripts[@]}" > $out_dir\/GTF_list.txt
        echo "* Merge trasncripts... "
        cuffmerge -p 10 -o $out_dir -g $annotation -s $reference_genome $out_dir\/GTF_list.txt

        # Remove GTF list
        rm $out_dir\/GTF_list.txt

        # Get input and output files and then print out metadata.json file
        input_files=("$annotation" "$reference_genome" "${array_transcripts[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/merged.gtf)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json
        echo "* Finished."
    else
        echo "* Merged transcripts have been made."
    fi
}

main "$@"
