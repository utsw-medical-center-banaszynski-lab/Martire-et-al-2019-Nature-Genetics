#!/bin/bash
# compare-expression-cufflinks.sh

script_name="compare-expression-cufflinks.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-o  --Path to the output directory"
  echo "-e  --The experiment names, comma seperated"
  echo "-g  --Path to the combined transcripts file"
  echo "-t  --Flag if series of experiments is a time series"
  echo "-a  --File path to alignments/cxb, comma seperated for replicates and space sperated for experiemnts"
  echo "-r  --File path to reference genome"
  echo "-v  --Version of script"
  echo "Example: $script_name -o 'experiment_152_comparison' -e 'foo,man' -a 'foo_rep1.bam,foo_rep2.bam man_rep1.bam,man_rep2.bam' -r '/path/to/genome.fa' -g '/path/to/merged.gtf' [ -t]"
  echo "Example: $script_name -o 'experiment_152_comparison' -e 'foo,man' -a 'foo_rep1.cxb,foo_rep2.cxb man_rep1.cxb,man_rep2.cxb' -r '/path/to/genome.fa' -g '/path/to/merged.gtf' [ -t]"
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
    while getopts :o:e:g:ta:r:vh opt
        do
            case $opt in
                o) out=$OPTARG;;
                e) experiments=$OPTARG;;
                g) transcripts=$OPTARG;;
                t) time_series="time series";;
                a) alignments=$OPTARG;;
                r) reference_genome=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $out ]] || [[ -z $experiments ]] || [[ -z $transcripts ]] || [[ -z $alignments ]] || [[ -z $reference_genome ]]; then
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

        echo "* Value of experiments: '$experiments'"
        # Check if time series is
        if [[ $time_series = 'time series' ]]; then
            echo "* Computing Differential Expression as Time Series..."
            cuffdiff -L $experiments -o $out_dir -p 10 --library-type fr-firststrand -b $reference_genome $transcripts $alignments -T
        else
            echo "* Computing Differential Expression..."
            cuffdiff -L $experiments -o $out_dir -p 10 --library-type fr-firststrand -b $reference_genome $transcripts $alignments
        fi
        # Get input and output files and then print out metadata.json file
        alignments_array=(${alignments//[,| ]/ })
        input_files=("$transcripts" "$reference_genome" "${alignments_array[@]}")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        option_list=("$time_series")
        printf -v option "\"%s\"," "${option_list[@]}"
        option=${option%,}

        # If time series is used output in metadata
        if [[ $time_series = 'time series' ]]; then
            printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s], "options": [%s]}' "$script_name" "$script_ver" "$input" "$output" "$option" | python -m json.tool > $out_dir/metadata.json
        else
            printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json
        fi
        echo "* Finished."
    else
        echo "* Expression comparisons have been made."
    fi
}

main "$@"
