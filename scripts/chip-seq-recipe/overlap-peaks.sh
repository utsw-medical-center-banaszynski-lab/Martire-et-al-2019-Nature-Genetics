#!/bin/bash
# overlap-peaks.sh

script_name="overlap-peaks.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  -- File path to Rep1 and 2 peaks."
  echo "-s  -- File path to Pooled Peak"
  echo "-r  -- File path to Pseudo Peak files"
  echo "-p  -- Peak Type (narrowPeak, gappedPeak, broadPeak )"
  echo "-u  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.narrowPeak foo2.narrowPeak' -s 'foo.narrowPeak' -r 'pseduofoo1.narrowPeak pseduofoo2.narrowPeak' -p 'narrowPeak' -u 'hg19' [-o '/path/to/output/dir/']"
  exit 1
}


main(){



    # Load required modules
    module load python/2.7.x-anaconda
    module load gcc/4.8.1
    module load bedtools/2.17.0
    module load UCSC_userApps/v317

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:r:p:u:o:vh opt
        do
            case $opt in
                f) exp_peak=$OPTARG;;
                s) pool_peak=$OPTARG;;
                r) pseudo_peak=$OPTARG;;
                p) peak_type=$OPTARG;;
                u) ucsc_reference=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $exp_peak ]] || [[ -z $pool_peak ]] || [[ -z $pseudo_peak ]] || [[ -z $peak_type ]]; then
        usage
    fi

    # Check if length of arguments in aln1, aln2 and exp are the same
    exp_peak_array=(${exp_peak//[,| ]/ })
    pseudo_peak_array=(${pseudo_peak//[,| ]/ })
    if [[(${#exp_peak_array[@]} -ne 2) && (${#pseudo_peak_array[@]} -ne 2)]]; then
        echo "The number of arguments are not enough."
        exit 1
    fi

    # Define the genome size to use
    if [[ $peak_type != 'narrowPeak' ]] && [[ $peak_type != 'broadPeak' ]] && [[ $peak_type != 'gappedPeak' ]]; then
        echo "Needs to be narrowPeak, broadPeak, or gappedPeak"
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        one_parent=$(dirname "${exp_peak_array[0]}")
        out_dir=$(dirname "${one_parent}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi


    # Align if file doesn't exist
    if [ ! -e $out_dir/metadata.json ]; then

        echo "* Downloading chrom.sizes..."
        fetchChromSizes $ucsc_reference > $out_dir/chrom.sizes
        chrom_sizes=$out_dir/chrom.sizes

        # Run overlap peaks
        python overlap-peaks.py ${exp_peak_array[0]} ${exp_peak_array[1]} $pool_peak ${pseudo_peak_array[0]} ${pseudo_peak_array[1]} $chrom_sizes $out_dir $peak_type

        # Remove chrom sizes
        rm $out_dir/chrom.sizes

        # Get input and output files and then print out metadata.json file
        input_files=("${exp_peak_array[@]}" "${pseudo_peak_array[@]}" "$pool_peak", "$peak_type", "$ucsc_reference" )
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json

    else
        echo "* Overlap peaks have already been called."
    fi
}


main "$@"
