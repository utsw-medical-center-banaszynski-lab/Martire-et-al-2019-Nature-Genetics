#!/bin/bash
# cross-correlation.sh

script_name="cross-correlation.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --File path to Replicate 1 alignments"
  echo "-s  --File path to Replicate 2 alignments"
  echo "-e  --The experiment name"
  echo "-a  --Path to the annotation file"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.bam' -s 'foo2.bam' -a '/path/to/gtf/annotation.gtf' -e 'foo' -o '/path/to/output/dir/' -r 'hg19'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){

    # Load required modules
    module load python/2.7.x-anaconda
    module load R/3.2.1-intel
    module load UCSC_userApps/v317

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:a:o:e:r:vh opt
        do
            case $opt in
                f) aln1=$OPTARG;;
                s) aln2=$OPTARG;;
                a) annotation=$OPTARG;;
                o) out=$OPTARG;;
                e) experiment=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln1 ]] || [[ -z $aln2 ]] || [[ -z $annotation ]] || [[ -z $out ]] || [[ -z $experiment ]] || [[ -z $ucsc_reference ]]; then
        usage
    fi

    # Define the output directory
    out_dir=$out\/$experiment

    # Make sure directories exist
    if [ ! -e $out/ ]; then
        mkdir $out
    fi

    if [ ! -e $out_dir/ ]; then
        mkdir $out_dir
    fi

    # Fetch chrom sizes
    echo "* Downloading chrom.sizes..."
    fetchChromSizes $ucsc_reference > $out_dir/chrom.sizes

    # Run cross-correlation if directory doesn't exist
    cross_correlation_dir=$out_dir\/$script_name-$script_ver
    if [ ! -d $cross_correlation_dir ]; then
        mkdir $cross_correlation_dir
        echo "* Running correlation..."
        Rscript cross-correlation.R --replicate1 $aln1 --replicate2 $aln2 --annotation $annotation --out $cross_correlation_dir

        rm $out_dir/chrom.sizes

        # Get input and output files and then print out metadata.json file
        input_files=("$aln1" "$aln2" "$annotation")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($cross_correlation_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $cross_correlation_dir/metadata.json
        echo "* Finished."
    else
        echo "* Cross-correlation has been run. "
    fi
}

main "$@"
