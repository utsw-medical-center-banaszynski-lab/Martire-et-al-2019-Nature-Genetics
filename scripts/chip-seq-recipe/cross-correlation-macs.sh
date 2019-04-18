#!/bin/bash
# cross-correlation-macs.sh

script_name="cross-correlation-macs.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --File path to alignment"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.bam' [-o '/path/to/output/dir/']"
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
    module load boost/gcc/1.54.0
    module load UCSC_userApps/v317
    module load macs/2.1.0-20151222
    module load bedtools/2.25.0

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:o:vh opt
        do
            case $opt in
                f) aln_file=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln_file ]]; then
        usage
    fi


    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        one_parent=$(dirname "${aln_file}")
        out_dir=$(dirname "${one_parent}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Define the output file name, based on the first file
    aln_fn=$(basename "${aln_file}")
    output_fp=${aln_fn%.bam}

    if [ ! -d $out_dir/$output_fp.tagAlign.gz ]; then

        # Create tagAlign file
        bamToBed -i $aln_file | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | tee $out_dir/$output_fp.tagAlign | gzip -c > $out_dir/$output_fp.tagAlign.gz

    	# Subsample tagAlign file

    	nreads=10000000
        nreads_per=`echo $(( nreads/1000000 )) `
    	subsampled_tagalign=$output_fp.filt.nodup.sample.$nreads_per.se.tagAlign.gz
    	grep -v "chrM" $out_dir/$output_fp.tagAlign | shuf -n $nreads | gzip -c > $out_dir/$subsampled_tagalign

    	# Calculate Cross-correlation QC scores
    	cc_score=$subsampled_tagalign.cc.qc


        macs2 predictd -i $out_dir/$subsampled_tagalign 2>&1 | awk -F: '{print $4}' | awk '{$1=$1}1'  >> $out_dir/$cc_score

        rm $out_dir/$subsampled_tagalign
        rm $out_dir/$output_fp.tagAlign

        # Get input and output files and then print out metadata.json file
        input_files=("$aln1")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($cross_correlation_dir\/$output_fp*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$output_fp.metadata.json
        echo "* Finished."
    else
        echo "* $aln_fn cross-correlation has been run. "
    fi
}

main "$@"
