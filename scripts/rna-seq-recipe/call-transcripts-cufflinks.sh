#!/bin/bash
# call-transcripts-cufflinks.sh

script_name="call-transcripts-cufflinks.sh"
script_ver="1.0.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-f  --File path to alignment"
  echo "-g  --Path to genome"
  echo "-a  --Path to the annotation file"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo_1.bam' -g '/path/to/genome.fa' -a '/path/to/gtf/gencode.v19.annotation.gtf'"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}


main(){

    # Load required modules
    module load bowtie2/intel/2.1.0
    module load samtools/0.1.19
    module load tophat/2.0.12
    module load cufflinks/2.2.1
    module load iGenomes/2013-03-25
    module load python/2.7.x-anaconda

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:g:a:vh opt
        do
            case $opt in
                f) alignment_file=$OPTARG;;
                g) reference_genome=$OPTARG;;
                a) transcriptome=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $alignment_file ]] || [[ -z $reference_genome ]] || [[ -z $transcriptome ]]; then
        usage
    fi

    # Compute transcripts if not already done so
    alignment_fn=$(basename "${alignment_file}")
    first_parent=$(dirname "${alignment_file}") # Get the first parent directory and go up
    out_dir=$(dirname "${first_parent}")\/$script_name-$script_ver/
    cwd=`pwd`
    if [ ! -e $out_dir/$alignment_fn.transcripts.gtf ]; then
        mkdir $out_dir
        echo "* Value of alignments: '$alignment_fn'"
        echo "* Value of annotation: '$transcriptome'"
        echo "* Value of genome: '$reference_genome'"
        cd $out_dir # patch fix workaround for cufflinks
        cufflinks -p 10 --library-type fr-firststrand -G $transcriptome -b $reference_genome $alignment_file
        cd $cwd
        # Rename and move files
        rename_files=("transcripts.gtf" "skipped.gtf" "isoforms.fpkm_tracking"  "genes.fpkm_tracking")
        for file in ${rename_files[@]}
        do
            mv $out_dir\/$file $out_dir\/$alignment_fn.$file
        done

        # Get input and output files and then print out metadata.json file
        input_files=("$transcriptome" "$reference_genome" "$alignment_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$alignment_fn*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$alignment_fn.metadata.json
        echo "* Finished."
    else
        echo "* Transcripts have been called for $alignment_fn"
    fi
}

main "$@"
