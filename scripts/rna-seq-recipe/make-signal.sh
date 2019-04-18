#!/bin/bash
# make-signal.sh

script_name="make-signal.sh"
script_ver="1.1.1"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --File path to alignment"
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-s  --Type of reads (e.g paired-end, single-end)"
  echo "-v  --Version of script"
  echo "Example: $script_name -f '/path/to/alignment.bam' -r 'hg19' -s 'single-end'"
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
    module load bedtools/2.25.0
    module load samtools/intel/1.3
    module load star/2.4.2a
    module load UCSC_userApps/v317

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:r:s:vh opt
        do
            case $opt in
                f) alignment=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                s) sequencing_type=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $alignment ]] || [[ -z $ucsc_reference ]] || [[ -z $sequencing_type ]]; then
        usage
    fi

    if [ $sequencing_type != 'single-end' ] && [ $sequencing_type != 'paired-end' ]; then
        usage
    fi

    # Compute transcripts if not already done so
    alignment_fn=$(basename "${alignment}")
    first_parent=$(dirname "${alignment}") # Get the first parent directory and go up
    out_dir=$(dirname "${first_parent}")\/$script_name-$script_ver

    alignment_prefix=${alignment_fn%.fastq.gz.accepted_hits.bam}
    # Make Out directory if doesn't exist
    if [ ! -e $out_dir ]; then
        mkdir $out_dir
    fi

    if [ ! -e $out_dir/$aligment_prefix.metadata.json ]; then
        echo "* Value of annotations: '$alignment_fn'"
        echo "* Value of genome: '$ucsc_reference'"

        # Fetch chrom sizes
        echo "* Downloading chrom.sizes..."
        fetchChromSizes $ucsc_reference > $out_dir/chrom.sizes

        # Make signal for Positive and Negative strand

        if [ $sequencing_type = 'single-end' ]; then
            echo "* Generating Positive Signal... "
            bedtools genomecov -ibam $alignment -g $out_dir/chrom.sizes -split -bg -strand + > $out_dir/$alignment_prefix.positive.bedGraph
            bedSort $out_dir/$alignment_prefix.positive.bedGraph $out_dir/$alignment_prefix.positive.sorted.bedgraph
            bedGraphToBigWig $out_dir/$alignment_prefix.positive.sorted.bedgraph $out_dir/chrom.sizes $out_dir/$alignment_prefix.positive.bw
            rm $out_dir/$alignment_prefix.positive.bedgraph
            rm $out_dir/$alignment_prefix.positive.sorted.bedgraph

            echo "* Generating Negative Signal... "
            bedtools genomecov -ibam $alignment -g $out_dir/chrom.sizes -split -bg -strand - > $out_dir/$alignment_prefix.negative.bedGraph
            bedSort $out_dir/$alignment_prefix.negative.bedGraph $out_dir/$alignment_prefix.negative.sorted.bedgraph
            bedGraphToBigWig $out_dir/$alignment_prefix.negative.sorted.bedgraph $out_dir/chrom.sizes $out_dir/$alignment_prefix.negative.bw
            rm $out_dir/$alignment_prefix.negative.bedgraph
            rm $out_dir/$alignment_prefix.negative.sorted.bedgraph
        elif [ $sequencing_type = 'paired-end' ]; then
            cd $out_dir
            STAR --runMode inputAlignmentsFromBAM --inputBAMfile $alignment --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $alignment_prefix --outWigReferencesPrefix chr
            bedSort $out_dir/$alignment_prefix\Signal.UniqueMultiple.str1.out.bg $out_dir/$alignment_prefix\Signal.UniqueMultiple.str1.out.sorted.bg
            bedSort $out_dir/$alignment_prefix\Signal.UniqueMultiple.str2.out.bg $out_dir/$alignment_prefix\Signal.UniqueMultiple.str2.out.sorted.bg
            bedGraphToBigWig $out_dir/$alignment_prefix\Signal.UniqueMultiple.str1.out.sorted.bg $out_dir/chrom.sizes $out_dir/$alignment_prefix.negative.bw
            bedGraphToBigWig $out_dir/$alignment_prefix\Signal.UniqueMultiple.str2.out.sorted.bg $out_dir/chrom.sizes $out_dir/$alignment_prefix.positive.bw
            rm $out_dir/$alignment_prefix*.bg
        fi
        rm $out_dir/chrom.sizes


        # Get input and output files and then print out metadata.json file
        input_files=("$alignment")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/$alignment_prefix*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$alignment_prefix.metadata.json
        echo "* Finished."
    else
        echo "* Signals have already been called for $alignment_fn"
    fi
}

main "$@"
