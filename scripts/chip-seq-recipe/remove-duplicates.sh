#!/bin/bash
# remove-duplicates.sh

script_name="remove-duplicates.sh"
script_ver="1.1.0"

#Help function
usage() {
  echo "-h Help documentation for $script_name"
  echo "-a  --Path to bam file"
  echo "-s  --Type of reads (e.g paired-end, single-end)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -a 'foo_1.bam' -s 'single-end' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}

main(){


    # Load required modules
    module load samtools/0.1.19
    module load bedtools/2.17.0
    module load python/2.7.x-anaconda
    module load picard/1.127

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :a:s:o:vh opt
        do
            case $opt in
                a) aln_file=$OPTARG;;
                s) sequencing_type=$OPTARG;;
                o) out=$OPTARG;;
                v) version;;
                h) usage;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln_file ]] || [[ -z $sequencing_type ]]; then
        usage
    fi


    if [ $sequencing_type != 'single-end' ] && [ $sequencing_type != 'paired-end' ]; then
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

    # Define the output file name, based on the file
    raw_fn=$(basename "${aln_file}")
    output_fp=${raw_fn%_*}

    # If filtered alignment file doesn't exist
    if [ ! -e $out_dir/$output_fp.filtered.no_dups.bam ]; then

		# Remove unmapped, mate unmapped
		# not primary alignment, reads failing platform
		# Remove low MAPQ reads
		# Obtain name sorted BAM file
        if [ $sequencing_type = 'single-end' ]; then
            samtools view -F 1804 -q 30 -b $aln_file > $out_dir/$output_fp.filtered.bam
        else
            samtools view -F 1804 -f 2 -q 30 -u $aln_file | samtools sort -n - -o $out_dir/tmp.$output_fp
            #fill in mate coordinates, ISIZE and mate-related flags
			#fixmate requires name-sorted alignment; -r removes secondary and unmapped (redundant here because already done above?)
			#- send output to stdout
            samtools fixmate -r $out_dir/tmp.$output_fp - | samtools view -F 1804 -f 2 -u - | samtools sort - -o $out_dir/$output_fp.filtered.bam
            rm $out_dir/tmp.$output_fp
        fi

        # Mark duplicates
        java -Xmx4G -jar $PICARD/picard.jar MarkDuplicates INPUT=$out_dir/$output_fp.filtered.bam OUTPUT=$out_dir/$output_fp.tmp.bam METRICS_FILE=$out_dir/$output_fp.duplicate_metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
        mv $out_dir/$output_fp.tmp.bam $out_dir/$output_fp.filtered.bam

        # Remove duplicates and Index final position sorted BAM
        if [ $sequencing_type = 'single-end' ]; then
            samtools view -F 1804 -b $out_dir/$output_fp.filtered.bam > $out_dir/$output_fp.filtered.no_dups.bam
            samtools index $out_dir/$output_fp.filtered.no_dups.bam $out_dir/$output_fp.filtered.no_dups.bai
        else
            samtools view -F 1804 -b $out_dir/$output_fp.filtered.bam > $out_dir/$output_fp.filtered.no_dups.bam
        fi

        # Generate mapping statistics
        samtools flagstat $out_dir/$output_fp.filtered.no_dups.bam > $out_dir/$output_fp.filtered.no_dups_stats.txt

        # Compute library complexity
        # Sort by name
        # convert to bedPE and obtain fragment coordinates
        # sort by position and strand
        # Obtain unique count statistics
        # PBC File output
    	# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
        if [ $sequencing_type = 'single-end' ]; then
            bamToBed -i $out_dir/$output_fp.filtered.bam  > $out_dir/$output_fp.filtered.bed
            printf "TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair\n" > $out_dir/$output_fp.pbc_metrics.txt
            awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' $out_dir/$output_fp.filtered.bed \
            | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> $out_dir/$output_fp.pbc_metrics.txt
        else
            samtools sort -n $out_dir/$output_fp.filtered.bam -o - | bamToBed -bedpe -i stdin > $out_dir/$output_fp.filtered.bed
            printf "TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF=Distinct/Total\tPBC1=OnePair/Distinct\tPBC2=OnePair/TwoPair\n" > $out_dir/$output_fp.pbc_metrics.txt
            awk 'BEGIN{OFS="\t"}{print $1,$2,$4,$6,$9,$10}' $out_dir/$output_fp.filtered.bed \
            | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' >> $out_dir/$output_fp.pbc_metrics.txt
        fi

        rm $out_dir/$output_fp.filtered.bed
        rm $out_dir/$output_fp.filtered.bam

        # Get input and output files and then print out metadata.json file
        input_files=("$aln_file")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir/$output_fp*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/$output_fp.metadata.json

    else
        echo "* $raw_fn has already removed duplicates"
    fi
}

main "$@"
