#!/bin/bash
# call-peaks-macs.sh

script_name="call-peaks-macs2.sh"
script_ver="2.0.0"

#Help function
usage() {
  echo "-h  --Help documentation for $script_name"
  echo "-f  --File path to experiment tagAlign files."
  echo "-s  --File path to control tagAlign files."
  echo "-x  --File path to experiment cross-correlation scores."
  echo "-r  --UCSC Reference genome (e.g. hg19, mm10)"
  echo "-o  --Path to output directory"
  echo "-v  --Version of script"
  echo "Example: $script_name -f 'foo1.tagAlign.gz foo2.tagAlign.gz' -s 'con1.tagAlign.gz con2.tagAlign.gz' -x 'foo1.cc foo2.cc' -r 'hg19' [-o '/path/to/output/dir/']"
  exit 1
}

# Version function
version(){
    echo "$script_name $script_ver"
    exit 1
}


# Peak calling function
call_peak() {
    # Establish variables

    experiment=$1
    control=$2
    xcor_scores_input=$3
    chrom_sizes=$4
    genomesize=$5
    out_dir=$6

    #Extract the fragment length estimate from cross-correlation scores file
    fraglen=`cat $xcor_scores_input |  grep "predicted" | cut -f6 -d ' '`
    echo $fraglen

    # Generate narrow peaks and preliminary signal tracks
    base_fn=$(basename "${experiment}")
    prefix=${base_fn%.tagAlign.gz}
    macs2 callpeak -t $experiment -c $control -f BED -n $out_dir/$prefix -g $genomesize -p 1e-2 --nomodel --shift 0 --extsize $fraglen --keep-dup all -B --SPMR

    # Rescale Col5 scores to range 10-1000 to conform to narrowPeak.as format (score must be <1000)
    python rescale.py $out_dir/$prefix\_peaks.narrowPeak

    # Sort by Col8 in descending order and replace long peak names in Column 4 with Peak_<peakRank>
    sort -k 8gr,8gr $out_dir/rescaled-$prefix\_peaks.narrowPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'| tee $out_dir/$prefix.narrowPeak
    sort -k1,1 -k2,2n  $out_dir/$prefix.narrowPeak >  $out_dir/$prefix.sorted.narrowPeak

    # Generate bigBeds from beds
    bedToBigBed $out_dir/$prefix.sorted.narrowPeak -type=bed6+4 -as=./resources/as/narrowPeak.as $chrom_sizes  $out_dir/$prefix.narrowPeak.bb
    rm $out_dir/sorted-$prefix\_peaks.narrowPeak
    rm $out_dir/rescaled-$prefix\_peaks.narrowPeak

    # Generate broad and gapped peaks
    macs2 callpeak -t $experiment -c $control -f BED -n $out_dir/$prefix -g $genomesize -p 1e-2 --broad --nomodel --shift 0 --extsize $fraglen --keep-dup all

    # Rescale Col5 scores to range 10-1000 to conform to broadPeak.as and gappedPeak.as format (score must be <1000)
    python rescale.py $out_dir/$prefix\_peaks.broadPeak
    python rescale.py $out_dir/$prefix\_peaks.gappedPeak

    # Sort by Col8 (for broadPeak) or Col 14(for gappedPeak)  in descending order and replace long peak names in Column 4 with Peak_<peakRank>
    sort -k 8gr,8gr $out_dir/rescaled-$prefix\_peaks.broadPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'| tee $out_dir/$prefix.broadPeak
    sort -k1,1 -k2,2n  $out_dir/$prefix.broadPeak >  $out_dir/$prefix.sorted.broadPeak
    sort -k 14gr,14gr $out_dir/rescaled-$prefix\_peaks.gappedPeak | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR ; print $0}'| tee $out_dir/$prefix.gappedPeak
    sort -k1,1 -k2,2n  $out_dir/$prefix.gappedPeak >  $out_dir/$prefix.sorted.gappedPeak

    # Generate bigBeds from beds
    bedToBigBed $out_dir/$prefix.sorted.broadPeak -type=bed6+3 -as=./resources/as/broadPeak.as $chrom_sizes  $out_dir/$prefix.broadPeak.bb
    bedToBigBed $out_dir/$prefix.sorted.gappedPeak -type=bed12+3 -as=./resources/as/gappedPeak.as $chrom_sizes  $out_dir/$prefix.gappedPeak.bb
    rm $out_dir/sorted-$prefix\_peaks.broadPeak
    rm $out_dir/rescaled-$prefix\_peaks.broadPeak
    rm $out_dir/sorted-$prefix\_peaks.gappedPeak
    rm $out_dir/rescaled-$prefix\_peaks.gappedPeak

    # Generate fold enrichment signal tracks
    macs2 bdgcmp -t $out_dir/$prefix\_treat_pileup.bdg -c $out_dir/$prefix\_control_lambda.bdg --outdir $out_dir -o $prefix\_FE.bdg -m FE

    # Genearte bigWigs from bedgraph to support vizualization
    bedtools slop -i $out_dir/$prefix\_FE.bdg -g /$chrom_sizes -b 0 | bedClip stdin $chrom_sizes $out_dir/$prefix.fc.signal.bedgraph
    bedSort $out_dir/$prefix.fc.signal.bedgraph $out_dir/$prefix.sorted.fc.signal.bedgraph
    bedGraphToBigWig $out_dir/$prefix.sorted.fc.signal.bedgraph $chrom_sizes $out_dir/$prefix.fc_signal.bw
    rm $out_dir/$prefix.fc.signal.bedgraph
    rm $out_dir/$prefix.sorted.fc.signal.bedgraph

    # Generate fold enrichment signal tracks

    # Compute sval = min(no. of reads in ChIP, no. of reads in control) / 1,000,000
    experiment_reads=`gzip -dc $experiment | wc -l | cut -f1`
    control_reads=`gzip -dc $control | wc -l | cut -f1`
    if [ $experiment_reads -ge $control_reads ]; then
        min=$control_reads
    else
        min=$experiment_reads
    fi
    sval=`echo $min/1000000 | bc -l`

    macs2 bdgcmp -t $out_dir/$prefix\_treat_pileup.bdg -c $out_dir/$prefix\_control_lambda.bdg --outdir $out_dir -o $prefix\_ppois.bdg -m ppois -S $sval

    # Genearte bigWigs from bedgraph to support vizualization
    bedtools slop -i $out_dir/$prefix\_ppois.bdg -g $chrom_sizes -b 0 | bedClip stdin $chrom_sizes $out_dir/$prefix.pval.signal.bedgraph
    bedSort $out_dir/$prefix.pval.signal.bedgraph $out_dir/$prefix.sorted.pval.signal.bedgraph
    bedGraphToBigWig $out_dir/$prefix.sorted.pval.signal.bedgraph $chrom_sizes $out_dir/$prefix.pval_signal.bw
    rm $out_dir/$prefix.pval.signal.bedgraph
    rm $out_dir/$prefix.sorted.pval.signal.bedgraph
}

main(){

    # Load required modules
    module load python/2.7.x-anaconda
    module load R/3.1.0-intel
    module load macs
    module load gcc/4.8.1
    module load bedtools/2.17.0
    module load UCSC_userApps/v317

    # Parsing options
    OPTIND=1 # Reset OPTIND
    while getopts :f:s:x:r:o:hv opt
        do
            case $opt in
                f) aln=$OPTARG;;
                s) control=$OPTARG;;
                x) aln_cross=$OPTARG;;
                r) ucsc_reference=$OPTARG;;
                o) out=$OPTARG;;
                h) usage;;
                v) version;;
            esac
        done

    shift $(($OPTIND -1))

    # Check for mandatory options
    if [[ -z $aln ]] || [[ -z $control ]] || [[ -z $ucsc_reference ]] || [[ -z $aln_cross ]]; then
        usage
    fi

    # Check if length of arguments in aln1, aln2 and exp are the same
    array_aln=(${aln//[,| ]/ })
    array_control=(${control//[,| ]/ })
    array_aln_cross=(${aln_cross//[,| ]/ })
    if [[(${#array_aln[@]} -ne ${#array_control[@]}) && (${#array_control[@]} -ne ${#array_aln[@]}) && (${#array_aln[@]} -ne ${#array_aln_cross[@]})]]; then
        echo "The number of arguments are not equal."
        exit 1
    fi

    # Define the genome size to use
    if [ $ucsc_reference = 'hg19' ]; then
        genome_size='hs'
    elif [ $ucsc_reference = 'mm10' ]; then
        genome_size='mm'
    elif [ $ucsc_reference = 'mm9' ]; then
        genome_size='mm'
    else
        usage
    fi

    # Define the output directory, if none defined make the location relative to first file
    if [ -z $out ]; then
        one_parent=$(dirname "${array_aln[0]}")
        out_dir=$(dirname "${one_parent}")\/$script_name-$script_ver/
    else
        out_dir=$out\/$script_name-$script_ver
    fi

    if [ ! -d $out_dir ]; then
        mkdir $out_dir
    fi

    # Align if file doesn't exist
    if [ ! -e $out_dir/metadata.json ]; then

        # split files
        rep1=${array_aln[0]}
        rep2=${array_aln[1]}
        con1=${array_control[0]}
        con2=${array_control[1]}
        rep1_xcor=${array_aln_cross[0]}
        rep2_xcor=${array_aln_cross[1]}

        # Count tags
        ntags_rep1=`gunzip $rep1 -c | wc -l | cut -f1`
        ntags_rep2=`gunzip $rep2 -c | wc -l | cut -f1`
        ntags_ctl1=`gunzip $con1 -c | wc -l | cut -f1`
        ntags_ctl2=`gunzip $con2 -c | wc -l | cut -f1`

        # pool all experiment and controls into single files
        gzip -dc $aln | gzip -c > $out_dir/experiment_pooled.tagAlign.gz
        pooled_replicates=`ls $out_dir/experiment_pooled.tagAlign.gz`
        gzip -dc $control | gzip -c > $out_dir/control_pooled.tagAlign.gz
        pooled_controls=`ls $out_dir/control_pooled.tagAlign.gz`

        # Use the pooled controls for the reps depending on the ratio of rep to control reads
        ratio_ctl_reads=`echo $ntags_ctl1/$ntags_ctl2 | bc -l`
        if [ $(echo " $ratio_ctl_reads < 1" | bc) -eq 1 ]; then
                ratio_ctl_reads=`echo 1/$ratio_ctl_reads | bc -l`
        fi

        ratio_cutoff=1.2
        if [ $(echo " $ratio_ctl_reads > $ratio_cutoff" | bc) -eq 1 ]; then
                echo "Number of reads in controls differ by > factor of $ratio_cutoff. Using pooled controls."
                con1=$pooled_controls
                con2=$pooled_controls
        else
                if [ $(echo " $ntags_ctl1 < $ntags_rep1" | bc) -eq 1 ]; then
                        echo "Fewer reads in control replicate 1 than experiment replicate 1.  Using pooled controls for replicate 1."
                        con1=$pooled_controls
                elif [ $(echo " $ntags_ctl2 < $ntags_rep2" | bc) -eq 1 ]; then
                        echo "Fewer reads in control replicate 2 than experiment replicate 2.  Using pooled controls for replicate 2."
                        con2=$pooled_controls
                else
                        echo "Using distinct controls for replicate 1 and 2."
                fi
        fi

        # Make Psudeo replicates
        (exec ./pseudoreplicator.sh -f $pooled_replicates -o $out_dir)
        wait

        # Fetch chrom sizes
        echo "* Downloading chrom.sizes..."
        fetchChromSizes $ucsc_reference > $out_dir/chrom.sizes
        chrom_sizes=$out_dir/chrom.sizes

        # Call peaks Replicates 1 and 2
        call_peak $rep1 $con1 $rep1_xcor $chrom_sizes $genome_size $out_dir

        call_peak $rep2 $con2 $rep2_xcor $chrom_sizes $genome_size $out_dir

        # Call Peaks Pooled replicates
        base_fn=$(basename "${pooled_replicates}")
        output_fp=${base_fn%.tagAlign.gz}
        (exec ./cross-correlation-only-macs.sh -f $pooled_replicates -o $out_dir)
        wait
        pooled_xcor=`ls $out_dir/cross-correlation-only-macs.sh-1.0.0/$output_fp*.cc.qc`
        call_peak $pooled_replicates $pooled_controls $pooled_xcor $chrom_sizes $genome_size $out_dir

        # Pseudoreplicates

        #Pooled Replicated
        #Psuedo Replicate 1
        base_fn=$(basename "${pooled_replicates}")
        output_fp=${base_fn%.tagAlign.gz}
        repppr1=`ls $out_dir/pseudoreplicator.sh-1.0.0/$output_fp.pr1.tagAlign.gz`
        (exec ./cross-correlation-only-macs.sh -f $repppr1 -o $out_dir)
        wait
        repppr1_xcor=`ls $out_dir/cross-correlation-only-macs.sh-1.0.0/$output_fp.pr1*.cc.qc`
        call_peak $repppr1 $pooled_controls $repppr1_xcor $chrom_sizes $genome_size $out_dir

        #Psuedo Replicate 2
        base_fn=$(basename "${pooled_replicates}")
        output_fp=${base_fn%.tagAlign.gz}
        repppr2=`ls $out_dir/pseudoreplicator.sh-1.0.0/$output_fp.pr2.tagAlign.gz`
        (exec ./cross-correlation-only-macs.sh -f $repppr2 -o $out_dir)
        wait
        repppr2_xcor=`ls $out_dir/cross-correlation-only-macs.sh-1.0.0/$output_fp.pr2*.cc.qc`
        call_peak $repppr2 $pooled_controls $repppr2_xcor $chrom_sizes $genome_size $out_dir

        # Remove chrom sizes
        rm $out_dir/chrom.sizes

        # Get input and output files and then print out metadata.json file
        input_files=("${array_aln[@]}" "${array_control[@]}" "${array_aln_cross[@]}" "$ucsc_reference")
        printf -v input "\"%s\"," "${input_files[@]}"
        input=${input%,}
        output_file=($out_dir\/*)
        printf -v output "\"%s\"," "${output_file[@]}"
        output=${output%,}
        printf '{"script name":"%s","script version":"%s", "input files": [%s], "output files": [%s]}' "$script_name" "$script_ver" "$input" "$output"  | python -m json.tool > $out_dir/metadata.json
    else
        echo "* Peaks have been called"
    fi
}

main "$@"
