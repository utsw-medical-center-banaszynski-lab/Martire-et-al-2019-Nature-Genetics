#!/bin/bash
#SBATCH --job-name=profile-heatmap                              # job name
#SBATCH --partition=super                                 # select partion from 128GB, 256GB, 384GB, GPU and super
#SBATCH --nodes=1                                         # number of nodes requested by user
#SBATCH --time=0-24:00:00                                 # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=job_prof_HM.%j.out                         # standard output file name
#SBATCH --error=job_prof_HM.%j.err                         # standard error output file name
#SBATCH --mail-user=email@utsouthwestern.edu           # specify an email address
#SBATCH --mail-type=ALL                                   # send email when job status change (start, end, abortion and etc.)

module load deeptools

#WT and H3.3 KO H3K27ac under all ESC Enhancers
computeMatrix reference-point --referencePoint center -R ESC_enhancers.bed -b 3000 -a 3000 --sortRegions no -S /path-to-biwig/WT_H3K27ac.bw /path-to-biwig/KO_H3K27ac.bw --skipZeros -o matrix_1.gz --outFileSortedRegions matrix_1.bed --numberOfProcessors max
plotHeatmap -m matrix_1.gz -out WT_H33KO_H3K27ac_under_ESCenhancers_heatmap.pdf --heatmapHeight 10 --refPointLabel center --samplesLabel 129_K27ac 282_K27ac --colorMap Blues --plotTitle 'WT & H3.3 KO H3K27ac under ESC Enhancers' --dpi 300 --whatToShow 'heatmap and colorbar' --boxAroundHeatmaps no
plotProfile -m matrix_1.gz -out WT_H33KO_H3K27ac_under_ESCenhancers_profile.pdf --perGroup --color blue red --samplesLabel 129_K27ac 282_K27ac --startLabel center --plotTitle "WT & H3.3 KO H3K27ac under ESC Enhancers"


#WT and H3.3 KO H3K27ac under Refseq genes
computeMatrix reference-point --referencePoint TSS -R mm10_RefSeqgenes.bed -b 3000 -a 3000 --sortRegions no -S /path-to-biwig/WT_H3K27ac.bw /path-to-biwig/KO_H3K27ac.bw --skipZeros -o matrix_2.gz --outFileSortedRegions matrix_2.bed --numberOfProcessors max
plotHeatmap -m matrix_2.gz -out WT_H33KO_H3K27ac_under_RefseqGenes_heatmap.pdf --heatmapHeight 10 --refPointLabel center --samplesLabel 129_K27ac 282_K27ac --colorMap Blues --plotTitle 'WT & H3.3 KO H3K27ac under Refseq genes' --dpi 300 --whatToShow 'heatmap and colorbar' --boxAroundHeatmaps no
plotProfile -m matrix_2.gz -out WT_H33KO_H3K27ac_under_RefseqGenes_profile.pdf --perGroup --color blue red --samplesLabel 129_K27ac 282_K27ac --startLabel TSS --plotTitle "WT & H3.3 KO H3K27ac under Refseq genes"
