#!/bin/bash

#'''Usage: ./annotate.sh
#Needs: SV_Simulation.py Iteration#.bed output files in same folder and downloaded annotations in ../../downloaded_annotations'''

for SAMPLE in $(seq 100)
do
  #bedtools intersect -wa -wb -a Iteration$SAMPLE.bed -b ../../downloaded_annotations/Homo_sapiens.GRCh38.96.gtf > Iteration_$SAMPLE.geneAnnotation.bed
  #bedtools intersect -wa -wb -a Iteration$SAMPLE.bed -b ../../downloaded_annotations/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20180918.gff > Iteration_$SAMPLE.regulatoryAnnotation.bed
  bedtools intersect -wa -wb -a Iteration$SAMPLE.bed -b ../../downloaded_annotations/E027_25_imputed12marks_segments.bed.rmCHR.bed > Iteration_$SAMPLE.myoepithelial_chromHMM.bed
  bedtools intersect -wa -wb -a Iteration$SAMPLE.bed -b ../../downloaded_annotations/E028_25_imputed12marks_segments.bed.rmCHR.bed > Iteration_$SAMPLE.vHMEC_chromHMM.bed
  bedtools intersect -wa -wb -a Iteration$SAMPLE.bed -b ../../downloaded_annotations/E119_25_imputed12marks_segments.bed.rmCHR.bed > Iteration_$SAMPLE.HMEC_chromHMM.bed
  echo annotated $SAMPLE
done
