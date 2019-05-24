#!/bin/bash

#bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/N_only_51_withTRA.bed -b ../downloaded_annotations/Homo_sapiens.GRCh38.96.gtf > 51N_only_withTRA_geneAnnotation.bed
#bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/N_only_51_withTRA.bed -b ../downloaded_annotations/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20180918.gff > 51N_only_withTRA_regulatoryAnnotation.bed

bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/N_only_51_withTRA.bed -b ../downloaded_annotations/E027_25_imputed12marks_segments.bed.rmCHR.bed > N_only.myoepithelial_chromHMM.bed
bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/N_only_51_withTRA.bed -b ../downloaded_annotations/E028_25_imputed12marks_segments.bed.rmCHR.bed > N_only.vHMEC_chromHMM.bed
bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/N_only_51_withTRA.bed -b ../downloaded_annotations/E119_25_imputed12marks_segments.bed.rmCHR.bed > N_only.HMEC_chromHMM.bed

#bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/T_only_51_withTRA.bed -b ../downloaded_annotations/Homo_sapiens.GRCh38.96.gtf > 51T_only_withTRA_geneAnnotation.bed
#bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/T_only_51_withTRA.bed -b ../downloaded_annotations/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20180918.gff > 51T_only_withTRA_regulatoryAnnotation.bed

bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/T_only_51_withTRA.bed -b ../downloaded_annotations/E027_25_imputed12marks_segments.bed.rmCHR.bed > T_only.myoepithelial_chromHMM.bed
bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/T_only_51_withTRA.bed -b ../downloaded_annotations/E028_25_imputed12marks_segments.bed.rmCHR.bed > T_only.vHMEC_chromHMM.bed
bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/T_only_51_withTRA.bed -b ../downloaded_annotations/E119_25_imputed12marks_segments.bed.rmCHR.bed > T_only.HMEC_chromHMM.bed

#bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/both_51_withTRA.bed -b ../downloaded_annotations/Homo_sapiens.GRCh38.96.gtf > 51Both_withTRA_geneAnnotation.bed
#bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/both_51_withTRA.bed -b ../downloaded_annotations/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20180918.gff > 51Both_withTRA_regulatoryAnnotation.bed

bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/both_51_withTRA.bed -b ../downloaded_annotations/E027_25_imputed12marks_segments.bed.rmCHR.bed > both.myoepithelial_chromHMM.bed
bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/both_51_withTRA.bed -b ../downloaded_annotations/E028_25_imputed12marks_segments.bed.rmCHR.bed > both.vHMEC_chromHMM.bed
bedtools intersect -wa -wb -a ../bedFiles_forBedtools_withIDs/both_51_withTRA.bed -b ../downloaded_annotations/E119_25_imputed12marks_segments.bed.rmCHR.bed > both.HMEC_chromHMM.bed
