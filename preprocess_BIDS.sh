#!/bin/sh 
#usage: bash preproc subjDIR outputDIR CO2 time [params]
VDIR=$PWD
DDIR=${VDIR}/$1
cd ${DDIR}/func/
#motion correction
mcflirt -in *bold.nii* -out BOLD_mcf -spline_final -bins 256 -cost leastsquares -dof 6 -plots -report
#calculate mean
fslmaths BOLD_mcf -Tmean BOLD_mean
bet BOLD_mean BOLD_mean_brain -R -f 0.30 -g 0 -m

#register mean bold to T1
cd ${VDIR}/$1/anat/
bet *T1w.nii* T1_brain.nii.gz -R -f 0.30 -g 0 -m
echo completed brain extraction of T1 image, moving on to segmentation
fast -t 1 -n 3 -H 0.1 -I 3 -b -l 20.0 -g --nopve -o T1_brain.nii
echo completed segmentation of T1, moving on to registration of mean BOLD to T1

cd ${DDIR}/func/
epi_reg --noclean -v --epi=BOLD_mean_brain --t1=${DDIR}/anat/*T1w.nii* --t1brain=${DDIR}/anat/T1_brain.nii.gz --wmseg=${DDIR}/anat/T1_brain_seg_2.nii.gz --out=BOLDtoT1.nii.gz
echo completed registration of BOLD to T1, inverting T1 segmentations
convert_xfm -omat ${DDIR}/func/BOLDtoT1_inv.mat -inverse ${DDIR}/func/BOLDtoT1.mat
echo applying inverse transformations to T1 segmentations
echo CSF seg
flirt -in ${DDIR}/anat/T1_brain_seg_0.nii.gz -applyxfm -init ${DDIR}/func/BOLDtoT1_inv.mat -out ${DDIR}/func/brain_seg_0.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref ${DDIR}/func/BOLD_mean.nii.gz
echo GM seg
flirt -in ${DDIR}/anat/T1_brain_seg_1.nii.gz -applyxfm -init ${DDIR}/func/BOLDtoT1_inv.mat -out ${DDIR}/func/brain_seg_1.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref ${DDIR}/func/BOLD_mean.nii.gz
echo WM seg
flirt -in ${DDIR}/anat/T1_brain_seg_2.nii.gz -applyxfm -init ${DDIR}/func/BOLDtoT1_inv.mat -out ${DDIR}/func/brain_seg_2.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref ${DDIR}/func/BOLD_mean.nii.gz
echo brain mask
flirt -in ${DDIR}/anat/T1_brain_mask.nii.gz -applyxfm -init ${DDIR}/func/BOLDtoT1_inv.mat -out ${DDIR}/func/brain_mask.nii.gz -paddingsize 0.0 -interp nearestneighbour -ref ${DDIR}/func/BOLD_mean.nii.gz


# Run Gas_CVR
echo "Running Gas_CVR analysis"
OUTDIR=$2
CO2=$3
TIME=$4
shift; shift; shift; shift

# Prepare data directory
echo "Preparing output directory..."
cd ${VDIR}

mkdir ${OUTDIR}
cp ${DDIR}/func/BOLD_mcf.nii.gz ${OUTDIR}/BOLD_mcf.nii.gz
cp ${DDIR}/func/BOLD_mcf.par ${OUTDIR}/BOLD_mcf.par
cp ${DDIR}/func/BOLD_mean_brain_mask.nii.gz ${OUTDIR}/BOLD_mean_brain_mask.nii.gz
cp ${DDIR}/func/brain_seg_0.nii.gz ${OUTDIR}/brain_seg_0.nii.gz
cp ${DDIR}/func/brain_seg_1.nii.gz ${OUTDIR}/brain_seg_1.nii.gz
cp ${DDIR}/func/brain_seg_2.nii.gz ${OUTDIR}/brain_seg_2.nii.gz

cp ${CO2} ${OUTDIR}/CO2.txt
cp ${TIME} ${OUTDIR}/time.txt

echo "Running gas_cvr..."

COMMAND="./run_gas_cvr.sh /opt/MCR-2021b/v911/ \
  ${OUTDIR} \
  time.txt \
  CO2.txt \
  BOLD_mcf.nii.gz \
  brain_seg_0.nii.gz \
  brain_seg_1.nii.gz \
  brain_seg_2.nii.gz \
  BOLD_mean_brain_mask.nii.gz \
  BOLD_mcf.par \
  "$@""

echo "Command is: ${COMMAND}"
eval ${COMMAND}

echo "Task completed."