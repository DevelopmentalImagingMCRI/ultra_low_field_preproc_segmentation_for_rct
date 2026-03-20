#!/bin/bash
#SBATCH --job-name=recon
#SBATCH --output=**REPLACE** path to /Logs/%x-%j-%a.out

# run a single restoration/reconstruction
hostname
date
echo $SLURM_CPUS_PER_TASK
# arguments are a a set of files and a destination folder and prefix
SANDBOX=**REPLACE** path to niftymic.simg
# **REPLACE** make a python environment manager available
module purge
module load miniconda3/24.7.1-q7s6uee

# activate the conda environment
source  activate **REPLACE** path to /Scripts/PySuperRes2

# Lets assume that the argument is a space separated list of arguments.
# slurm array jobs will provide an index into the file
# the first entry will be destination folder
ARGLIST=${1}
FILEIDX=${SLURM_ARRAY_TASK_ID:-1}

ARGSARRAY=( $(sed "${FILEIDX}q;d" ${ARGLIST}) )

TARGETFOLDER=${ARGSARRAY[0]}
SUBID=${ARGSARRAY[1]}
SESSION=${ARGSARRAY[2]}
FAST_ACQUISITION=${ARGSARRAY[3]}

spack unload -a

ALPHA=${2:-0.001}

# note that we won't support denoise and unwarp together
# Don't unwarp by default
UNWARP=${3:-0}

# optionally modify TARGETFOLDER to allow parameter sweeps
TARGETFOLDER=${TARGETFOLDER}/estimated_psf_${ALPHA}

if [ ${UNWARP} -gt 0 ]
   then
       # test denoising
       TARGETFOLDER=${TARGETFOLDER}_unwarp/
fi


mkdir -p ${TARGETFOLDER}
unset ARGSARRAY[0]
unset ARGSARRAY[1]
unset ARGSARRAY[2]
unset ARGSARRAY[3]
# to deal with the empty array element
ARGSARRAY=( ${ARGSARRAY[*]} )

P=sub-${SUBID}_ses-${SESSION}
if [ -e ${TARGETFOLDER}/${P}_acq-isoMIC_T2w.nii.gz ] ; then
    echo "Already run"
    exit
fi

#WDIR=$(mktemp -d  --tmpdir=/scratch/) || exit 1
WDIR=$(mktemp -d ) || exit 1

trap 'rm -rf $WDIR; exit 0' 0 1 2 3 14 15

export TMPDIR=${WDIR}

SCRIPTDIR=**REPLACE** path to /OSF/Scripts/

MASKDIR=${WDIR}/Masks
COMBDIR=${WDIR}/Comb
COMBDIR2=${WDIR}/Comb2

UNWARPDIR=${WDIR}/Unwarp
THRESHDIR=${WDIR}/LowThresh
mkdir -p ${MASKDIR} ${COMBDIR} ${COMBDIR2} ${UNWARPDIR} ${THRESHDIR}

# setting up names of images at various stages
declare -a MASKARRAY
declare -a PREPROC_MASK_ARRAY
declare -a PREPROC_IMAGE_ARRAY
declare -a UNWARP_IMAGE_ARRAY
declare -a WARP_IMAGE_ARRAY
declare -a LOWTHRESH_IMAGE_ARRAY
IMAGES=${#ARGSARRAY[@]}

for (( j=0;j<${IMAGES};j++ )); do
    B=$(basename ${ARGSARRAY[$j]})
    M=${B/.nii.gz/_mask.nii.gz}
    UW=${B/.nii.gz/_uw.nii.gz}
    PMZ=${B/.nii.gz/_mask_ppz.nii.gz}
    PIZ=${B/.nii.gz/_ppz.nii.gz}
    M=${MASKDIR}/${M}
    MASKARRAY[$j]=$M
    PREPROC_MASK_ARRAY[$j]=/GIFT/${PMZ}
    PREPROC_IMAGE_ARRAY[$j]=/GIFT/${PIZ}
    UNWARP_IMAGE_ARRAY[$j]=${COMBDIR}/${UW}
    WARP_IMAGE_ARRAY[$j]=${UNWARPDIR}/Midway/${B}
    LOWTHRESH_IMAGE_ARRAY[$j]=${THRESHDIR}/${B}
done

# Could set this based on FAST_ACQUISITION
ISORES=1.5

# Set up covariance for deconvolution
# Always the two in-plane FWHM then the throughplane
# These are the values direct from my psf experiments
# Set up covariance for deconvolution
if [ ${FAST_ACQUISITION} = "TRUE" ] ; then
   echo "Using settings for FAST acquisition"
   COV_AXI="0.678 0.678 3.81"
   COV_COR="0.686 0.686 3.11"
   COV_SAG="0.734 0.734 4.93"
else
   echo "Using settings for original acquisitions"
   COV_AXI="0.549 0.549 4.31"
   COV_COR="0.5 0.5 3.37"
   COV_SAG="0.5 0.5 4.24"
fi


if [ ${UNWARP} -gt 0 ]
then
    echo "Unwarping"
    AXUW=**REPLACE** path to /MinihybridPhantomUnWarpBam/T2w_fast_ax.nii.gz
    CORUW=**REPLACE** path to /MinihybridPhantomUnWarpBam/T2w_fast_cor.nii.gz
    SAGUW=**REPLACE** path to /MinihybridPhantomUnWarpBam/T2w_fast_sag.nii.gz

    # **REPLACE**
    module load ants/2.5.1-773ggwn

    IM=${ARGSARRAY[0]}
    OP=$(basename ${IM})
    OP1=${UNWARPDIR}/${OP}
    antsApplyTransforms -d 3 -f 0 -i ${IM} -r ${IM} -o ${OP1} -t ${AXUW}

    IM=${ARGSARRAY[1]}
    OP=$(basename ${IM})
    OP2=${UNWARPDIR}/${OP}
    antsApplyTransforms -d 3 -f 0 -i ${IM} -r ${IM} -o ${OP2} -t ${CORUW}

    IM=${ARGSARRAY[2]}
    OP=$(basename ${IM})
    OP3=${UNWARPDIR}/${OP}
    antsApplyTransforms -d 3 -f 0 -i ${IM} -r ${IM} -o ${OP3} -t ${SAGUW}

    CLEANFILES=( ${OP1} ${OP2} ${OP3} )
    module unload ants
fi

echo ${CLEANFILES[*]}
export PYTHONUNBUFFERED=1
date
# Unwarp first - needs masks - we end up running it twice - hopefully the overwrite is OK.
${SCRIPTDIR}/mkMasks2.py --workdir ${MASKDIR} --images ${CLEANFILES[*]} --spacing ${ISORES}

date
echo "Combine Planes"
# if no target is provided the halfway space is used
${SCRIPTDIR}/combinePlanes2.py --workdir ${COMBDIR} --images ${CLEANFILES[*]} --masks ${MASKARRAY[*]} --spacing ${ISORES}
date


# at the end of this we have an ants template in $COMBDIR/build_template.nii.gz
# as well as a set of images that have been scaled and arranged so that the
# third axis is the thick slice.
# **REPLACE
module load apptainer/1.1.9-pjwiqeh

# multiple runs with different templates
# this line for custom covariance

apptainer run --bind ${COMBDIR}:/GIFT/ ${SANDBOX} niftymic_reconstruct_volume_RB --filenames ${PREPROC_IMAGE_ARRAY[*]}  \
--reference /GIFT/build_template.nii.gz --reference-mask /GIFT/build_template_mask.nii.gz  \
--output=/GIFT/iso_psfA.nii.gz --isotropic-resolution 1.5 --alpha ${ALPHA}  --v2v-reg-typeB Affine --v2v-reg-typeA Affine --v2v-method FLIRT --alpha-first 1 --predefined-covariance-axi ${COV_AXI} --predefined-covariance-cor ${COV_COR} --predefined-covariance-sag ${COV_SAG} 

#
# Move the results
mkdir -p ${TARGETFOLDER}
P=sub-${SUBID}_ses-${SESSION}
mv ${WDIR}/Comb/avstep0.nii.gz ${TARGETFOLDER}/${P}_acq-isoHalf_T2w.nii.gz
mv ${WDIR}/Comb/build_template.nii.gz ${TARGETFOLDER}/${P}_acq-isoANTS_T2w.nii.gz
mv ${WDIR}/Comb/build_template_nl.nii.gz ${TARGETFOLDER}/${P}_acq-isoANTSNL_T2w.nii.gz
mv ${WDIR}/Comb/nl_mrr_axial_sharpened.nii.gz ${TARGETFOLDER}/${P}_acq-isoAxMRRsharp_T2w.nii.gz
mv ${WDIR}/Comb/nl_mrr_axial_notsharpened.nii.gz ${TARGETFOLDER}/${P}_acq-isoAxMRRnosharp_T2w.nii.gz
#mv ${WDIR}/Comb/iso_psfB.nii.gz ${TARGETFOLDER}/${P}_acq-isoMIC_T2w.nii.gz
mv ${WDIR}/Comb/iso_psfA.nii.gz ${TARGETFOLDER}/${P}_acq-isoMIC_T2w.nii.gz
#mv ${WDIR} ${TARGETFOLDER}/${P}_tmp
