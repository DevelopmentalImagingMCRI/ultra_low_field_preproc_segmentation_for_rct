# Methods for preprocessing and segmentation of 64mT MRI of infants

This folder contains scripts for running the super-resolution reconstruction
of Hyperfine scans of babies and infants using a modified version of the
NiftyMIC tools.

`recon_one_slurm.sh` is the script to perform preprocessing, unwarping and super-resolution
reconstruction. It is not yet in a portable form. Below are setup, customization and execution instructions.

Parts of the script requiring customization are marked by the text `**REPLACE**`

Most are paths that are system dependent.


## Prerequisites

### Build the modified NiftyMIC docker and convert to apptainer/singularity

```
cd ../Deconvolution
docker build -t niftymic . 
# convert to apptainer/singularity - replace "singularity" with "apptainer" if
# the latter is installed on your system
# step 1 - save the image archive (could push to dockerhub instead)
docker save $(docker image ls niftymic -q) -o ./niftymic.docker
singularity build niftymic.simg docker-archive://niftymic.docker
```

Place `niftymic.simg` in a convenient location and ensure that the SANDBOX variable in the
slurm script matches it

### ANTS

The slurm script uses `ants` registration tools, assumed to be available via a `module load` command. 

### Custom python environment

```
# creates an environment named PySuperRes2
conda env create -f environment.yml
```

## Execution

The script is structured to be submitted as an array job using the following syntax:
```
batch --time=4:00:00 --partition=prod_med --array=1-$(wc -l < command_file) --cpus-per-task=4 --
mem=16G   ./recon_one_slurm_mk4.sh command_file 0.01 0 
```
Where, the `0.01` value corresponds to the regularization in deconvolution and the following `0` indicates no unwarping.

Use 1 if to include unwarping, which requires an estimate of the unwarp fields (see ...)

The `command_file` has the following structure:

```
path/to/target/folder SUBJECT_ID SESSIONID TRUE/FALSE /path/to/axial_scan.nii.gz /path/to/coronal_scan.nii.gz /path/to/sagittal_scan.nii.gz
```

Where the TRUE/FALSE value is to indicate whether a FAST style acquisition was used. Different point spread function 
estimates and warp fields will be selected based on this value.

The outputs will be written to the target folder in approximately BIDS format, assuming that the input was approximately in BIDS.
