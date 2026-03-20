#!/usr/bin/env python

# do combine planes for arbitary numbers of scans with ants tools only, ready
# for super resolution

# plan should be for a maskless combination first to allow a better mask
# to be computed, then rerun

import argparse
import sys
import re
import os.path
import SimpleITK as sitk
from tempfile import mktemp
#import local_ants
import numpy as np

def is_readable_file(parser, arg):
    try:
        f = open(arg, 'r')
        f.close()
    except Exception:
        raise argparse.ArgumentTypeError("{0} does not exist or is not readable".format(arg))

    return(arg)


parser = argparse.ArgumentParser(description="Preprocessing of images for superresolution")

parser.add_argument("-i", "--images", nargs="+",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="One or more input nifti files")

parser.add_argument("-w", "--workdir",
                    required=True,
                    help="folder where temp files are written")

parser.add_argument("-s", "--spacing",
                    required=False,
                    default=None,
                    type=float,
                    help= "Target spacing - otherwise minimum of inputs")

# type of transform for final reg to template
#TFTYPE="Rigid"
TFTYPE="Affine"

args = parser.parse_args()

workdir = args.workdir
workdir = os.path.join(workdir, "")

antsworkdir = os.path.join(workdir, "antsstuff", "")

os.makedirs(antsworkdir, exist_ok = True)
    

####################
# preprocess with simpleitk - N4 (without mask) + intensity standardize + histogram normaliztion
preprocnames = [ re.sub("\\.nii\\.gz", "_mask.nii.gz", s) for s in args.images ]
preprocnames = [ os.path.join(workdir, os.path.basename(x)) for x in preprocnames ]

wnames = [ re.sub("\\.nii\\.gz", "_w.nii.gz", s) for s in args.images ]
wnames = [ os.path.join(workdir, os.path.basename(x)) for x in wnames ]

debugants = [ re.sub("\\.nii\\.gz", "_dbg.nii.gz", s) for s in args.images ]
debugants = [ os.path.join(workdir, os.path.basename(x)) for x in debugants ]

import sitkPreprocess as sP

METRIC="GC"
RS = 334
####################
# some of this stuff could be done with antspy
# I started writing it before playing with ants
# probably still need sitk for transforms at least,
# and there's a bigger choice of filters.

images =  [ sitk.ReadImage(x, sitk.sitkFloat32) for x in args.images ]

# get the voxel spacing to figure out our final target space
if args.spacing is None:
    spacings = [min(x.GetSpacing()) for x in images]
    spacingtarget = min(spacings)
    spacingtarget = (spacingtarget, ) * 3
else:
    spacingtarget = (args.spacing, ) * 3

images = None

## Now for ants registration
## aim is to register everything to first image in the list
import ants
import ants_utils
# haven't dealt with potential cropping issues.
# simplest way for now is to make sure that the two most
# similar images are the first two in the list

antsimages = [ ants.image_read(x) for x in args.images ]
# remove negative values, which sometimes mess up registration
antsimages = [ x * (x>0) for x in antsimages]
# Note - antsregistration has outprefix for the transforms, but it does need to be unique within an
# application because there's some sort of file listing that occurs in order to return the filenames.
# Thus there are problems when two registrations are run with the same outprefix 

def mkAntsTmp():
    # this isn't particularly nice. Probably OK if we
    # recreate the workdir every time
    return os.path.join(antsworkdir, mktemp(prefix=antsworkdir))

print("Starting Registration", flush=True)

regresults = [ ants.registration(fixed=antsimages[0], moving=antsimages[idx+1],
                                 type_of_transform=TFTYPE, outprefix=mkAntsTmp(), verbose=True, aff_random_sampling_rate=1.0, aff_metric = METRIC,
                                 random_seed = RS)
               for idx in range(len(antsimages)-1)  ]

# Perhaps I should include a shape update here?
print("Regresults")
print(regresults)
TF01 = sitk.ReadTransform(regresults[0]['fwdtransforms'][0])
halfway = sP.sitkHalfway(TF01)
hname = os.path.dirname(regresults[0]['fwdtransforms'][0])
hname = os.path.join(hname, "halfway.mat")
sP.sitkWriteAffine(halfway, hname, np.array(TF01.GetCenter()))

# create a high res sample image - note that this is to control the resampling, not part
# of registration
HighRes = ants.resample_image(antsimages[0], spacingtarget, interp_type = 0)

#W0 = ants.apply_transforms(fixed=antsimages[0], moving=antsimages[1], transformlist = regresults[0]['fwdtransforms'][0], whichtoinvert=[False])
W1 = ants.apply_transforms(fixed=antsimages[0], moving=antsimages[1], transformlist = hname, whichtoinvert=[False])

tflist = [ reg['fwdtransforms'][0] for reg in regresults]
print(tflist)

# No transform composition required
W0 = ants.apply_transforms(fixed=HighRes, moving=antsimages[0], transformlist = hname, whichtoinvert=[True])
ants.image_write(W0, wnames[0])

Accumulated = W0.clone()
for idx in range(len(regresults)):
    fixed=HighRes
    moving=antsimages[idx+1]

    tlist=[hname, tflist[idx]]
    invlist=[True, False]
    W = ants.apply_transforms(fixed=fixed, moving=moving, transformlist = tlist, whichtoinvert=invlist)
    ants.image_write(W, wnames[idx+1])
    Accumulated += W

Accumulated /= float(len(antsimages))
ants.image_write(Accumulated, os.path.join(workdir, "avstep0.nii.gz"))

import functools

# make a mask and transform back
hires = sitk.ReadImage(os.path.join(workdir, "avstep0.nii.gz"))

# basic image processing - almost certainly
# too simple
# small noise reduction opening first
hiresO = sitk.GrayscaleMorphologicalOpening(hires, kernelRadius = (3,3,3), kernelType = sitk.sitkBox)
#sitk.WriteImage(hiresO, "opened.nii.gz")

# Now thresholding - we need to ignore the zero'ed
# background because it messes with distributions
nonzero = hiresO > 0
otsuF = sitk.OtsuThresholdImageFilter()
otsuF.SetInsideValue(0)
otsuF.SetOutsideValue(1)
otsuF.SetMaskValue(1)
hiresM = otsuF.Execute(hiresO, nonzero)
#sitk.WriteImage(hiresM, "hiresM.nii.gz")

radius =  (5,5,5)
hiresM = sitk.BinaryDilate(hiresM, radius, sitk.sitkBox)
hiresM = sitk.BinaryFillhole(hiresM)
hiresM = sitk.BinaryErode(hiresM, radius, sitk.sitkBox)
# keep biggest
sitk.WriteImage(hiresM, os.path.join(workdir, "mask_template.nii.gz"))

# transform the masks back to inputs
mask = ants.image_read(os.path.join(workdir, "mask_template.nii.gz"))

mask0 = ants.apply_transforms(fixed=antsimages[0], moving=mask, transformlist = hname, whichtoinvert=[False])
ants.image_write(mask0, preprocnames[0])

for idx in range(len(regresults)):
    moving=mask
    fixed=antsimages[idx+1]

    tlist=[hname, tflist[idx]]
    invlist=[False, True]
    W = ants.apply_transforms(fixed=fixed, moving=moving, transformlist = tlist, interpolator='nearestNeighbor', whichtoinvert=invlist)
    ants.image_write(W, preprocnames[idx+1])
