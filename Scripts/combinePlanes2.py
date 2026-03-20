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
import ants
import local_ants
import ants_utils
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

parser.add_argument("-m", "--masks", nargs="+",
                    type=lambda x: is_readable_file(parser, x),
                    required=True,
                    help="One or more input nifti mask files (same length as iamges")

parser.add_argument("-t", "--target", nargs=1,
                    type=lambda x: is_readable_file(parser, x),
                    required=False,
                    help="A high res target image. One is created if not supplied")

parser.add_argument("-s", "--spacing",
                    required=False,
                    default=None,
                    type=float,
                    help= "Target spacing - otherwise minimum of inputs")

# type of transform for final reg to template
#TFTYPE="Rigid"
TFTYPE="Affine"
METRIC="GC"
AFF_SAMPLING=1.0
# Random seed
RS = 334

args = parser.parse_args()
if len(args.images) != len(args.masks):
    print("Masks and images must match")
    exit()


workdir = args.workdir
workdir = os.path.join(workdir, "")

antsworkdir = os.path.join(workdir, "antsstuff", "")

os.makedirs(antsworkdir, exist_ok = True)
    
print(workdir)
#if len(args.images) != len(args.out):
#    print("Out and images must match")
#    exit()


####################
# preprocess with simpleitk - N4 (without mask) + intensity standardize + histogram normaliztion
preprocnames = [ re.sub("\\.nii\\.gz", "_pp.nii.gz", s) for s in args.images ]
preprocnames = [ os.path.join(workdir, os.path.basename(x)) for x in preprocnames ]

# same, but dimensions reordered so that thick slices are dimension 3
preprocnamesZ = [ re.sub("\\.nii\\.gz", "_ppz.nii.gz", s) for s in args.images ]
preprocnamesZ = [ os.path.join(workdir, os.path.basename(x)) for x in preprocnamesZ ]

# masks with the same dimension ordering as the data
masknamesZ = [ re.sub("\\.nii\\.gz", "_mask_ppz.nii.gz", s) for s in args.images ]
masknamesZ = [ os.path.join(workdir, os.path.basename(x)) for x in masknamesZ ]



wnames = [ re.sub("\\.nii\\.gz", "_w.nii.gz", s) for s in args.images ]
wnames = [ os.path.join(workdir, os.path.basename(x)) for x in wnames ]

debugants = [ re.sub("\\.nii\\.gz", "_dbg.nii.gz", s) for s in args.images ]
debugants = [ os.path.join(workdir, os.path.basename(x)) for x in debugants ]

uwnames =  [ re.sub("\\.nii\\.gz", "_uw.nii.gz", s) for s in args.images ]
uwnames = [ os.path.join(workdir, os.path.basename(x)) for x in uwnames ]

import sitkPreprocess as sP

####################
# some of this stuff could be done with antspy
# I started writing it before playing with ants
# probably still need sitk for transforms at least,
# and there's a bigger choice of filters.

images =  [ sitk.ReadImage(x, sitk.sitkFloat32) for x in args.images ]
masks = [ sitk.ReadImage(x, sitk.sitkUInt8) for x in args.masks ]

# get the voxel spacing to figure out our final target space
if args.spacing is None:
    spacings = [min(x.GetSpacing()) for x in images]
    spacingtarget = min(spacings)
    spacingtarget = (spacingtarget, ) * 3
else:
    spacingtarget = (args.spacing, ) * 3

# could save these if debug flags are set
imagesN4 = [ sP.N4(*pair) for pair in zip(images, masks) ]

stdimages = sP.IntStd(imagesN4, masks)

normalised = sP.HistNorm(stdimages, masks)
normalised = sP.IntStd2(normalised, masks, targetMean=150, targetSD=50)

[sitk.WriteImage(p[0], p[1]) for p in zip(normalised, preprocnames)]

## permute axes so that the thick slice is the last dimension.
permuted = [sP.thickSliceLast(x) for x in normalised]
[sitk.WriteImage(p[0], p[1]) for p in zip(permuted, preprocnamesZ)]

permutedmask = [ sP.thickSliceLast(x) for x in masks ]
[sitk.WriteImage(p[0], p[1]) for p in zip(permutedmask, masknamesZ)]

normalised = permuted

def createMask(inname, outname):
    # duplication of mkMask - could refactor
    hires = sitk.ReadImage(inname)

    # basic image processing - almost certainly
    # too simple
    # small noise reduction opening first
    hiresO = sitk.GrayscaleMorphologicalOpening(hires, kernelRadius = (3,3,3), kernelType = sitk.sitkBox)

    # Now thresholding - we need to ignore the zero'ed
    # background because it messes with distributions
    nonzero = hiresO > 0
    otsuF = sitk.OtsuThresholdImageFilter()
    otsuF.SetInsideValue(0)
    otsuF.SetOutsideValue(1)
    otsuF.SetMaskValue(1)
    hiresM = otsuF.Execute(hiresO, nonzero)
    
    radius =  (5,5,5)
    hiresM = sitk.BinaryDilate(hiresM, radius, sitk.sitkBox)
    hiresM = sitk.BinaryFillhole(hiresM)
    hiresM = sitk.BinaryErode(hiresM, radius, sitk.sitkBox)
    # keep biggest
    sitk.WriteImage(hiresM, outname)

## Now for ants registration
## aim is to register everything to first image in the list

# Local versions of build_template and possibly registratrion - not
# all the arguments are available
#import local_ants
# haven't dealt with potential cropping issues.
# simplest way for now is to make sure that the two most
# similar images are the first two in the list

antsimages = [ ants.image_read(x) for x in preprocnamesZ ]
print("antsimages")
print(preprocnamesZ)
print(antsimages)

# Note - antsregistration has outprefix for the transforms, but it does need to be unique within an
# application because there's some sort of file listing that occurs in order to return the filenames.
# Thus there are problems when two registrations are run with the same outprefix 

def mkAntsTmp():
    # this isn't particularly nice. Probably OK if we
    # recreate the workdir every time
    return os.path.join(antsworkdir, mktemp(prefix=antsworkdir))

if args.target is None:
    print("Creating target template")
    # we need to create the target template
    regresults = [ ants.registration(fixed=antsimages[0], moving=antsimages[idx+1],
                                 type_of_transform=TFTYPE, outprefix=mkAntsTmp(), 
                                 verbose=True, aff_metric = METRIC,
                                 aff_random_sampling_rate = AFF_SAMPLING,
                                 random_seed = RS)
               for idx in range(len(antsimages)-1)  ]

    # Perhaps I should include a shape update here?
    #print(regresults)
    TF01 = sitk.ReadTransform(regresults[0]['fwdtransforms'][0])
    halfway = sP.sitkHalfway(TF01)
    hname = os.path.dirname(regresults[0]['fwdtransforms'][0])
    hname = os.path.join(hname, "halfway.mat")
    sP.sitkWriteAffine(halfway, hname)

    # create a high res sample image - note that this is to control the resampling, not part
    # of registration. We could do some extra work
    # to make sure we don't miss some edges

    HighRes = ants.resample_image(antsimages[0], spacingtarget, interp_type = 0)

    #W0 = ants.apply_transforms(fixed=antsimages[0], moving=antsimages[1], transformlist = regresults[0]['fwdtransforms'][0], whichtoinvert=[False])
    #W1 = ants.apply_transforms(fixed=antsimages[0], moving=antsimages[1], transformlist = hname, whichtoinvert=[False])

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
        # just compose transforms for later use


    Accumulated /= float(len(antsimages))
    fname = os.path.join(workdir, "avstep0.nii.gz")
    fnameM = os.path.join(workdir, "avstep0_mask.nii.gz")

    ants.image_write(Accumulated, fname)

    createMask(fname, fnameM)
else:
    Accumulated = ants.resample_image(ants.image_read(args.target[0]), spacingtarget, interp_type = 0)
    HighRes = ants.resample_image(antsimages[0], spacingtarget, interp_type = 0)

import functools

# build_template doesn't work with affine/rigid transforms
# I suspect it is to do with the transforms being read back in as images
# There aren't any for affines. Approximate by turning iterations
# down low
print("* Build linear template **")
template = local_ants.build_template(initial_template=Accumulated, image_list=antsimages, 
                               iterations=4, aff_iterations=(1000,500,250,0),
                               aff_shrink_factors = (6,4,2,1),
                               aff_smoothing_sigmas=(4,2,1,0),
                               reg_iterations=[2,2,2,2], useNoRigid=True, 
                               blending_weight = 1, type_of_transform='SyN', verbose = False,
                               winsorize_image_intensities = (0.01,0.99),
                               aff_metric = METRIC, syn_metric = METRIC,
                               random_seed = RS)
#                               type_of_transform=TFTYPE)
#                               outprefix=mkAntsTmp())

fname = os.path.join(workdir, "build_template.nii.gz")
fnameM = os.path.join(workdir, "build_template_mask.nii.gz")
ants.image_write(template,  fname)
# save a mask too
createMask(fname, fnameM)
print("* Finished linear template **")

print("* Build nonlinear template **")
templateNL = local_ants.build_template(initial_template = Accumulated, image_list=antsimages,
                                     iterations=4, aff_iterations=(1000,500,250,0),
                                     aff_shrink_factors = (6,4,2,1),
                                     aff_smoothing_sigmas=(4,2,1,0),
                                     reg_iterations=[100,100,70,20],
                                     grad_step = 0.1,
                                     flow_sigma=3,
                                     total_sigma=0,
                                     useNoRigid=True,
                                     blending_weight = 1,
                                     type_of_transform='SyN',
                                     verbose=True,
                                     aff_metric = METRIC, syn_metric = METRIC,
                                     random_seed = RS)

fname = os.path.join(workdir, "build_template_nl.nii.gz")
fnameM = os.path.join(workdir, "build_template_nl_mask.nii.gz")
ants.image_write(templateNL, fname)
createMask(fname, fnameM)
# finished nonlinear template
def mrrTemplate(imlist):
    print("* Build nonlinear template - pure CISO style**")
    # Finally tracked it down. The largest difference in the appearance of this vs
    # CISO/antsMultivariateTemplate is the default sharpening step, that can't be disabled.
    # The equivalent in python is via blending_weight = 0, which is entirely sharpened
    #
    templateS1 = local_ants.build_template(initial_template = imlist[0], image_list=imlist, 
                               iterations=4, aff_iterations=(1000,500,250,0),
                               reg_iterations=[100,100,70,20],
                               useNoRigid=False, 
                               blending_weight = 0, type_of_transform='SyN', verbose = True,
                               aff_metric = METRIC, syn_metric = METRIC, random_seed = RS)
    templateS1 = ants.resample_image(templateS1, resample_params=spacingtarget, interp_type = 0)

    regresults = [ ants.registration(fixed=templateS1, moving=imlist[idx],
                                     type_of_transform="Rigid", outprefix=mkAntsTmp(), verbose=False,
                                     aff_metric = METRIC, syn_metric = METRIC, random_seed = RS)
                 for idx in range(len(antsimages))  ]
    
    rr = [ x["warpedmovout"] for x in regresults]
    [ants.image_write(im, os.path.join(workdir, "nl_src_" + str(idx) + ".nii.gz")) for idx, im in enumerate(rr)   ]
    templateNL = local_ants.build_template(initial_template = templateS1, image_list=rr,
                                     iterations=4, aff_iterations=(1000,500,250,0),
                                     aff_shrink_factors = (6,4,2,1),
                                     aff_smoothing_sigmas=(4,2,1,0),
                                     reg_iterations=[100,100,70,20],
                                     grad_step = 0.1,
                                     flow_sigma=3,
                                     total_sigma=0,
                                     useNoRigid=True,
                                     blending_weight = 0,
                                     type_of_transform='SyN',
                                     verbose=True,
                                     aff_metric = METRIC, syn_metric = METRIC,
                                     random_seed = RS)
    fnameNL = os.path.join(workdir, "nl_build_template.nii.gz")
    ants.image_write(templateNL,  fnameNL)

def mkCISOMRR():
    images_reloaded = [ants.image_read(x) for x in args.images]
    print("* N4 correct")
    antsimagesN4 = [ ants.n4_bias_field_correction(im, verbose=False) for im in images_reloaded]
    antsimagesN4 = [ im/im.mean() for im in antsimagesN4]

    mrrTemplate(antsimagesN4)

# only if we want a duplication of CISO
#mkCISOMRR()

def mkHalfwayMRR(tmplte, images, oname, blend=0):
    # MRR, but using the halfway space as a template
    templateNL = local_ants.build_template(initial_template = tmplte, image_list=images,
                                    iterations=4, aff_iterations=(1000,500,250,0),
                                    aff_shrink_factors = (6,4,2,1),
                                    aff_smoothing_sigmas=(4,2,1,0),
                                    reg_iterations=[100,100,70,20],
                                    grad_step = 0.1,
                                    flow_sigma=3,
                                    total_sigma=0,
                                    useNoRigid=True,
                                    blending_weight = blend,
                                    type_of_transform='SyN',
                                    verbose=True,aff_metric = METRIC, syn_metric = METRIC, random_seed = RS)
    fnameNL = os.path.join(workdir, oname)
    ants.image_write(templateNL,  fnameNL)


# MRR style template based using the halfway space as a starting point - should use Accumulated for that
# this uses the axial, which is what flywheel MRR does
print("Creating standard MRR")
mkHalfwayMRR(HighRes, antsimages, "nl_mrr_axial_sharpened.nii.gz")
print("Creating unsharpened MRR")
mkHalfwayMRR(HighRes, antsimages, "nl_mrr_axial_notsharpened.nii.gz", blend=1)
print("Done MRR")
# don't think there is any point doing more than one iteration
# We need to account for template drift if we do
# Use build_template with 0 iterations on the nonlinear bit to get
# our affine reg with multiple iterations
# for iterations in range(4):
#    # Now reregister everything to this, with potentially more complex transforms
#    regresults2 = [ants.registration(fixed=Accumulated, moving=moving, type_of_transform="Affine",
#                                     outprefix=mkAntsTmp())
#                   for moving in antsimages ]
#
#    warped = [x["warpedmovout"] for x in regresults2]
#    Accumulated = functools.reduce(lambda x,y: x+y, warped)
#    Accumulated /= float(len(warped))
#    outname = "avstep" + str(iterations + 1) + ".nii.gz"
#    ants.image_write(Accumulated,  os.path.join(workdir, outname))

# get transforms to the final template
regresults2 = [ants.registration(fixed=templateNL, moving=moving, type_of_transform="SyN",
                                 outprefix=mkAntsTmp(), aff_metric = METRIC, syn_metric = METRIC, random_seed = RS)
               for moving in antsimages ]

# these can be used to unwarp the original.
def unwarp(im, reg):
    # apply the affine, thn nonlinear transform, then inverted affine
    print(reg['fwdtransforms'])
    rr = [reg['fwdtransforms'][-1]] + reg['fwdtransforms']
    print(rr)
    invert = (1,0,0)
    return ants.apply_transforms(fixed = im, moving = im, transformlist = rr, whichtoinvert = invert)

#unwarpres = [unwarp(param[0], param[1]) for param in zip(antsimages, regresults2)]
#print(regresults2)
# save for inspection
#[ ants.image_write(param[0]['warpedmovout'], param[1]) for param in zip(regresults2, debugants )]
#[ ants.image_write(param[0], param[1]) for param in zip(unwarpres, uwnames )]

#print(uwnames)

quit()

# Now we can begin super resolution
# First step for testing - convert the rigid transforms to versor3d for mial to work
# later we'll change mial to support affine
# No longer required now that we're using NiftyMIC

import transforms3d

def sitkAff2Vers(tf):
    vtf=sitk.VersorRigid3DTransform()
    vtf.SetCenter(tf.GetCenter())
    vtf.SetTranslation(tf.GetTranslation())
    # Directly set the matrix? Original ought to be rotation only, if we requested rigid
    # decompose
    mat = sP.sitkTFtoNP(tf)
    decomp = transforms3d.affines.decompose44(mat)
    M = decomp[1]
    mat = M.T.ravel().tolist()
    vtf.SetMatrix(mat)
    return vtf

sitkTFs = [ sitk.ReadTransform(f["fwdtransforms"][0]) for f in regresults2 ]

if TFTYPE == "Affine":
    affnames = [ "tf_" + str(x+1) + ".mat" for x in range(len(sitkTFs)) ]
    affnames = [ os.path.join(workdir, v) for v in affnames ]
    [ sitk.WriteTransform(p[0], p[1]) for p in zip(sitkTFs, affnames) ]
else:
    sitkVersor3D = [ sitkAff2Vers(x) for x in sitkTFs ]
    versornames = [ "tf_" + str(x+1) + ".mat" for x in range(len(sitkVersor3D)) ]
    versornames = [ os.path.join(workdir, v) for v in versornames ]
    [ sitk.WriteTransform(p[0], p[1]) for p in zip(sitkVersor3D, versornames) ]
