#!/usr/bin/env python

# Use simple ITK to implement the preprocessing steps required of mialsrtk, but
# on volumes not slices

import argparse
import sys

def is_readable_file(parser, arg):
    try:
        f = open(arg, 'r')
        f.close()
    except Exception:
        raise argparse.ArgumentTypeError("{0} does not exist or is not readable".format(arg))

    return(arg)


#parser = argparse.ArgumentParser(description="Preprocessing of images for superresolution")
#
#parser.add_argument("-i", "--images", nargs="+",
#                    type=lambda x: is_readable_file(parser, x),
#                    required=True,
#                    help="One or more input nifti files")
#
#parser.add_argument("-m", "--masks", nargs="+",
#                    type=lambda x: is_readable_file(parser, x),
#                    required=True,
#                    help="One or more input nifti mask files (same length as iamges")
#
#parser.add_argument("-o", "--out", nargs="+",
#                    required=True,
#                    help="One or more output nifti files (same length as images")
#
#
#args = parser.parse_args()
#
#if len(args.images) != len(args.masks):
#    print("Masks and images must match")
#    exit()
#print(args.out)
#if len(args.images) != len(args.out):
#    print("Out and images must match")
#    exit()


################
##
import SimpleITK as sitk
import re

# helper functions
def N4(image, mask):
    # just in case we need to get complicated and not use defaults
    return sitk.N4BiasFieldCorrection(image, mask)
##


def IntStd(images, masks, targetMax = 255):
    ## standardize intensity across images
    ## This one retains relative intensity
    def MaxMask(i, m):
        statsfilt = sitk.LabelIntensityStatisticsImageFilter()
        statsfilt.SetBackgroundValue(0)
        statsfilt.Execute(m, i)
        return statsfilt.GetMaximum(1)
    

    def ScaleOne(i, thismax, globalmax):
        newmax = (thismax/globalmax)*targetMax
        rescaler = sitk.RescaleIntensityImageFilter()
        rescaler.SetOutputMinimum(0)
        rescaler.SetOutputMaximum(newmax)
        res = rescaler.Execute(i)
        return res

    maxvals = [ MaxMask(*pair) for pair in zip(images, masks) ]
    glbmax = max(maxvals)
    imscaled = [ ScaleOne(params[0], params[1], glbmax) for params in zip(images, maxvals) ]
    return imscaled

def IntStd2(images, masks, targetMean=100, targetSD=50):
    ## standardize intensity across images
    ## This one sets a common median
    ## removes negative values
    def MeanSDMask(i, m):
        statsfilt = sitk.LabelIntensityStatisticsImageFilter()
        statsfilt.SetBackgroundValue(0)
        statsfilt.Execute(m, i)
        return {'Mean' : statsfilt.GetMean(1), "SD" : statsfilt.GetStandardDeviation(1) }
    

    def ScaleOne(i, dist):
        MN = dist["Mean"]
        SD = dist["SD"]

        res = targetSD * ((i - MN)/SD) + targetMean
        # remove negative values
        res = res * sitk.Cast(res > 0, res.GetPixelID())
        return res

    maxvals = [ MeanSDMask(*pair) for pair in zip(images, masks) ]
    print(maxvals)
    imscaled = [ ScaleOne(params[0], params[1]) for params in zip(images, maxvals) ]
    return imscaled


###########################################################
# Histogram normalization stuff from mialsrtkHistogramNormalization.py
import numpy as np
import scipy.ndimage
import copy

####################
def thickSliceLast(im):
    # permute axes so the thick slice is in the last dimension
    sp = np.array(im.GetSpacing())
    print(sp)
    thickest = sp.argmax()
    print(thickest)
    if thickest == 0:
        neworder=(1,2,0)
    elif thickest == 1:
        neworder=(0,2,1)
    else:
        return(im)

    return sitk.PermuteAxes(im, neworder)

def percentile_nonzero(image, percentile_nonzero):
    # pdb.set_trace()
    if len(image) < 1:
        value = None
    elif percentile_nonzero >= 100:
        sys.stderr.write(
            "ERROR: percentile_nonzero must be < 100.  you supplied: %s\n"
            % percentile_nonzero
        )
        value = None
    else:
        image_nonzero = image[image != 0]
        element_idx = int(len(image_nonzero) * (percentile_nonzero / 100.0))
        image_nonzero.sort()
        value = image_nonzero[element_idx]
    return value


def mean_nonzero(image):
    image_nonzero = image[image != 0]
    mean = np.sum(image_nonzero) / len(image_nonzero)
    return mean


def intensityNormalization(image, landmarks):
    print("min =" + str(landmarks["p1"]))
    print("max (99.8%) =", str(landmarks["p2"]))
    # print 'mean ='+str(landmarks['mean'])
    print("quartiles [25%,50%,75%] =" + str(landmarks["quartiles"]))
    return 1


def displayHistogram(image, image_name, loffset, roffset):
    bins = np.round(np.arange(loffset, np.max(image) - roffset, 40))
    histo, bins = np.histogram(image, bins=bins)
    bins_center = 0.5 * (bins[1:] + bins[:-1])
    # pyplot.plot(bins_center,histo,alpha=0.5)
    ##pyplot.hist(fit(np.random.uniform(x[0],x[-1],len(image))),bins=y)
    ##pyplot.hist(image,bins,histtype='step',alpha=0.5,label=image_name)
    return 1


def extractImageLandmarks(image):
    landmarks = {}
    landmarks["p1"] = percentile_nonzero(image, 0)
    landmarks["p2"] = percentile_nonzero(image, 99.8)
    # landmarks['mean']=mean_nonzero(image)
    landmarks["quartiles"] = [
        percentile_nonzero(image, 25),
        percentile_nonzero(image, 50),
        percentile_nonzero(image, 75),
    ]
    # landmarks['quartiles']=[percentile_nonzero(image,10),percentile_nonzero(image,20),percentile_nonzero(image,30),percentile_nonzero(image,40),percentile_nonzero(image,50),percentile_nonzero(image,60),percentile_nonzero(image,70),percentile_nonzero(image,80),percentile_nonzero(image,90)]
    # pdb.set_trace()
    return landmarks


def trainImageLandmarks(list_landmarks, verbose=False):
    mup_l = []
    mup_L = []
    mup_r = []
    mup_R = []
    maxLR = []
    index = 0
    while index < len(list_landmarks):
        landmarks = list_landmarks[index]["quartiles"]
        mup_l.append(np.min(landmarks - list_landmarks[index]["p1"]))
        mup_L.append(np.max(landmarks - list_landmarks[index]["p1"]))
        mup_r.append(np.min(list_landmarks[index]["p2"] - landmarks))
        mup_R.append(np.max(list_landmarks[index]["p2"] - landmarks))
        maxLR.append(
            np.max(
                [
                    float(mup_L[index]) / mup_l[index],
                    float(mup_R[index]) / mup_r[index],
                ]
            )
        )
        # print 'mup_l  =  '+str(mup_l[index])
        # print 'mup_L  =  '+str(mup_L[index])
        # print 'mup_r  =  '+str(mup_r[index])
        # print 'mup_R  =  '+str(mup_R[index])
        # print 'maxLR  =  '+str(maxLR[index])
        index += 1
    ymax = np.max(maxLR)
    ymax_index = maxLR.index(max(maxLR))
    dS = float(ymax * (mup_L[ymax_index] + mup_R[ymax_index]))
    if verbose:
        print(
            "Ymax  =  ",
            str(ymax)
            + "  at position "
            + str(ymax_index)
            + "  ,  dS = "
            + str(dS)
            + " (=s2 when s1=0)",
        )
    return list_landmarks, dS


def mapImageLandmarks(list_landmarks, s1, s2, verbose=False):
    list_landmarks_mapped = copy.deepcopy(list_landmarks)
    index = 0
    while index < len(list_landmarks):
        land_index = 0
        if verbose:
            print("Image index:", str(index))
        while land_index < len(list_landmarks[index]["quartiles"]):
            if verbose:
                print(
                    "old landmark:",
                    str(list_landmarks_mapped[index]["quartiles"][land_index]),
                )
            list_landmarks_mapped[index]["quartiles"][land_index] = s1 + float(
                (
                    list_landmarks_mapped[index]["quartiles"][land_index]
                    - list_landmarks_mapped[index]["p1"]
                )
                / float(
                    list_landmarks_mapped[index]["p2"]
                    - list_landmarks_mapped[index]["p1"]
                )
            ) * float((s2 - s1))
            if verbose:
                print(
                    "new landmark:",
                    str(list_landmarks_mapped[index]["quartiles"][land_index]),
                )
            land_index += 1
        if verbose:
            print(
                "p1, p2 = ",
                str(list_landmarks_mapped[index]["p1"]),
                ",",
                str(list_landmarks_mapped[index]["p2"]),
            )
        index += 1
    return list_landmarks_mapped


def verifyOne2OneMapping(s1, s2, list_landmarks, lmap_mean):
    landmarks = list_landmarks["quartiles"]
    mup_L = np.max(landmarks - list_landmarks["p1"])
    mup_R = np.max(list_landmarks["p2"] - landmarks)
    # print 'mup_L  =  '+str(mup_L)
    # print 'mup_R  =  '+str(mup_R)

    land_index = 0
    while land_index < len(lmap_mean):
        # pdb.set_trace()
        if np.logical_and(
            (lmap_mean[str(0)] - s1) >= mup_L,
            (s2 - lmap_mean[str(len(lmap_mean) - 1)]) >= mup_R,
        ):
            cond = 1
        else:
            cond = 0
        land_index += 1
    return cond


def mapImage(image, lmap_mean, list_landmarks, s1, s2, p1, p2):
    image_out = image.copy().astype("float")
    tmp = image.copy().astype("float")
    index = 0
    ##pyplot.figure(2)
    # pdb.set_trace()
    while index < len(lmap_mean) + 1:
        if index == 0:
            x = np.array([int(p1), int(list_landmarks[index])])
            y = np.array([int(s1), int(lmap_mean[str(index)])])
            coefs = np.polyfit(x, y, 1)
            mask = np.logical_and(image > 0, image <= list_landmarks[index])
            image_out[mask] = coefs[0] * image[mask] + coefs[1]
            ##pyplot.plot(x,y,marker='o', linestyle='--');
        else:
            if index == (len(lmap_mean)):
                x = np.array([int(list_landmarks[index - 1]), int(p2)])
                y = np.array([int(lmap_mean[str(index - 1)]), int(s2)])
                coefs = np.polyfit(x, y, 1)
                mask = image > list_landmarks[index - 1]
                image_out[mask] = coefs[0] * image[mask] + coefs[1]
                ##pyplot.plot(x,y,marker='o', linestyle='--');
            else:
                x = np.array(
                    [
                        int(list_landmarks[index - 1]),
                        int(list_landmarks[index]),
                    ]
                )
                y = np.array(
                    [
                        int(lmap_mean[str(index - 1)]),
                        int(lmap_mean[str(index)]),
                    ]
                )
                coefs = np.polyfit(x, y, 1)
                mask = np.logical_and(
                    image > list_landmarks[index - 1],
                    image <= list_landmarks[index],
                )
                image_out[mask] = coefs[0] * image[mask] + coefs[1]
                ##pyplot.plot(x,y,marker='o', linestyle='--');
        index += 1
    ##pyplot.show()
    return image_out


def computeMeanMapImageLandmarks(list_landmarks, verbose = False):
    mean_landmarks = {}
    index = 0
    while index < len(list_landmarks):
        land_index = 0
        while land_index < len(list_landmarks[index]["quartiles"]):
            if index == 0:
                mean_landmarks[str(land_index)] = list_landmarks[index][
                    "quartiles"
                ][land_index]
            else:
                mean_landmarks[str(land_index)] += list_landmarks[index][
                    "quartiles"
                ][land_index]
            land_index += 1
        index += 1

    land_index = 0
    while land_index < len(mean_landmarks):
        mean_landmarks[str(land_index)] = mean_landmarks[
            str(land_index)
        ] / len(list_landmarks)
        land_index += 1
    if verbose:
        print("Final landmark average : ")
        print(mean_landmarks)
    return mean_landmarks

def HistNorm(images, masks, verbose = False):
    list_landmarks = []

    s1 = 1
    # pyplot.figure(1)
    # pyplot.subplot(211)
    for index in range(len(images)):
        image = sitk.GetArrayFromImage(images[index])
        # image = scipy.ndimage.filters.gaussian_filter(image,1.0)
        mask = sitk.GetArrayFromImage(masks[index])
        maskedImage = np.reshape(
            image * mask, image.shape[0] * image.shape[1] * image.shape[2]
        )
        list_landmarks.append(extractImageLandmarks(maskedImage))

    list_landmarks, dS = trainImageLandmarks(list_landmarks, verbose)
    s2 = np.ceil(dS - s1)
    list_landmarks_mapped = mapImageLandmarks(list_landmarks, s1, s2, verbose)

    mean_landmarks = computeMeanMapImageLandmarks(
        list_landmarks_mapped, verbose
    )

    result = list()
    for index in range(len(images)):
        image_data = sitk.GetArrayFromImage(images[index])
        mask_data = sitk.GetArrayFromImage(masks[index])
        dimY = image_data.shape[0]
        dimX = image_data.shape[1]
        dimZ = image_data.shape[2]

        maskedImage = np.reshape(image_data, dimX * dimY * dimZ)
        # pdb.set_trace()
        maskedImageMapped = mapImage(
            maskedImage,
            mean_landmarks,
            list_landmarks[index]["quartiles"],
            s1,
            s2,
            list_landmarks[index]["p1"],
            list_landmarks[index]["p2"],
        )

        o2o = verifyOne2OneMapping(
            s1, s2, list_landmarks[index], mean_landmarks
        )
        new_image = sitk.GetImageFromArray(np.reshape(maskedImageMapped,
                                                      np.array([dimY, dimX, dimZ])))
        
        new_image.CopyInformation(images[index])
        result.append(new_image)
        #sitk.WriteImage(new_image, outputs[index])
    return result


    
import transforms3d

def mkHalfway(A):    
    # This uses a different convention to fsl (order of
    # composition)
    # fsl is rotation translation scale skew
    # this one is translation rotation scale skew
    # rotation, scale and translation components come
    # out the same as avscale, but skew is very different
    # could be to do with 3x3 vs 4x4
    decomp = transforms3d.affines.decompose44(A)
    #return (decomp, A)
    qt = transforms3d.quaternions.mat2quat(decomp[1])
    # We want a rotation that is half the original.
    # which corresponds to a quaternion square root
    qth = transforms3d.quaternions.qpow(qt, 0.5)
    newzoom = np.sqrt(decomp[2])
    newzoommat = np.matrix(np.eye(4))
    newskewmat = newzoommat.copy()
    newrot = newzoommat.copy()
    newrot[0:3, 0:3] = transforms3d.quaternions.quat2mat(qth)
    
    newzoommat[0:3, 0:3] = np.diag(newzoom)
    
    newskew = 0.5 * decomp[3]
    newskewmat[0:3,0:3] = transforms3d.affines.striu2mat(newskew)

    # translation
    X = A @ newskewmat.I @ newzoommat.I @ newrot.I + np.matrix(np.eye(4))
    Xsub = X[0:3, 0:3]

    trans =  decomp[0]
    # turn these into matrices
    newtrans = Xsub.I @ trans

    result = newrot @ newzoommat @ newskewmat
    result[0:3,3] = newtrans.T
    return result
    #return (newrot, newzoommat, newskewmat, decomp, A)

def sitkTFtoNP(tf):
    #A = np.asarray(tf.GetParameters()).reshape(4,3).T
    #A = np.append(A, np.array([0,0,0,1]).reshape(1,4), 0)
    M = np.array(tf.GetMatrix()).reshape((3,3))
    tr = np.array(tf.GetTranslation())
    cent = np.array(tf.GetCenter())
    adjust_trans = tr + cent - np.dot(M, cent)
    a44 = np.eye(4)
    a44[0:3,0:3] = M
    a44[0:3,3] = adjust_trans.flatten()
    return a44

def sitkHalfway(tf):
    A = sitkTFtoNP(tf)
    return mkHalfway(A)

def NPtoSITK(affmatrix44, cent=np.array((0.0,0.0,0.0))):
    cent = cent.reshape(3,1)
    M = affmatrix44[0:3,0:3]
    tr = affmatrix44[:3, 3].reshape(3,1)
    adjusted = tr - cent + np.dot(M, cent)
    tf = sitk.AffineTransform(3)
    tf.SetMatrix(np.array(M).flatten())
    tf.SetTranslation(np.array(adjusted).flatten())
    tf.SetCenter(np.array(cent).flatten())
    return tf

def sitkWriteAffine(affmatrix44, fname, cent=np.array((0.0,0.0,0.0))):
    aff = NPtoSITK(affmatrix44, cent)
    sitk.WriteTransform(aff, fname)

# fsl does some fancy stuff to produce the square root.
# rotation is via quaternions
# translation is less obvious.

# th = (mat*scale.i()*skew.i()*rot.i() + id4).SubMatrix(1,3,1,3).i() 
#	* trans.SubMatrix(1,3,1,1);
# The stuff in brackets (mat * ) produces a scale component, that we chop off
# this means we're only in the halfway space for non transform components
# Adding identity before inverting means we "double" the halfway transform,
# and we therefore get something smaller when we invert. Composing that
# with the original translation gives us the answer
# The remaining 3x3 is the forward half transform, which is then inverted again
##########################################################
#n4names = [ re.sub("\\.nii\\.gz", "_N4.nii.gz", s) for s in args.images ]
#print(n4names)

#masks = [ sitk.ReadImage(x, sitk.sitkUInt8) for x in args.masks ]
#images =  [ sitk.ReadImage(x, sitk.sitkFloat32) for x in args.images ]

#imagesN4 = [ N4(images[idx], masks[idx]) for idx in range(len(images)) ]

#imagesN4 = [ N4(*pair) for pair in zip(images, masks) ]

#stdimages = IntStd(imagesN4, masks)

#HistNorm(stdimages, masks, args.out)

