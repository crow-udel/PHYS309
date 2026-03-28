#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
#

""" 
TIP: if you edit this script in hipe:
- all lines of code actually executed by this script will appear in black,
- all comments and explanations will appear in green,
- all line of codes that are not executed by default but that could be useful
  in some particular cases, will appear in red.

If you want to edit the script and remove comments and commented lines of code,
be careful in preserving the indentation. Otherwise if statements and loops
might not work.

This script invokes the complete Scanamorphos port for the PACS photometer 
pipeline. 

Description:
- This script takes your data from Level 1 to 2.5, starting from a Level 1 HSA 
  product.
- This script is written so that you can run it entirely in one go or 
  line-by-line.  
- Comments are included in this script, but for a full explanation see the PACS 
  Data Reduction Guide (PDRG): 
  Help Menu -> Help Contents -> PACS Data Reduction Guide: Photometry

Inputs:
- camera: Photometer camera ("red" or "blue"). Be careful with memory 
  allocation, blue requires about 4 times the memory used in red.
- obsids: List of observation ids with the scans and cross scans needed to 
  combine. The cross scan direction should be perpendicular to the scan 
  direction.
- solarSystemObject: Set it to True if the source is a Solar System object 
  (i.e. is a moving target). Note that if there is extended emission in the
  map, the main JScanam assumptions will fail, because the scan and cross 
  scan will be shifted to the object reference system.
- galactic: The galactic option should be set to True if the maps extended 
  emission is not confined to a small area. The safest is to set it always 
  True. 
- calculateRaDec: Calculate the coordinates of every pixel with time. The 
  tasks will run faster if this option is set to True, but they will require 
  a factor of 2-3 more memory.
- showMapsAfterTasks: If True, it will produce a map after each step in the 
  pipeline.
- debug: Debug mode. Use wisely, debug=True will clutter your desktop with 
  graphs and displays. It will also consume more memory.
- deglitch: Set it to True if you want to deglitch the frames with a lower
  threshold than the one used by default (nSigmaDeglitch = 5).
- nSigmaDeglitch: The new deglitch threshold to use in case the deglitch 
  parameter is set to True.
- makeFinalMap: Generate a final map.
- outputPixelsize: Pixel size to use for the final map. It will use 1.6" for 
  the blue camera and 3.2" for the red camera if it's set to -1. 
- pixfrac: This parameter is used for the final map. It fixes the ratio of the
  dropsize (to be used in the Drizzle projection algorithm) to the input detector 
  pixel size.

########################### SECTION 0 ##########################################
################################################################################
########################### SETTINGS ###########################################
"""
from java.lang import System
System.setProperty("https.protocols","TLSv1.2");
"""Note: these are 100 micron scans"""
camera = "blue"
obsid1 = 1342213171
obsid2 = 1342213172
solarSystemObject = False
galactic = True
calculateRaDec = False
showMapsAfterTasks = False
debug = False
deglitch = True
nSigmaDeglitch = 3
makeFinalMap = True
outputPixelSize = -1
pixfrac = 0.1

"""
########################### SECTION 1 ##########################################
################################################################################
########################### PROCESSING #########################################
"""
print " -- Executing JScanam -- " 

"""
 Load Level 1 data from the HSA
"""
print " Loading Level 1 data "
scansObs = getObservation(obsid1, useHsa=True, instrument="PACS")
level1 = PacsContext(scansObs.level1)
scans = level1.averaged.getCamera(camera).product.selectAll() 
blueFilter1 = scansObs.meta["blue"].value

cscansObs = getObservation(obsid2, useHsa=True, instrument="PACS")
level1 = PacsContext(cscansObs.level1)
cscans = level1.averaged.getCamera(camera).product.selectAll()
blueFilter2 = cscansObs.meta["blue"].value

calTree = getCalTree(obs=cscansObs)

if (camera == "blue" and blueFilter1 != blueFilter2):
  print ""
  print "  ------------------------------------------------------------------"
  print "  ALERT!!! You are trying to combine two blue observations obtained "  
  print "  with different filters, which is most probably wrong!             "
  print "  ------------------------------------------------------------------"
  print ""

"""
 Set the scans and cross scans coordinates to the object reference system if it's a SSO
"""
if(solarSystemObject):
  print " Setting the scans and crossScans coordinates to the object reference system "
  pp = scansObs.auxiliary.pointing
  orbitEphem = scansObs.auxiliary.orbitEphemeris
  horizons = scansObs.auxiliary.horizons
  cal = getCalTree(obs=scansObs)
  scans = photAddInstantPointing(scans, pp, calTree=cal, orbitEphem=orbitEphem, horizonsProduct=horizons)
  timeOffset = scans.getStatus("FINETIME")[0]
  scans = correctRaDec4Sso(scans, timeOffset=timeOffset, orbitEphem=orbitEphem, horizonsProduct=horizons, linear=0)
  #
  pp = cscansObs.auxiliary.pointing
  orbitEphem = cscansObs.auxiliary.orbitEphemeris
  horizons = cscansObs.auxiliary.horizons
  cal = getCalTree(obs=cscansObs)
  cscans = photAddInstantPointing(cscans, pp, calTree=cal, orbitEphem=orbitEphem, horizonsProduct=horizons)
  cscans = correctRaDec4Sso(cscans, timeOffset=timeOffset, orbitEphem=orbitEphem, horizonsProduct=horizons, linear=0)
  del pp, orbitEphem, horizons, cal, timeOffset

"""
 Calculate the coordinates of every pixel with time
"""
if(calculateRaDec):
   print " Calculating the coordinates of every pixel with time "
   scans = photAssignRaDec(scans, calTree=calTree)
   cscans = photAssignRaDec(cscans, calTree=calTree)

"""
 Reduce the product size removing unnecessary status information
"""
# print " Reducing the product size "
# scans = scanamorphosReduceFramesSize(scans)
# cscans = scanamorphosReduceFramesSize(cscans)

"""
 Remove turn arounds
"""
print " Removing turn arounds "
scans = scanamorphosRemoveTurnarounds(scans, limit=50.0, debug=debug)
cscans = scanamorphosRemoveTurnarounds(cscans, limit=50.0, debug=debug)

"""
 Mask long term glitches. 
 This task produces a mask called Scanam_LongTermGlitchMask. You should check 
 this mask and if the results are not optimal (some glitches are not detected 
 or some sources are flagged as glitches), you can try to modify the parameter 
 stepAfter in order to get better results.
"""
print " Masking long term glitches "
scans = scanamorphosMaskLongTermGlitches(scans, stepAfter=20, galactic=galactic, calTree=calTree, debug=debug)
cscans = scanamorphosMaskLongTermGlitches(cscans, stepAfter=20, galactic=galactic, calTree=calTree, debug=debug)

if(showMapsAfterTasks):
   map, mi = photProject(scans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Scans after masking long term glitches")
   map, mi = photProject(cscans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Cross scans after masking long term glitches")

"""
 Save the scans and cross scans for later use
"""
print " Saving a copy of the scans and cross scans in a temporal pool "
from herschel.pacs.share.util import PacsProductSinkWrapper
scansRef = PacsProductSinkWrapper.getInstance().saveAlways(scans)
cscansRef = PacsProductSinkWrapper.getInstance().saveAlways(cscans)

"""
 Subtract the baseline of every scanleg in every pixel
"""
print " Initial baseline subtraction / per scanleg "
scans = scanamorphosScanlegBaselineFitPerPixel(scans, nSigma=2)
cscans = scanamorphosScanlegBaselineFitPerPixel(cscans, nSigma=2)

if(showMapsAfterTasks):
   map, mi = photProject(scans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Scans after baseline subtraction / per scanleg")
   map, mi = photProject(cscans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Cross scans after baseline subtraction / per scanleg")

"""
 Create the source mask.
 Modify nSigma until you get an optimal mask. The mask should cover only a small
 fraction of the map (<~30%). It's not necessary that all the faint emission is
 masked, only the brightest regions.
"""
print " Creating a source mask "
scans.join(cscans)
del(cscans)
sourceImage, scans = scanamorphosCreateSourceMask(scans, nSigma=4.0, createMask=False, galactic=galactic, calTree=calTree, debug=debug)

if(showMapsAfterTasks):
   d = Display(sourceImage, title="Masked sources")

"""
 Replace the scans and cross scans by the saved ones
"""
print " Loading the saved scans and cross scans from the temporal pool "
del(scans)
scans = scansRef.product
cscans = cscansRef.product
del(scansRef, cscansRef)
System.gc()

"""
 Add the mask to the saved scans and cross scans
"""
print " Adding the source mask to the scans and cross scans "
maskImage, scans = scanamorphosCreateSourceMask(scans, inputImage=sourceImage, createMask=True, calTree=calTree, debug=debug)
maskImage, cscans = scanamorphosCreateSourceMask(cscans, inputImage=sourceImage, createMask=True, calTree=calTree, debug=debug)

"""
 Baseline subtraction. 
 Here we use galactic=True, because we only want to remove an offset
"""
print " Baseline subtraction "
scans = scanamorphosBaselineSubtraction(scans, galactic=True, debug=debug)
cscans = scanamorphosBaselineSubtraction(cscans, galactic=True, debug=debug)

if(showMapsAfterTasks):
   map, mi = photProject(scans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Scans after baseline subtraction")
   map, mi = photProject(cscans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Cross scans after baseline subtraction")

"""
 Baseline pre-processing
"""
print " Baseline pre-processing "
scans = scanamorphosBaselinePreprocessing(scans, debug=debug)
cscans = scanamorphosBaselinePreprocessing(cscans, debug=debug)

if(showMapsAfterTasks):
   map, mi = photProject(scans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Scans after baseline pre-processing")
   map, mi = photProject(cscans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Cross scans after baseline pre-processing")

"""
 Baseline subtraction
"""
print " Baseline subtraction "
scans = scanamorphosBaselineSubtraction(scans, galactic=galactic, debug=debug)
cscans = scanamorphosBaselineSubtraction(cscans, galactic=galactic, debug=debug)

if(showMapsAfterTasks):
   map, mi = photProject(scans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Scans after baseline subtraction")
   map, mi = photProject(cscans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Cross scans after baseline subtraction")

"""
 Destriping
"""
print " Destriping "
scans, cscans = scanamorphosDestriping(scans, cscans, iterations=6, calTree=calTree, debug=debug)

if(showMapsAfterTasks):
   map, mi = photProject(scans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Scans after destriping")
   map, mi = photProject(cscans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Cross scans after destriping")

"""
 Merge the scans and cross scans
"""
print " Merging scans and cross scans in a frame product "
System.gc()
mergedScans = scans
mergedScans.join(cscans)
del(scans, cscans)
System.gc()

"""
 Deglitch the merged scans
 
 We use nSigma=5, and not a lower value, to avoid masking signal with strong drifts 
 that will be corrected by the scanamorphosIndividualDrifts task
"""
print " Deglitching the merged scans "
scanamorphosDeglitch(mergedScans, nSigma=5, calTree=calTree, debug=debug)

"""
 Remove individual drifts
"""
print " Removing individual drifts "
mergedScans = scanamorphosIndividualDrifts(mergedScans, calTree=calTree, debug=debug)

if(showMapsAfterTasks):
   map, mi = photProject(mergedScans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(map, title="Merged scans after individual drifts correction")

"""
 Deglitch the merged scans again with the user provided nSigmaDeglitch

 Be careful setting nSigmaDeglitch to very small values, because it can
 mask bright sources. You should always check the mask called Scanamorphos_GlitchMask
 to make sure that the task is masking only real glitches.
"""
if(deglitch):
    print " Deglitching the merged scans "
    scanamorphosDeglitch(mergedScans, nSigma=nSigmaDeglitch, calTree=calTree, debug=debug)

"""
 Create a new source mask, and remove the background median offset 
"""
# print " Obtaining a new source mask "
# sourceImage, mergedScans = scanamorphosCreateSourceMask(mergedScans, nSigma=2.0, createMask=True, galactic=galactic, calTree=calTree, debug=debug)
#
# if(showMapsAfterTasks):
#   d = Display(sourceImage, title="Masked sources")
#
# print " Removing the background median offset "
# mergedScans = scanamorphosSubtractOffset(mergedScans, calTree=calTree)
#
# if(showMapsAfterTasks):
#   map, mi = photProject(mergedScans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
#   d = Display(map, title="Merged scans after the background median offset subtraction")


"""
########################### SECTION 2 ##########################################
################################################################################
########################### SEE THE RESULTS#####################################

 Save the frames - NOT HELPFUL
"""
# outputFitsFile = "/Users/sdr/Dropbox/HST/HD166/jscanam_100micron.fits"
# FitsArchive().save(outputFitsFile, mergedScans)

"""
 Final projection with drizzle
"""
if(outputPixelSize < 0):
   if(camera == "blue"):
      from herschel.pacs.spg.phot import PhotHelper
      if(PhotHelper.isParallelObs(mergedScans)):
         outputPixelSize = 3.2
      else:
         outputPixelSize = 1.6
   else:
      outputPixelSize = 3.2

if(makeFinalMap):
   print " Projecting the merged scans onto the final map "
   print " Projection with drizzle using: " + str(outputPixelSize) + " arcsec"
   finalMap, mi = photProject(mergedScans, outputPixelsize=outputPixelSize, pixfrac=pixfrac, calTree=calTree, useMasterMaskForWcs=False)
   d = Display(finalMap, title="Final map")
   simpleFitsWriter(product=finalMap, file='58Eri_jscanam_100micron.fits')
   # outputMapFile = "/your/Home/Dir/finalMap_"+camera+".jpg"
   # d.saveAsJPG(outputMapFile)
