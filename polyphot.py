#!/usr/bin/env python
from __future__ import print_function
#=============================================================================#
#                                                                             #
# NAME:     polyphot.py                                                       #
#                                                                             #
# PURPOSE:  Script to to load a FITS image and allow the user to measure the  #
#           flux of a source by drawing a polygon around the emission.        #
#                                                                             #
# MODIFIED: 30-Apr-2018 by C. Purcell                                         #
#                                                                             #
#=============================================================================#

# Fractional amplitude calibration uncertainty (CORNISH = 0.089)
errCalib_frac = 0.0

# Absolute error in the position
errPosition_asec = 0.1

# Clean bias for radio data, e.g. for CORNISH CB = -0.94 mJy/beam
cleanBias_sigma = 0.0

# Scaling factor to convert the data units to Jy/beam
scaleUnitsJy = 8.461595e-3   # IRAC Band 4 MJy/SR to Jy
scaleUnitsJy = 1.0

# Directory in which older versions of the save files are stored. The script
#  will attempt to load the polygons from these if they exist.
oldDatDir = '.'

# Size of the figure (x, y) in inches
figSize_inch = (10.0, 8.0)

# Plot positions in a kvis annotation file in this directory (if it exists)?
doPlotAnn=True

#=============================================================================#

import os
import sys
import copy
import re
import math as m
from matplotlib.collections import RegularPolyCollection
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.widgets import Button, Slider
from matplotlib.pyplot import figure, show
import matplotlib.cm as cm
from pylab import *
import astropy.io.fits as pf
import astropy.wcs.wcs as pw
import numpy as np
from scipy import stats
from scipy.stats import norm
from scipy import median
from scipy.stats.stats import nanmedian, _nanmedian

#from util_misc import 
from util_fits import strip_fits_dims
from util_fits import mkWCSDict
from util_fits import get_beam_from_history
from util_fits import calc_beam_area
#from util_fits import *
from util_phot import polylst_parse
from util_phot import get_pix_in_poly
from util_phot import measure_pix_props
from util_phot import poly_area


#-----------------------------------------------------------------------------#
def main():

    oldSrcPoly = []
    oldSkyPoly = []
    oldExcPoly = []
    
    # Get the filename from the command line
    args = sys.argv[1:]
    if len(args)<1: usage()
    else: filename = args[0]
    
    # Read the FITS header and image data from the file
    header=pf.getheader(filename)
    xydata = pf.getdata(filename)
    xydata, header = strip_fits_dims(xydata, header, 2, 5)
    w = mkWCSDict(header)
    wcs = pw.WCS(w['header2D'])
    
    # Search for beam keywords in the FITS file (for radio data)
    doDivideBeamArea = True
    if "BMAJ" in header and "BMIN" in header:
        header.update('BPA', 0.0)
        print("Found standard BMAJ and BMIN header keywords.")
    elif "CLEANBMAJ" in header and "CLEANBMIN" in header:
        print("Found AIPS style BMAJ and BMIN header keywords.")
        header.update('BMAJ', header['CLEANBMJ'])
        header.update('BMIN', header['CLEANBMN'])
        header.update('BPA', 0.0)
    else:
        bmaj, bmin, bpa = get_beam_from_history(header)
        if bmaj and bmin:
            print("Found AIPS beam header in history.")
            header.update('BMAJ', bmaj)
            header.update('BMIN', bmin)
            header.update('BPA', 0.0)
        else:
            doDivideBeamArea = False
            header.update('BMAJ', 0.0)
            header.update('BMIN', 0.0)
            header.update('BPA', 0.0)
            print("No beam header keywords located.")
            print("WARNING: will not divide by beam area.")
    
    # Determine the max and min of the data array
    zmin = np.nanmin(xydata)
    zmax = np.nanmax(xydata)

    # Load the KVIS annotation file if it exists
    fitsRoot, ext = os.path.splitext(filename)
    annfilename = fitsRoot + '.ann'
    if os.path.exists(annfilename) and doPlotAnn:
        ann_list = kvis_ann_parse(annfilename)
        ann_list = wcs.wcs_world2pix(ann_list, 0)
            
    # Load a pre-existing polygon file
    oldDat = None
    if os.path.exists(oldDatDir):
        if os.path.exists(oldDatDir +'/'+ fitsRoot + '.polyflux.dat'):
            oldDat = oldDatDir + '/'+ fitsRoot + '.polyflux.dat'
    if len(args)>=2:
        if os.path.exists(args[1]):
            oldDat = args[1]
    if oldDat:
        print("Found previous dat file '%s'." % oldDat)            
        oldMeasureLst = datfile_read(oldDat)
        if len(oldMeasureLst) == 1:
            mIndx = 0
        elif len(oldMeasureLst) > 1:
            print("Found %i individual measurements" % len(oldMeasureLst))
            mIndx = int(raw_input("Choose a measurement to load [1]:") or 1)-1
            if mIndx < 0 or mIndx > len(oldMeasureLst):
                print("Please choose between 1 and %s." % len(oldMeasureLst))
                sys.exit()                
        if "POLYSRC" in oldMeasureLst[mIndx]:
            print("Loading source polygon ...")
            oldSrcPoly =  oldMeasureLst[mIndx]['POLYSRC']
        if "POLYSKY" in oldMeasureLst[mIndx]:
            print("Loading sky polygon ...")
            oldSkyPoly =  oldMeasureLst[mIndx]['POLYSKY']

    # Define a figure on which to plot the image
    fig = figure(figsize=figSize_inch)
    axplot = fig.add_subplot(1,1,1)
    axplot.set_title("Click to outline the source emission.")
    subplots_adjust(bottom=0.12)

    # Plot the image data and colorbar
    cax = axplot.imshow(xydata,interpolation='nearest',origin='lower',
                        cmap=cm.jet,vmin=zmin,vmax=zmax)
    cbar=colorbar(cax,pad=0.0)

    # Plot the fits annotation
    if os.path.exists(annfilename) and doPlotAnn:
        
        # Create a circle for each Gaussian fit
        ann = RegularPolyCollection(
            10,
            rotation=0,
            sizes=(5,),
            facecolors = 'none',
            edgecolors = 'black',
            linewidths = (1,),
            offsets = ann_list,
            transOffset = axplot.transData
            )
        axplot.add_collection(ann)
        
    # Add the buttons to the bottom edge
    axreset = axes([0.49, 0.025, 0.1, 0.04])
    axmeasure = axes([0.60, 0.025, 0.1, 0.04])
    axsave = axes([0.71, 0.025, 0.1, 0.04])
    breset = Button(axreset, 'Reset')
    bmeasure = Button(axmeasure, 'Measure')
    bsave = Button(axsave, 'Save')
    
    # Polygon editor class 
    class ThreePolyEditor:
        """
        Left-click to add a point to the current polygon. Middle-click to delete
        the last point and right-click to close the polygon. Three polygons
        added in sequence: source-aperture, sky aperture and exclusion aperture.
        """
    
        def __init__(self):
            self.mode = 'src' 
        
            # Lists to store the vertices of the source, sky background and
            # exclusion-zone polygons
            self.offsets_src = []
            self.offsets_sky = []
            self.offsets_exc = []        
            self.offsets = []        # general working array
            self.apertures = []      # completed polygons for plotting
            self.measurements = {}
            
            # Working polygon collection (points plotted as small polys)
            self.points = RegularPolyCollection(
                10,
                rotation=0,
                sizes=(50,),
                facecolors = 'white',
                edgecolors = 'black',
                linewidths = (1,),
                offsets = self.offsets,
                transOffset = axplot.transData
                )

        # Handle mouse clicks
        def onclick(self,event):
            
            # Disable click events if using the zoom & pan modes.
            if fig.canvas.widgetlock.locked():
                return

            if event.button==1:
            
                # Left-click on plot = add a point
                if axplot.contains(event)[0] and self.mode != 'done':
                    x = event.xdata
                    y = event.ydata
                    self.offsets.append((x,y))
                    self.points.set_offsets(self.offsets)
                    if len(self.offsets)==1 : axplot.add_collection(self.points)
                    self.update()
                
                # Right-click on wedge = halve the lower colour clip
                if cbar.ax.contains(event)[0]:
                    clims = cax.get_clim()
                    cax.set_clim(clims[0]/2.0,clims[1])
                    self.update()
                
            elif event.button==2:
            
                # Middle-click = delete last created point
                if axplot.contains(event)[0] and len(self.offsets)>0 \
                       and self.mode != 'done':
                    self.offsets.pop(-1)
                    if len(self.offsets)==0:
                        self.points.remove()
                    else:
                        self.points.set_offsets(self.offsets)
                    self.update()
                
                # Middle-click on wedge = reset the colour clip
                if cbar.ax.contains(event)[0]:
                    clims = cax.get_clim()
                    cax.set_clim(zmin,zmax)
                    self.update()
                
            if event.button==3:
            
                # Right-click on plot = complete the polygon
                if  axplot.contains(event)[0] and len(self.offsets)>2 \
                       and self.mode != 'done':
                    cpoly = Polygon(self.offsets, animated=False,linewidth=2.0)
                    if self.mode == 'src':
                        cpoly.set_edgecolor('white')
                        cpoly.set_facecolor('none')
                    elif self.mode == 'sky':
                        cpoly.set_edgecolor('yellow')
                        cpoly.set_facecolor('none')
                    elif self.mode == 'exc':
                        cpoly.set_edgecolor('black')
                        cpoly.set_facecolor('none')
                    self.apertures.append(axplot.add_patch(cpoly))
                    self.update('complete')
            
                # Right-click on wedge = halve the upper colour clip
                if cbar.ax.contains(event)[0]:
                    clims = cax.get_clim()
                    cax.set_clim(clims[0],clims[1]/2.0)
                    self.update()
                
        # Store completed polygons in the relevant lists
        def update(self,action=''):

            # When a polygon is complete switch to the next mode
            if action == 'complete':
                if self.mode == 'src':
                    self.offsets_src = self.offsets
                    self.mode = 'sky'
                elif self.mode == 'sky':
                    self.offsets_sky = self.offsets
                    self.mode = 'exc'
                elif self.mode == 'exc':
                    self.offsets_exc = self.offsets
                    self.mode = 'done'
                self.offsets=[]
                
            # When the reset button is complete clear the polygon lists & plot
            elif action == 'reset':
                self.offsets_src = []
                self.offsets_sky = []
                self.offsets_exc = []        
                self.offsets = []
                self.mode = 'src'
                for i in range(len(axplot.collections)):
                    axplot.collections.pop()
                for aperture in self.apertures: aperture.remove()
                self.apertures = []
                self.cent_pt = None
                cax.set_clim(zmin,zmax)
                self.measurements = {}

            # Update the title
            if self.mode == 'src':
                title_str = "Click to outline the source emission."
            elif self.mode == 'sky':
                title_str = 'Click to outline a region of background sky.'
            elif self.mode == 'exc':
                title_str = 'Click to define a source exclusion zone.'
            else:
                title_str = 'Click Reset or Measure when done'
            axplot.set_title(title_str)
            fig.canvas.draw()

        # Reset the plot and clear polygons
        def doreset(self, event):
            self.update('reset')

        # Measure the polygon definitions
        def domeasure(self, event):
            
            if self.offsets_src == [] or self.offsets_sky == []:
                print("Please complete both the source and sky polygons.")
                return
        
            # Measure the source properties
            # [summ,avg,median,stdev,madfm,max,min,npix,centroid_px,cent,centmax]
            # [0   1   2      3     4     5   6   7    8           9    10     ]
            # [radWCentPix,rSum,rSum2,rAmpSum,xSum,xSum2,xAmpSum,ySum,ySum2]
            # [11          12   13    14      15   16    17      18   19   ]
            # [yAmpSum, srcArr, summ2]
            # [20,      21,     22   ]
            print("\nGetting pixels in source polygon ...", end="")
            sys.stdout.flush() 
            src_indices = get_pix_in_poly(self.offsets_src)
            print("done.")
            sys.stdout.flush() 
            print("Measuring properties of the source ...", end="")
            sys.stdout.flush() 
            xi = np.array(src_indices)[:,0]
            yi = np.array(src_indices)[:,1]
            msSrc = measure_pix_props(xydata * scaleUnitsJy, (yi, xi), True)
            print("done.")
            sys.stdout.flush() 
            
            # Calculate the source polygon area and equivalent diameter
            area_px = poly_area(self.offsets_src)
            area_degsq = area_px*(w['pixscale']**2.0)
            diameter_px = 2.0*m.sqrt(area_px/m.pi)
            diameter_deg = diameter_px*w['pixscale']

            # Calculate the brightness weighted diameter in pixels and deg
            diamWtPix = msSrc["radWCentPix"] * 2.0
            diamWtDeg = diamWtPix * w['pixscale']
            
            # Measure the sky properties
            print("Getting pixels in sky polygon ...", end="")
            sys.stdout.flush() 
            sky_indices = get_pix_in_poly(self.offsets_sky)
            print("done.")
            sys.stdout.flush() 
            print("Measuring properties of the sky ...", end="")
            sys.stdout.flush() 
            xi = np.array(sky_indices)[:,0]
            yi = np.array(sky_indices)[:,1]
            msSky = measure_pix_props(xydata * scaleUnitsJy, (yi, xi), True)
            print("done.")
            sys.stdout.flush() 

            # Correlated RMS noise
            diamWcent_asec = diamWtDeg * 3600.0
            if doDivideBeamArea:
                bmaj_asec = header['BMAJ'] * 3600.0
                bmin_asec = header['BMIN'] * 3600.0
                bpa_asec = header['BPA']
                bm_asec = m.sqrt(bmaj_asec*bmin_asec)
                term1 = (4 * bm_asec * bm_asec)/(diamWcent_asec**2.0)
                term2 = m.pow((1.0 + m.pow(bm_asec / diamWcent_asec, 2.0)),
                              -3.0)
                pMADFMskyEff_Jybm = m.sqrt(term1 * term2 * msSky["madfm"]**2.)
                pStdevSkyEff_Jybm = m.sqrt(term1 * term2 * msSky["stdev"]**2.)
            else:
                bm_asec = 1.0
                pMADFMskyEff_Jybm = msSky["madfm"]
                pStdevSkyEff_Jybm = msSky["stdev"]
            
            # Clean bias
            cleanBias_Jybm = cleanBias_sigma * msSky["madfm"]
            
            # Peak pixel corrected for clean bias
            maxPixSrcAbs_Jybm = msSrc["max"] - cleanBias_Jybm
            
            # Calculate the beam area in pixels
            if doDivideBeamArea:
                pix_per_beam = calc_beam_area(header)
            else:
                pix_per_beam = 1.0
                
            # Calculate the integrated flux
            integ_flux = (msSrc["summ"]-(msSky["median"]*msSrc["npix"]))\
                         /pix_per_beam
            
            # Calculate the uncertainty in the integrated flux
            # F.Masci, IPAC: 'Flux-Uncertainty from Aperture Photometry'
            t1 = (msSrc["npix"] * pMADFMskyEff_Jybm)**2.0
            t2 = ( (m.pi * msSrc["npix"]**2.0 * pMADFMskyEff_Jybm**2.0) /
                   (2.0 * msSky["npix"]) )
            dinteg_flux = m.sqrt( t1 + t2) / pix_per_beam
            
            # Calculate the integrated flux corrected for clean bias
            nBeams = max(diamWcent_asec**2.0/bm_asec**2.0, 1.0)
            cleanBiasFlux_Jy = cleanBias_Jybm * nBeams
            integFluxAbs = integ_flux - cleanBiasFlux_Jy            
            
            # Absolute Err in photometric peak (noise + calib err)
            dMaxPixSrcAbs_Jybm = m.sqrt(pMADFMskyEff_Jybm *
                                           pMADFMskyEff_Jybm +
                                           maxPixSrcAbs_Jybm * errCalib_frac * 
                                           maxPixSrcAbs_Jybm * errCalib_frac)
                    
            # Calculate the absolute uncertainty on the integrated flux
            errSumSrc = np.sqrt( (msSrc["srcArr"] * errCalib_frac)**2.0 +
                                 pMADFMskyEff_Jybm**2.0 ).sum()
            t1 = errSumSrc**2.0 
            t2 = ( (m.pi * msSrc["npix"]**2.0 * pMADFMskyEff_Jybm**2.0) /
                   (2.0 * msSky["npix"]) )
            dintegFluxAbs = m.sqrt( t1 + t2) / pix_per_beam
            
            #errSumSrc = np.sqrt( (src_prop_raw[22] * errCalib_frac)**2.0 +
            #                     pMADFMskyEff_Jybm**2.0 ).sum()
            #t1 = errSumSrc**2.0 
            #t2 = ( (m.pi * src_prop_raw[7]**2.0 * pMADFMskyEff_Jybm**2.0) /
            #       (2.0 * sky_prop_raw[7]) )
            #dintegFluxAbs = m.sqrt( t1 + t2) / pix_per_beam
            
            # Calculate the centroid in world coordinates
            [wcentroid_deg] = wcs.wcs_pix2world([msSrc["wcent"]], 0)
            [centroid_deg] = wcs.wcs_pix2world([msSrc["cent"]], 0)
            [centmax_deg] = wcs.wcs_pix2world([msSrc["centmax"]], 0)
            
            # Calculate the error on the X position
            t1 = (msSrc["xSum2"] * pMADFMskyEff_Jybm**2.0 /
                  msSrc["xAmpSum"]**2.0)
            t2 = (msSrc["npix"] * pMADFMskyEff_Jybm**2.0 /
                  msSrc["summ"]**2.0)
            dX_frac = m.sqrt(t1 + t2)
            dX_pix = msSrc["wcent"][0] * dX_frac
            dX_deg = w['pixscale'] * dX_pix

            # Calculate the error on the Y position
            t1 = (msSrc["ySum2"] * pMADFMskyEff_Jybm**2.0 /
                  msSrc["yAmpSum"]**2.0)
            t2 = (msSrc["npix"] * pMADFMskyEff_Jybm**2.0 /
                  msSrc["summ"]**2.0)
            dY_frac = m.sqrt(t1 + t2)
            dY_pix = msSrc["wcent"][1] * dY_frac
            dY_deg = w['pixscale'] * dY_pix
            dCentroid_deg = [dX_deg, dY_deg]
            
            # Absolute positional error:
            errPosition_deg = errPosition_asec / 3600.0
            dXabs_deg = m.sqrt(dX_deg**2.0 + errPosition_deg**2.0)
            dYabs_deg = m.sqrt(dY_deg**2.0 + errPosition_deg**2.0)
            dCentroidAbs_deg = [dXabs_deg, dYabs_deg]
            
            # Calculate the error on the weighted diameter
            t1 = (msSrc["rSum2"] * pMADFMskyEff_Jybm**2.0 /
                  msSrc["rAmpSum"]**2.0)
            t2 = (msSrc["npix"] * pMADFMskyEff_Jybm**2.0 /
                  msSrc["summ"]**2.0)
            dDiamWt_frac = m.sqrt (t1 + t2)
            dDiamWt_deg = diamWtDeg * dDiamWt_frac
            
            # Convert the polygons to world coordinates            
            offsets_src_wld = wcs.wcs_pix2world(self.offsets_src, 0)
            offsets_sky_wld = wcs.wcs_pix2world(self.offsets_sky, 0) 
            try:
                offsets_exc_wld = wcs.wcs_pix2world(self.offsets_exc, 0)
            except Exception:
                offsets_exc_wld = np.array([])

            # Calculate final numbers for the user
            self.measurements['FLUX_JY'] = integ_flux
            self.measurements['FLUXABS_JY'] = integFluxAbs
            self.measurements['dFLUXABS_JY'] = dintegFluxAbs
            self.measurements['dFLUX_JY'] = dinteg_flux
            self.measurements['MAX_JYBM'] = msSrc["max"]
            self.measurements['MAXABS_JYBM'] = maxPixSrcAbs_Jybm
            self.measurements['dMAXABS_JYBM'] = dMaxPixSrcAbs_Jybm
            self.measurements['MIN_JYBM'] = msSrc["min"]
            self.measurements['WCENT_DEG'] = wcentroid_deg
            self.measurements['dWCENT_DEG'] = dCentroid_deg
            self.measurements['dWCENTABS_DEG'] = dCentroidAbs_deg
            self.measurements['NPIXSRC'] = msSrc["npix"]
            self.measurements['BGSTD_JYBM'] = msSky["stdev"]
            self.measurements['BGSTDCOR_JYBM'] = pStdevSkyEff_Jybm
            self.measurements['BGMAD_JYBM'] = msSky["madfm"]
            self.measurements['BGMADCOR_JYBM'] = pMADFMskyEff_Jybm
            self.measurements['BGMED_JYBM'] = msSky["median"]
            self.measurements['BGAVG_JYBM'] = msSky["mean"]
            self.measurements['BGMAX_JYBM'] = msSky["max"]
            self.measurements['BGMIN_JYBM'] = msSky["min"]
            self.measurements['NPIXSKY'] = msSky["npix"]
            self.measurements['CLEANBIAS_JYBM'] = cleanBias_Jybm
            self.measurements['CLEANBIASFLUX_JY'] = cleanBiasFlux_Jy
            self.measurements['CENT_DEG'] = centroid_deg
            self.measurements['CENTMAX_DEG'] = centmax_deg
            self.measurements['DIAM_DEG'] = diameter_deg
            self.measurements['DIAM_WT_DEG'] = diamWtDeg
            self.measurements['dDIAM_WT_DEG'] = dDiamWt_deg
            self.measurements['AREA_DEGSQ'] = area_degsq
            self.measurements['BEAMAREA_PX'] = pix_per_beam
            self.measurements['NBEAMS'] = nBeams
            self.measurements['CALIBERR'] = errCalib_frac
            self.measurements['POLYSRC'] = offsets_src_wld.flatten().tolist()
            self.measurements['POLYSKY'] = offsets_sky_wld.flatten().tolist()
            self.measurements['POLYEXC'] = offsets_exc_wld.flatten().tolist()
            self.measurements['SUMSQ']   = msSrc["summ2"]
            self.measurements['SUM']   = msSrc["summ"]
            
            # Feeback to user
            print("")
            keys = ['SUMSQ',
                    'SUM',
                    'FLUX_JY', 
                    'FLUXABS_JY', 
                    'dFLUXABS_JY', 
                    'dFLUX_JY', 
                    'MAX_JYBM', 
                    'MAXABS_JYBM', 
                    'dMAXABS_JYBM', 
                    'MIN_JYBM', 
                    'WCENT_DEG', 
                    'dWCENT_DEG', 
                    'dWCENTABS_DEG', 
                    'NPIXSRC', 
                    'BGSTD_JYBM', 
                    'BGSTDCOR_JYBM', 
                    'BGMAD_JYBM', 
                    'BGMADCOR_JYBM', 
                    'BGMED_JYBM', 
                    'BGAVG_JYBM', 
                    'BGMAX_JYBM', 
                    'BGMIN_JYBM', 
                    'NPIXSKY', 
                    'CLEANBIAS_JYBM',
                    'CLEANBIASFLUX_JY',
                    'CENT_DEG', 
                    'CENTMAX_DEG', 
                    'DIAM_DEG', 
                    'DIAM_WT_DEG', 
                    'dDIAM_WT_DEG', 
                    'AREA_DEGSQ', 
                    'BEAMAREA_PX', 
                    'NBEAMS', 
                    'CALIBERR']
            for key in keys:
                 #print(key, "        \t=\t", self.measurements[key])
                print("{:<20s} = {:}".format(key, self.measurements[key]))
            print("")

            # Plot the centroid of the pixels
            self.pcentre = RegularPolyCollection(
                10,
                rotation=0,
                sizes=(50,),
                facecolors = 'white',
                edgecolors = 'black',
                linewidths = (1,),
                offsets = [msSrc["wcent"]],
                transOffset = axplot.transData
                )
            axplot.add_collection(self.pcentre)
            self.update()
            
        # Save the measurements to a file
        def dosavefile(self, event):

            # Strip the '.fits' to form the datfile name
            suffix = '.fits'
            if filename.endswith(suffix):
                datfilename = filename[:-len(suffix)] + '.polyflux.dat'
            else:
                datfilename = filename + '.polyflux.dat'

            if not self.measurements:
                print("Measurements empty! Click measure before saving.")
                return

            # Open the datfile for appending
            if os.path.exists(datfilename):
                print("\nAppending measurement to file '%s' ... " %
                      datfilename, end="")
            else:
                print("\nSaving measurement to NEW file '%s' ... " %
                datfilename, end="")
            sys.stdout.flush()
            datHandle = open(datfilename, 'a')
            datHandle.write('START\n')
            for key,val in self.measurements.iteritems():
                datHandle.write("%s=%s\n" % (key,val))
            datHandle.write('END\n\n')
            datHandle.close()
            print("done.")

        # Plot pre-existing apertures
        def plot_apertures(self):
            if not self.offsets_src == []:
                cpoly = Polygon(self.offsets_src, animated=False,linewidth=2.0)
                cpoly.set_edgecolor('white')
                cpoly.set_facecolor('none')
                self.apertures.append(axplot.add_patch(cpoly))
                self.mode = 'sky'
                
            if not self.offsets_sky == []:
                cpoly = Polygon(self.offsets_sky, animated=False,linewidth=2.0)
                cpoly.set_edgecolor('yellow')
                cpoly.set_facecolor('none')
                self.apertures.append(axplot.add_patch(cpoly))
                self.mode = 'exc'
            self.update()
            
#-----------------------------------------------------------------------------#

    # Make an instance of the poly-editor and bind events to methods
    editor = ThreePolyEditor()
    fig.canvas.mpl_connect('button_press_event', editor.onclick)
    breset.on_clicked(editor.doreset)
    bmeasure.on_clicked(editor.domeasure)
    bsave.on_clicked(editor.dosavefile)

    # Load up the old polygons if they exist
    if not (oldSrcPoly == [] and oldSkyPoly == []):

        editor.offsets_src = wcs.wcs_world2pix(polylst_parse(oldSrcPoly),0)
        editor.offsets_sky = wcs.wcs_world2pix(polylst_parse(oldSkyPoly),0)
        editor.plot_apertures()

    # Draw the plot to the screen
    show()


#-----------------------------------------------------------------------------#
def kvis_ann_parse(filename):
    """
    Parse the ellipses from a Kis .ann file
    """
    
    coord_list=[];
    ell = re.compile('^ELLIPSE\s+W\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')
    
    ANNFILE = open(filename, "r")
    for line in ANNFILE:
        entry={}
        line = line.rstrip("\n\r")
        mch = ell.match(line)
        if mch:
            entry = (float(mch.group(1)),float(mch.group(2)))
            coord_list.append(entry)
    return coord_list

        
#-----------------------------------------------------------------------------#
def datfile_read(filename):
    """
    Read a configuration file and output a 'KEY=VALUE' dictionary
    """

    # Open the file for reading & initialise variables 
    print("Reading ", filename)
    measurements = []               # List of measurements (param_tables)
    param_table = dict()            # Dictionary to hold keyword-value pairs
    DATFILE = open(filename, "r")

    # Compile a few useful regular expressions
    spaces = re.compile('\s+')
    comma_and_spaces = re.compile(',\s+')
    comma_or_space = re.compile('[\s|,]')
    brackets = re.compile('[\[|\]\(|\)|\{|\}]')
    comment = re.compile('#.*')
    quotes = re.compile('\'[^\']*\'')
    key_val = re.compile('^.+=.+')
    polykey = re.compile('^POLY\.+')
    centrkey = re.compile('^CENT\.+')

    # Read in the input file, line by line
    for line in DATFILE:
        line = line.rstrip("\n\r")

        # Filter for comments and blank lines
        if not comment.match(line) and key_val.match(line):
            # Weed out internal comments & split on 1st '='
            line = comment.sub('',line)
            (keyword, value) = line.split('=',1)
            keyword = keyword.strip()          # kill external whitespace
            keyword = spaces.sub('', keyword)  # kill internal whitespaces
            value = brackets.sub('',value)     # kill brackets
            value = value.strip()              # kill external whitespace
            value = spaces.sub(' ', value)     # shrink internal whitespace
            value = value.replace("'", '')     # kill quotes
            value = comma_and_spaces.sub(',', value) # kill ambiguous spaces

            # If the line contains a value
            if value:
                valuelist = comma_or_space.split(value)
                param_table[keyword] = valuelist

        # If 'END' is encountered then append param_table to measurements
        if(line == 'END') :
            measurements.append(param_table)
            param_table={}
    return measurements
    

#-----------------------------------------------------------------------------#
def usage():
    """
    Print usage instructions and exit.
    """

    print("\n\tUSAGE: polyphot.py <file.fits>\n")
    sys.exit(1)


#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    main()

