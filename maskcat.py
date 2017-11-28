#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Reprocess a catalog to flag objects.

Read the GALEX and CFHTLS region filenames from the header and it to recompute
the FLAG column. The FLAG column is the combinaison of the following values:

- 1 for GALEX masks
- 2 for CFHTLS masks
- 4 for bad stamps
- 8 for CFHTLS conservative masks (big halos, color cyan or magenta in reg maks files)

For the binary mask:

- 128 = outside of tiles
- 2 for CFHTLS masks
- 1 for GALEX masks

Test data::

    basepath = '/data/v3.0'
    fieldpath = basepath + '/XMMLSS/XMMLSS_01'
    catpath = fieldpath + '/OUTPUT_XMMLSS_01_W1_T0007_r984/XMMLSS_01-nd-prior_v3-bgfromdiff-flux.out.fits'
    galpath = fieldpath + '/INPUT_XMMLSS_01/'
    optpath = basepath + '/CFHTLS_T07/W1/'

"""

import numpy as np
import os
import pyregion
import shutil
import sys
import time

from argh import ArghParser
from astropy import wcs
from astropy.io import fits
from PIL import Image, ImageDraw
from scipy import ndimage

from fitsdate import fitsdate

def timeit(method):
    "Decorator to measure the execution time of a function."

    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()

        if kw.get('debug', False):
            # print '%r (%r, %r) %2.2f sec' % \
            #     (method.__name__, args, kw, te-ts)
            print '%r - %2.2f sec' % (method.__name__, te-ts)
        return result

    return timed


@timeit
def readfile(filename, debug=False):
    print "Reading " + filename + " ...",
    sys.stdout.flush()
    reg = pyregion.open(filename)
    print "{0} regions".format(len(reg))
    return reg


@timeit
def regions_mask(regions, header, structure, dilate=False, debug=False, colors=None):
    print "Masking with {0} regions".format(len(regions))
    shape = (header.get('NAXIS1'), header.get('NAXIS2'))

    # Create the image mask. With PIL (0,0) is the upper left corner.
    img = Image.new('L', shape, 0)
    w = wcs.WCS(header)

    for reg in regions:
#        print(reg.attr[1].get("color"))
        if (not colors or (reg.attr[1].get("color") in colors) ):
#            print(reg.attr[1].get("color"))
            if reg.name == 'polygon':
                world = np.array([reg.coord_list[::2], reg.coord_list[1::2]])
                pix = w.wcs_sky2pix(world.transpose(), 0).flatten() + 0.5
                # +0.5 as ImageDraw use integer value as left/bottom borders (same
                # as IDL's polyfillv)
                # print pix.tolist()
                ImageDraw.Draw(img).polygon(pix.tolist(), outline=1, fill=1)
            elif reg.name == 'circle':
                rad  =  reg.coord_list[2]
                xy   =  np.array([ [reg.coord_list[0]+rad, reg.coord_list[0]-rad],
                                   [reg.coord_list[1]-rad, reg.coord_list[1]+rad] ])
                xypix = w.wcs_sky2pix(xy.transpose(), 0).flatten() + 0.5
#                print xypix.tolist()
                ImageDraw.Draw(img).ellipse(xypix.tolist(), outline=1, fill=1)            
            elif reg.name == 'box':
                dist =  np.array([ reg.coord_list[2],reg.coord_list[3] ])/2.0
                xy   =  np.array([ [reg.coord_list[0]+dist[0], reg.coord_list[0]-dist[0]],
                                   [reg.coord_list[1]-dist[1], reg.coord_list[1]+dist[1]] ])
                xypix = w.wcs_sky2pix(xy.transpose(), 0).flatten() + 0.5
#                print xypix.tolist()
                ImageDraw.Draw(img).rectangle(xypix.tolist(), outline=1, fill=1)            

            else:
                print "Region type not managed: " + reg.name

    mask = np.array(img)
    # mask = regions.get_mask(header=header, shape=shape)
    if dilate:
        mask = ndimage.binary_dilation(mask, structure=structure)

    return mask.astype(np.bool)


def prior_mask(priorfile, imhdr, struct):
    "Create a binary mask with the prior limits."

    opthdr = fits.getheader(priorfile, 1)
    if opthdr.get('RA_MIN'):
        ramin, ramax, decmin, decmax = [opthdr.get(k) for k in
                                        ['RA_MIN', 'RA_MAX', 'DEC_MIN', 'DEC_MAX']]
        opt_limits = [ramin, decmin, ramin, decmax, ramax, decmax,
                      ramax, decmin, ramin, decmin]
        opt_limits = ','.join([str(o) for o in opt_limits])
        region = pyregion.parse('fk5;polygon({0})'.format(opt_limits))
        mask = regions_mask(region, imhdr, struct, dilate=False)
    else:
        mask = None

    return mask


def binary_disc(radius):
    x, y = np.arange(2 * radius + 1), np.arange(2 * radius + 1)
    dist = np.sqrt((x - radius) ** 2 + (y[:, np.newaxis] - radius) ** 2)
    return (dist <= radius).astype(int)


def diff(new_file, orig_file):
    """Compute the diff images for the CFHTLS and GALEX masks."""

    header = fits.getheader(new_file)
    new = fits.getdata(new_file)
    orig = fits.getdata(orig_file)

    diff = (new & 1).astype(float) - (orig & 1).astype(float)
    fits.writeto(new_file.replace('.fits', '-diff-galex.fits'), diff, header,
                 clobber=True)

    diff = (new & 2).astype(float) - (orig & 2).astype(float)
    fits.writeto(new_file.replace('.fits', '-diff-cfhtls.fits'), diff, header,
                 clobber=True)


def backup_file(filename):
    """Copy the file to keep a backup.

    The backup has the same name as the original file with ".bak" appended. If
    "file.bak" already exists then "file.bak.1" is used, and so on. This comes
    from `pyfits.hdulist`.

    """
    backup = filename + '.bak'
    idx = 1
    while os.path.exists(backup):
        backup = filename + '.bak.' + str(idx)
        idx += 1

    if os.path.exists(filename):
        try:
            shutil.copy(filename, backup)
        except IOError, e:
            raise IOError('Failed to save backup to destination %s: '
                          '%s' % (filename, str(e)))


def flag(catpath, galpath=None, optpath=None,simu=False,debug=False):
    """Flag a catalog with DS9 region files and reprocess the mask.

    Read regions files to recompute the FLAG column and the mask file. The
    catalog and image files (mask, masked residuals) are overwritten but a
    backup of the original files is saved before.

    """

    backup_file(catpath)
    maskpath = catpath.replace('flux.out.fits', 'mask.fits')
    backup_file(maskpath)
    # for simu
    if simu:
        maskpath1 =  maskpath.replace('total-mask.fits', '1-mask.fits')
        try:
            shutil.copy(maskpath1, maskpath)
        except IOError, e:
            raise IOError('Failed to copy mask file from simu %s: '
                          '%s' % (maskpath1, str(e)))

    hdulist = fits.open(catpath, mode='update')
    cat = hdulist[1].data
    hdr = hdulist[1].header
       
    galpath = galpath or hdr['GALPATH']
    optpath = optpath or hdr['OPTPATH']
    
    if not simu:
        galimg =  hdr['GALIMG']
    else:
        # get the original image (not the simu one)
        outname = hdr['OUTNAME']
        simuname  = hdr['SIMUNAME']
        outpath = os.path.dirname(catpath)
        ref = os.path.join(outpath,outname.rstrip(simuname)+'-flux.out.fits')
        refhdulist = fits.open(ref)
        refhdr = refhdulist[1].header
        galimg = refhdr['GALIMG']
        refhdulist.close()

    imhdr = fits.getheader(os.path.join(galpath,galimg))
    w = wcs.WCS(imhdr)
    shape = (imhdr.get('NAXIS1'), imhdr.get('NAXIS2'))

    # convert world coords to pixel
    world = np.array([cat['RA'], cat['DEC']])
    pix = w.wcs_sky2pix(world.transpose(), 0)
    ipix = pix.round().astype(int)
    x, y = ipix[:, 0], ipix[:, 1]

    # prepare the structuring element: disc with fwhm radius
    struct = binary_disc(int(np.ceil(hdr['FWHM'])))
 
    if hdr.get('GALMASK', ''):
        gregfile = hdr['GALMASK']
    else: 
        gregfile = galimg.replace('.fits', '.mask.reg')
            
            
    gmskused = [] 
    fullgregfile = os.path.join(galpath,gregfile)
    if os.path.isfile(fullgregfile):
        gmskused = gregfile
        reglist = readfile(fullgregfile, debug=debug)
        galmask = regions_mask(reglist, imhdr, struct, dilate=True)
        galmasked = galmask[y, x]

        # remove the existing flag and add the new one
        cat['FLAG'][:] &= ~1 # reset all bit 0
        cat['FLAG'][galmasked] |= 1 # set bit 0 for galmasked object
    else:
        print "Region file {0} not found".format(fullgregfile)
                        
    optcat = hdr.get('PRIORLST').split(',')
    if hdr.get('MASKLST', ''):
        optreg = hdr['MASKLST'].split(',')
    elif hdr.get('OPTMASK', ''):
        optreg = hdr['OPTMASK'].split(',')
    else:
        optreg = [f.replace('.prior.fits', '.mask.reg') for f in optcat]
        optreg = [f.replace('.fits', '.reg') for f in optreg]

    # default mask (all masked)
    optmask = np.ones(shape, dtype=np.bool)
    optmask2 = np.ones(shape, dtype=np.bool)
    omskused = []
    
    # this loop will unmask good regions
    for opt, reg in zip(optcat, optreg):
        regfile = os.path.join(optpath, reg.strip())
        optfile = os.path.join(optpath, opt.strip())

        if not os.path.isfile(regfile):
            print "Region file {0} not found".format(reg)
            regmask = None
            regmask2 = None
        else:
            omskused.append(reg)
            reglist = readfile(regfile, debug=debug)
            regmask  = regions_mask(reglist, imhdr, struct, dilate=True, colors='green')
            # read optional big masks
            regmask2 = regions_mask(reglist, imhdr, struct, dilate=True, colors=['cyan','magenta'])

        priormask = prior_mask(optfile, imhdr, struct)
        
        if priormask != None:
            if regmask != None:
                regmask = ~priormask | regmask
                regmask2 = ~priormask | regmask2
            else:
                regmask = ~priormask
                regmask2 = ~priormask
        
        if regmask != None: # priormask = None
            optmask  &= regmask
            optmask2 &= regmask2
 
    # correct the default mask if nothing was masked
    #  (ie accept all instead of reject all)       
    if optmask.all():  
        optmask = np.zeros(shape, dtype=np.bool)
    if optmask2.all():  
        optmask2 = np.zeros(shape, dtype=np.bool)
        
    # Python is 0-indexed and the indices have the slowest axis first and fast
    # axis last, i.e. for a 2-D image, the fast axis (X-axis) which
    # corresponds to the FITS NAXIS1 keyword, is the second index
    optmasked = optmask[y, x]
    optmasked2 = optmask2[y, x]

    # remove the existing flag and add the new one
    cat['FLAG'][:] &= ~2 # reset all bit 1
    cat['FLAG'][optmasked] |= 2 # set bit 1 for optmasked object
    cat['FLAG'][:] &= ~8 # reset all bit 3
    cat['FLAG'][optmasked2] |= 8 # set bit 3 for optmasked object

    # Recompute the mask
    hdumask = fits.open(maskpath, mode='update')
    mask = hdumask[0].data
    maskhdr = hdumask[0].header
    mask &= 128    # keep only the tiles mask

    # eventually create/update the GMSKUSED keyword
    if gmskused:
        hdr.set('GMSKUSED',gmskused,'galex mask used',after='DATE')
        maskhdr.set('GMSKUSED',gmskused,'galex mask used',after='DATE')

        # add the new galmask
        try:
            mask += galmask.astype(int)
        except UnboundLocalError: 
            pass

    # eventually create the OMSKUSED keyword
    if omskused:
        hdr.set('OMSKUSED',','.join(omskused),'Optical masks used',after='DATE')
        maskhdr.set('OMSKUSED',','.join(omskused),'Optical masks used',after='DATE')

        # add the new optmask (green only)
        try:
            mask += 2 * optmask.astype(int)
        except UnboundLocalError: 
            pass
    
        # add the extended optmask (magenta &  cyan)
        try:
            mask += 8 * optmask2.astype(int)
        except UnboundLocalError: 
            pass

        
    # Recompute the processed surface and update the keyword
    surface = np.count_nonzero((mask & (7+128)) == 0) * imhdr['CDELT1'] ** 2
    maskhdr['SURFACE'] = surface
    hdr['SURFACE'] = surface
    print "Non masked surface: {0} square degree".format(surface)

    # compute the conservative surface (with cyan and magenta colors)
    surface2 = np.count_nonzero((mask & (15+128)) == 0) * imhdr['CDELT1'] ** 2
    hdr.set('SURFACE2',surface2,'surface excluding all flags',after='SURFACE')
    maskhdr.set('SURFACE2',surface2,'surface excluding all flags',after='SURFACE')
        
    # Recompute the masked residual
    if not simu:
        imask = ((mask & (7+128))== 0).astype(int)
        diffpath = catpath.replace('flux.out.fits', 'diff.fits')
        diff = fits.getdata(diffpath)
        
        bkg = fits.getdata(os.path.join(galpath, hdr['GALBKG']))

        diff_masked = diffpath.replace('diff.fits', 'diff-masked.fits')
        diff_bkg_masked = diffpath.replace('diff.fits', 'diff-bkg-masked.fits')

        if os.path.isfile(diff_masked):
            backup_file(diff_masked)
        if os.path.isfile(diff_bkg_masked):
            backup_file(diff_bkg_masked)

        fits.writeto(diff_masked, diff * imask, maskhdr, clobber=True)
        fits.writeto(diff_bkg_masked, (diff - bkg) * imask, maskhdr, clobber=True)

    # flush the new flaged catalog on disk
    hdulist.flush()

    # flush the new mask on disk
    hdumask.flush()

    return {'fitsdate':fitsdate(catpath)}

if __name__ == '__main__':
    p = ArghParser()
    p.add_commands([diff, flag])
    p.dispatch()
