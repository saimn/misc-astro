#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Merge the matched catalogs for one field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Merge the matched catalogs for stamp & dirac, NUV and FUV.
The merged catalog is based on the stamp NUV catalog for all common field.

The variables below :data:`NFSD_COLS`, :data:`NF_COLS` and :data:`MIX_COLS`
define the columns that will be copied and how they will be renamed.

Code that was used to test the routine versus a catalog generated with
``merge_catalogs.pro``:

::

    os.chdir(os.path.expanduser('XMMLSS/XMMLSS_00/OUTPUT_XMMLSS_00_DEEP'))

    from astropy.io.fits import getdata
    test = getdata('XMMLSS_00-xd-prior_v3-merged-cat.fits', 1)
    ref = getdata('XMMLSS_00-xd-prior_v3-merged-cat.fits.bak', 1)

    ref_names = np.array(ref.columns.names)
    test_names = np.array(test.columns.names)
    iref, itest = matchid(ref_names, test_names)

    print "Non matching columns in test: {}".format(test_names[itest == False])
    print "Non matching columns in reference : {}".format(ref_names[iref == False])

    max_length = max(np.concatenate((test_names, ref_names)), key=len)
    for col in test_names[itest]:
        if col not in ['PRIOR_CAT', 'CFHTLS_CAT']:
            print " - {} : {}".format(col, np.max(np.abs(ref[col] - test[col])))

"""

from __future__ import absolute_import

import argparse
import os
import numpy as np
import sys
from astropy.io import fits

from merge_stats import stats

#: Default EM mag/flux/snr columns
MAG_COLS = ['FLUX_EM', 'MAG_EM', 'MAG_ERR_EM', 'SNR_EM']

#: Columns that will be copied with ``NUV_`` / ``FUV_`` prefix and ``_STAMP``
#: / ``_DIRAC`` suffix. The ``_EM`` suffix will be removed.
NFSD_COLS = MAG_COLS + ['FLUX_ERR_EM', 'BACKGD_EM']

#: Columns that will be copied with ``NUV_`` / ``FUV_`` prefix, using stamp
#: inputs.
NF_COLS = ['FLAG', 'NC_R1', 'NC_R2', 'NC_R3', 'NN_DIST', 'DX_MAXLIKE',
           'DY_MAXLIKE', 'TILE_LIKHD', 'TILE_MNBKGD', 'TILE_NOPT',
           'TILE_QUIETS_RATE', 'TILE_TTLFLUX']

#: Columns that will be copied with ``NUV_`` / ``FUV_`` prefix, using dirac
#: inputs.
MIX_COLS = ['MAG_MIX', 'MAG_ERR_MIX', 'FLAG_MIX', 'NFLAG_MIX', 'IDFLAG_MIX']

#: Columns that should not be copied
SKIP_COLS = ['IDENT', 'TILEX', 'TILEY', 'FLUX_EM_DBIAS', 'MAG_EM_DBIAS',
             'MAG_ERR_EM_DBIAS', 'SNR_EM_DBIAS']


def get_column_values(zeros, indices, data, column):
    """ Get a column's data filled inside another one.

    Return an array with the size and default values from `zeros`, with the
    values from `data[column]` put at `indices`.
    """

    coldat = zeros.copy()
    coldat[indices] = data.field(column)
    return coldat


def readfits(filename, ext=1):
    if filename:
        return fits.getdata(filename, ext)
    else:
        return None


def merge_catalogs(stamp_nuv, outname, stamp_fuv=None, dirac_nuv=None,
                   dirac_fuv=None, debug=False):
    "Merge the 4 catalogs"

    # use stamp_nuv as reference
    header = fits.getheader(stamp_nuv, 0)

    # read catalogs
    snuv = readfits(stamp_nuv)
    sfuv = readfits(stamp_fuv)
    dnuv = readfits(dirac_nuv)
    dfuv = readfits(dirac_fuv)

    cols = snuv.columns
    use_dbias = False
    bias_cols = [col + '_DBIAS' for col in MAG_COLS]

    if set(bias_cols).issubset(set(cols.names)):
        print "Using _DBIAS columns"
        use_dbias = True

    if not stamp_fuv and not dirac_fuv:
        SKIP_COLS.extend(['FUV_MAG_GALEX',
                          'FUV_MAGERR_GALEX',
                          'FUV_MAG_APER7_GALEX',
                          'FUV_MAGERR_APER7_GALEX',
                          'FUV_ALPHA_J2000_GALEX',
                          'FUV_DELTA_J2000_GALEX',
                          'FUV_X_IMAGE_GALEX',
                          'FUV_Y_IMAGE_GALEX',
                          'FUV_CN_MAG_GALEX',
                          'FUV_SKYBG_GALEX'])

    ident = snuv['IDENT']
    if sfuv is not None:
        ident = np.union1d(ident, sfuv['IDENT'])

    # match the catalogs with the common ident array
    ind_snuv = np.in1d(ident, snuv['IDENT'])

    if sfuv is not None:
        ind_sfuv = np.in1d(ident, sfuv['IDENT'])

    if dnuv is not None:
        ind_dnuv = np.in1d(ident, dnuv['IDENT'])

    if dfuv is not None:
        ind_dfuv = np.in1d(ident, dfuv['IDENT'])

    nbpriors = len(ident)
    outcols = [ident]
    outnames = [('IDENT', np.dtype(int))]

    def add_column(name, dtype, values):
        if debug:
            print "add column : {0} ({1}, {2})".format(name, dtype, values[0])
        outcols.append(values)
        outnames.append((name, dtype))

    for i, col in enumerate(cols.names):
        dtype = snuv.field(col).dtype

        # create empty columns with the good type and add -99 for non string
        # columns
        zeros = np.zeros(nbpriors, dtype=dtype)
        if dtype.kind != 'S':
            zeros -= 99

        # add columns and fill it with the good inputs
        if col in NFSD_COLS:
            c = col.replace('_EM', '')

            snuv_col = col
            if use_dbias and col in MAG_COLS:
                snuv_col = col + '_DBIAS'

            add_column('NUV_' + c + '_STAMP', dtype,
                       get_column_values(zeros, ind_snuv, snuv, snuv_col))

            if dnuv is not None:
                add_column('NUV_' + c + '_DIRAC', dtype,
                           get_column_values(zeros, ind_dnuv, dnuv, col))

            if sfuv is not None:
                add_column('FUV_' + c + '_STAMP', dtype,
                           get_column_values(zeros, ind_sfuv, sfuv, col))

            if dfuv is not None:
                add_column('FUV_' + c + '_DIRAC', dtype,
                           get_column_values(zeros, ind_dfuv, dfuv, col))
        elif col in NF_COLS:
            add_column('NUV_' + col, dtype,
                       get_column_values(zeros, ind_snuv, snuv, col))

            if sfuv is not None:
                add_column('FUV_' + col, dtype,
                           get_column_values(zeros, ind_sfuv, sfuv, col))
        elif col in MIX_COLS:
            add_column('NUV_' + col, dtype,
                       get_column_values(zeros, ind_dnuv, dnuv, col))

            if dfuv is not None:
                add_column('FUV_' + col, dtype,
                           get_column_values(zeros, ind_dfuv, dfuv, col))
        elif col not in SKIP_COLS:
            add_column(col, dtype,
                       get_column_values(zeros, ind_snuv, snuv, col))

    outarr = np.rec.fromarrays(outcols, dtype=outnames)

    if sfuv is not None:
        common_cols = set(cols.names) - set(NFSD_COLS + NF_COLS +
                                            MIX_COLS + SKIP_COLS)

        # Fill common columns for priors used if FUV but not in NUV. Currently
        # priors for stamp and dirac are the same for a given band, so we just
        # need to check stamp_fuv

        # Find sources that are in the fuv cat but not in the nuv one
        ind = ~np.in1d(sfuv['IDENT'], snuv['IDENT'])

        # For these sources, copy the values of common cols to outarr
        if ind.any():
            id_sfuv = np.in1d(ident, sfuv['IDENT'][ind])
            for col in common_cols:
                outarr[col][id_sfuv] = sfuv[col][ind]

    # Check that it is now ok
    if outarr['ALPHA'].min() == -99 or outarr['DELTA'].min() == -99:
        print ("Some sources have a radec = -99 which should not happen."
               "Check this !")

    # Find objects detected in at least one of the 4 catalogs
    # cols = [col for col in ['NUV_MAG_STAMP', 'NUV_MAG_DIRAC',
    #                         'FUV_MAG_STAMP', 'FUV_MAG_DIRAC']
    #         if col in outarr.dtype.names]

    # mix = np.array([outarr.field(col) > 0 for col in cols])

    # detected = mix.any(axis=0)
    # print "{0} objects detected in at least one catalog".format(
    #     np.count_nonzero(detected))

    # tbhdu = fits.BinTableHDU(outarr[detected])

    tbhdu = fits.BinTableHDU(outarr)

    if debug:
        d = tbhdu.data[0]
        for i, n in enumerate(d.array.names):
            print '{0} : {1}'.format(n, d[i])

    # compute the number of detected sources in each catalog
    header.set('NNBSTAMP', np.count_nonzero(outarr.field('NUV_MAG_STAMP') > 0),
               'Number of NUV sources detected with stamp', after='DATE')

    if dnuv is not None:
        header.set('NNBDIRAC',
                   np.count_nonzero(outarr.field('NUV_MAG_DIRAC') > 0),
                   'Number of NUV sources detected with dirac', after='DATE')

    if sfuv is not None:
        header.set('FNBSTAMP',
                   np.count_nonzero(outarr.field('FUV_MAG_STAMP') > 0),
                   'Number of FUV sources detected with stamp', after='DATE')

    if dfuv is not None:
        header.set('FNBDIRAC',
                   np.count_nonzero(outarr.field('FUV_MAG_DIRAC') > 0),
                   'Number of FUV sources detected with dirac', after='DATE')

    print "Writing file " + outname
    hdu = fits.PrimaryHDU(header=header)
    hdulist = fits.HDUList([hdu, tbhdu])
    hdulist.writeto(outname, clobber=True)


def main():
    "main program"

    parser = argparse.ArgumentParser(
        description='Merge the matched catalogs for one field')
    parser.add_argument('field', help='Name of the field')
    parser.add_argument('path', help='Path of the catalogs')
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Show all message, including debug messages')
    parser.add_argument('-s', '--stats', action='store_true',
                        help='Compute some stats on the merged catalog')

    args = parser.parse_args()

    files = {'stamp_nuv': '-nd-prior_v3-bgfromdiff',
             'stamp_fuv': '-fd-prior_v3-bgfromdiff',
             'dirac_nuv': '-nd-prior_v3-dirac',
             'dirac_fuv': '-fd-prior_v3-dirac'}

    for k, v in files.iteritems():
        files[k] = os.path.join(args.path,
                                args.field + v + '-EM_PRIOR_GALEX.fits')

    files['outname'] = os.path.join(args.path, args.field +
                                    '-xd-prior_v3-merged-cat.fits')

    merge_catalogs(debug=args.debug, **files)

    if args.stats:
        stats(files['outname'])


if __name__ == '__main__':
    status = main()
    sys.exit(status)
