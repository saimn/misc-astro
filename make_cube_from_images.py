# -*- coding: utf-8 -*-

from __future__ import absolute_import, print_function

import argparse
import numpy as np
from astropy.io import fits
from mpdaf.obj import Cube


def make_cube_from_images(flist, outf):
    data = [np.copy(fits.getdata(f, ext=1)) for f in flist]
    shape = np.max(np.array([d.shape for d in data]), axis=0)
    for d in data:
        d.resize(shape)

    cube = Cube(data=np.dstack(data), copy=False)
    cube.write(outf, savemask='nan')


def main():
    parser = argparse.ArgumentParser(
        description='Create a cube from a list of images.')
    parser.add_argument('output_cube', help='Output cube (FITS file).')
    parser.add_argument('images', nargs='*', help='Input images (FITS files).')
    args = parser.parse_args()
    make_cube_from_images(args.flist, args.output_cube)


if __name__ == '__main__':
    main()
