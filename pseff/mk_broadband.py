#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mk_broadband_psf
----------------
Create a series of broadband PSFs given input various SED and monochromatic
PSF files as well as an instrumental thruput

"""
from __future__ import print_function, division

import os
import time
import argparse
import configparser
import astropy.units as u
import pseff


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('psf_cube_files', type=str, nargs='+',
                        help='Monochromatic PSF FITS files')
    parser.add_argument('sed_files', type=str, nargs='+',
                        help='Source energy distribution of stars')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help="Configuration file")

    return parser.parse_args()


def extract_config(configfile):
    """
    Extract typed information from config parser or return default

    Parameters
    ----------
    configfile: str
        Path to configuration text file

    Returns
    -------
    dictionary of useful parameters

    """
    confparams = configparser.SafeConfigParser()
    confparams.read(configfile)

    # PSF
    n_wl = confparams.getint('psf', 'n_wavelength')
    psf_size = confparams.getint('psf', 'size')
    cube_shape = (n_wl, psf_size, psf_size)
    pixel_size = confparams.getfloat('psf', 'pixel_size')
    psf_wlunit = getattr(u, confparams.get('psf', 'wlunit').lower())

    psf_dict = dict(shape=cube_shape, pixel_size=pixel_size, wlunit=psf_wlunit)

    # Telescope
    thruput_filename = confparams.get('thruput', 'filename')
    thruput_cols = [confparams.getint('thruput', 'wl_col'),
                    confparams.getint('thruput', 'thruput_col')]
    thruput_wlunit = getattr(u, confparams.get('thruput', 'wlunit').lower())

    thruput_dict = dict(filename=thruput_filename, cols=thruput_cols,
                        wlunit=thruput_wlunit)

    # SED
    sed_type = confparams.get('sed', 'type')
    sed_wlunit = getattr(u, confparams.get('sed', 'wlunit').lower())

    sed_dict = dict(wlunit=sed_wlunit)

    # Broadband
    wl_step = confparams.getfloat('broadband', 'wl_step')
    niemi = confparams.getboolean('broadband', 'use_niemi')
    aocs = confparams.getboolean('broadband', 'use_aocs')
    prf = confparams.getboolean('broadband', 'use_prf')
    down_fact = confparams.getfloat('broadband', 'down_factor')

    broadband_dict = dict(wl_step=wl_step, niemi=niemi, aocs=aocs, prf=prf)

    outputfile = os.path.join(confparams.get('output', 'path'),
                              confparams.get('output', 'nameformat'))

    return dict(psf=psf_dict,
                thruput=thruput_dict,
                sed=sed_dict,
                sed_type=sed_type,
                broadband=broadband_dict,
                down_fact=down_fact,
                outputfile=outputfile)


def main():
    args = parse_args()

    conf = extract_config(args.config)

    start_time = time.time()

    for psf_file in args.psf_cube_files:
        psf = pseff.PSFcube.from_fits(psf_file, **conf['psf'])
        thruput = pseff.TelescopeThroughput.from_file(**conf['thruput'])
        for sed_file in args.sed_files:
            if conf['sed_type'].lower() == 'fits':
                sed = pseff.SED.from_fits(sed_file, **conf['sed'])
            else:
                sed = pseff.SED.from_ascii(sed_file, **conf['sed'])

            psf.to_broadband(sed=sed, thruput=thruput, **conf['broadband'])

            psf_basename = os.path.basename(os.path.splitext(psf_file)[0])
            sed_basename = os.path.basename(os.path.splitext(sed_file)[0])
            outfile = conf['outputfile'].format(psf=psf_basename,
                                                sed=sed_basename)
            psf.write_effective_psf(outfile)

            if conf['down_fact'] > 1:
                psf.downsample_detector(conf['down_fact'])
                psf.write_effective_psf(outfile.replace('.fits',
                                                        '_downsampled.fits'))

    dtime = int(time.time() - start_time)
    time_text = ("\n---\n"
                 "Elapsed time: {0:02d}h{1:02d}m{2:02d}s")
    print(time_text.format(dtime // 3600, (dtime % 3600) // 60, dtime % 60))


if __name__ == '__main__':
    main()