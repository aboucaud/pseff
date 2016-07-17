#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Broadband PSF
=============
Convert a series of monochromatic PSF images into an instrumental
broadband PSF given a stellar spectral energy distribution

"""
from __future__ import print_function, division

import os
import copy
import astropy.units as u
import astropy.io.fits as pyfits
import numpy as np
import numpy.fft as npf
import scipy.ndimage as spn
import scipy.interpolate as spi


class SED(object):
    """Stellar Spectral Energy Distribution class"""

    def __init__(self, wavelength, spectrum, wlunit, flambda=True):
        """
        Parameters
        ----------
        wavelength: array_like
            Array of wavelength values corresponding to the SED
        spectrum: array_like
            Array of SED values
        wlunit: `astropy.unit`
            Unit of the given wavelength
        flambda: bool, optional
            Format of SED values (default is `True`)
                * `True` corresponds to energy flux per unit wavelength
                * `False` orresponds to energy flux per unit frequency

        """
        self.wavelength = wavelength
        self.spectrum = spectrum
        self.wlunit = wlunit
        self.flambda = flambda

    @classmethod
    def from_basel(cls, filename):
        """
        Construct SED object from BaSeL spectral library

        Parameters
        ----------
        filename: str
            Path to a BaSeL SED FITS table

        """
        with pyfits.open(filename) as hdu:
            wl = hdu[1].data
            sed = hdu[2].data

        return cls(wl, sed, wlunit=u.angstrom, flambda=False)

    @classmethod
    def from_ascii(cls, filename, cols=[0, 1], wlunit=u.angstrom):
        """
        Construct SED object from ASCII table

        Parameters
        ----------
        filename: str
            Path to ASCII SED file
        cols: list of int, optional
            Columns to use in the file (default is 0 and 1)

        """
        wl, sed = np.loadtxt(filename, usecols=cols, unpack=True)

        return cls(wl, sed, wlunit=wlunit, flambda=True)

    @property
    def photon_count(self):
        """
        Conversion from spectral energy distribution to photon number

        A conversion to photon number per unit wavelength is needed.
        For a single photon, Energy = h*f   where h is Planck's constant
        and f is frequency, or E = h*c/lambda where c is speed of light
        and lambda is wavelength.

        If the spectra is given in "f-lambda", which means energy flux
        per unit wavelength interval, one needs to divide by the energy
        of each photon, which means multiplying the f-lambda values by
        lambda.

        Otherwise it is usually given in "f-nu", which means energy flux
        per unit frequency interval, one needs to multiply by the energy
        of each photon, which means dividing the f-nu values by lambda.

        """
        if not self.flambda:
            return self.spectrum / self.wavelength

        return self.spectrum * self.wavelength

    def __call__(self, wavelength_range, ext_wlunit):
        """
        Return the photon count values for a wavelength range.

        Computes an interpolation over the available SED data to
        provide the photon count on a particular wavelength range.

        Parameters
        ----------
        wavelength_range: array_like
            Array of wavelength values to provide the SED for
        ext_wlunit: `astropy.units`
            Unit of the given wavelength range to perform interpolation

        Returns
        -------
        photon_count: array_like
            Photon count values for the given wavelengths

        """
        conversion_rate = self.wlunit.to(ext_wlunit)

        return np.interp(wavelength_range,
                         self.wavelength * conversion_rate,
                         self.photon_count)


class FlatSED(object):
    """Flat energy distribution"""

    def __call__(self, wavelength_range, wlunit=None):
        """
        The photon count is proportional to the wavelength for a flat sed
        """
        return wavelength_range


class TelescopeThroughput(object):
    """Telescope throughput class"""
    def __init__(self, wavelength, throughput, wlunit):
        """
        Parameters
        ----------
        wavelength: array_like
            Array of wavelength values corresponding to the SED
        throughput: array_like
            Array of telescope throughput
        wlunit: `astropy.units`
            Unit of the given wavelength

        """
        self.wavelength = wavelength
        self.throughput = throughput
        self.wlunit = wlunit

    @classmethod
    def from_file(cls, filename, wlunit, cols=[0, 1]):
        """
        Constructor for the throughput class from ASCII file

        Parameters
        ----------
        filename: str
            Path to throughput file
        wlunit: `astropy.units`
            Unit of wavelength
        cols: list of int, optional
            Columns to read in file (default is 0 and 1)

        """
        wlgth, thruput = np.loadtxt(filename, usecols=cols, unpack=True)

        return cls(wlgth, thruput, wlunit)

    def __call__(self, wavelength_range, ext_wlunit):
        """
        Return the telescope throughput values for a wavelength range.

        Computes an interpolation over the available data to provide
        the telescope throughput  on a particular range of wavelengths.

        Parameters
        ----------
        wavelength_range: array_like
            Array of wavelength values to provide the throughput for
        ext_wlunit: `astropy.unit`
            Unit of the given wavelength range to perform interpolation

        Returns
        -------
        throughput: array_like
            Telescope throughput values for the given wavelengths

        """
        conversion_rate = self.wlunit.to(ext_wlunit)

        return np.interp(wavelength_range,
                         self.wavelength * conversion_rate,
                         self.throughput)


class PSFcube(object):
    "Class that defines a wavelength array of monochromatic PSFs"
    def __init__(self, wavelength, data, pixel_size, header=None, wlunit=u.um):
        """
        Parameters
        ----------
        wavelength: array_like
            Wavelength range corresponding to the monochromatic PSFs
        data: array_like
            Stack of monochromatic PSFs
        pixel_size: float
            Pixel size of the data in microns
        header: FITS header, optional
            Header associated to FITS file
        wlunit: `astropy.unit`, optional
            Unit of the wavelength array (default um)

        """
        self.wavelength = wavelength
        self.psf_cube = data
        self.pixel_size = pixel_size
        self.header = header
        self.wlunit = wlunit

        self._psf_computed = False

    @classmethod
    def from_fits(cls, filename, shape, pixel_size, wlunit):
        """
        Constructor for the PSFcube from a FITS image file

        Parameters
        ----------
        filename: str
            Path to FITS image file of the PSFs
        shape: tuple of int
            Shape of the provided cube as (n_wavelength, psf_size, psf_size)
        pixel_size: float
            Pixel size of the data in microns
        wlunit: `astropy.unit`
            Unit of the wavelength array

        """
        data = np.empty(shape, dtype=float)
        n_psf = len(data)
        wavelength = np.empty(n_psf, dtype=float)

        for i in range(n_psf):
            img, hdr = pyfits.getdata(filename, ext=i+1, header=True)
            data[i] = img
            wavelength[i] = float(hdr['WLGTH0'])

        return cls(wavelength, data, pixel_size, header=hdr, wlunit=wlunit)

    @property
    def psf(self):
        """Accessor for the broadband PSF"""
        if not self._psf_computed:
            print("Broadband PSF not computed yet")
            return None

        return self._psf

    def expand_psf_cube(self):
        """
        Interpolate the wavefront pixels to expand the data along lambda
        """
        rtf_cube = npf.rfftn(npf.ifftshift(self.psf_cube), axes=(-2, -1))

        _, psf_size, rtf_size = rtf_cube.shape

        expanded_rtf_cube = np.empty((self.expanded_range.size,
                                      psf_size, rtf_size),
                                     dtype=complex)

        for x in range(rtf_size):
            for y in range(psf_size):
                # Extract the real and imaginary functions at the current point
                expanded_pixel = self.expand_carray(rtf_cube[:, y, x],
                                                    self.wavelength,
                                                    self.expanded_range)
                # Store the expanded (x, y) pixel point acorss the full cube.
                expanded_rtf_cube[:, y, x] = expanded_pixel

        return expanded_rtf_cube

    def to_broadband(self, sed, thruput, wl_step=1,
                     niemi=False, aocs=False, prf=False):
        """
        Compute the broadband PSF from the cube of monochromatic PSFs

        Parameters
        ----------
        sed: `SED`
            Instance of a stellar SED
        thruput: `TelescopeThroughput`
            Instance of a telescope throughput
        wl_step: float or int, optional
            Linear step between wavelength in expanded range
            (default is 1)
        niemi: bool, optional
            Apply the Niemi effect (default `False`)
        aocs: bool, optional
            Apply the AOCS effect (default `False`)
        prf: bool, optional
            Apply the  pixel response effect (default `False`)

        Returns
        -------
        broadband_psf: array_like
            Image of the broadband PSF

        """
        # Create the range of wavelengths for interpolation
        self.expanded_range = np.arange(self.wavelength[0],
                                        self.wavelength[-1] + 1,
                                        step=wl_step)

        # Retrieve the expanded wavefront
        expanded_wf = self.expand_psf_cube()

        # Apply SED + Telescope throughput weighting
        total_thruput = (sed(self.expanded_range, ext_wlunit=self.wlunit) *
                         thruput(self.expanded_range, ext_wlunit=self.wlunit))
        expanded_wf *= total_thruput[:, None, None]

        # Apply the Niemi Effect (WL Dependant)
        if niemi:
            expanded_wf *= self.niemi_effect(expanded_wf.shape,
                                             self.expanded_range,
                                             self.pixel_size)

        # Integrate of wavelength range
        broadband_wf = np.sum(expanded_wf, axis=0)

        # Apply AOCS Effect to the coadded PSF (Not WL Dependant)
        if aocs:
            broadband_wf *= self.aocs_effect()

        # Apply SubPixel Response Effect to the coadded PSF (Not WL Dependant)
        if prf:
            broadband_wf *= self.pixelresponse_effect()

        broadband_psf = np.abs(npf.ifftshift(npf.irfft2(broadband_wf)))

        self._psf = self.normalize(broadband_psf)
        self._psf_computed = True

    def downsample_detector(self, ds_factor=1):
        """
        Downsampled the broadband PSF array to detector resolution.

        Parameters
        ----------
        ds_factor: float, optional
            Downsampling factor to apply to the broadband PSF
            (default is 1)

        """
        if ds_factor > 1:
            psf_shape = self.psf.shape

            y, x = np.indices(psf_shape) // ds_factor
            new_shape = np.array(psf_shape, dtype=int) // ds_factor

            regions = x + y * psf_shape[1] // ds_factor

            psf_sum = spn.sum(self.psf,
                              labels=regions,
                              index=np.arange(regions.max() + 1))
            downsampled_psf = psf_sum.reshape(*new_shape)

            self._psf = self.normalize(downsampled_psf)

    def write_effective_psf(self, outfits, clobber=False):
        """
        Write broadband PSF to disk

        Parameters
        ----------
        outfits: str
            Output FITS file
        clobber: bool, optional
            Overwrite file if already existing (default `False`)

        """
        self._update_header_info()

        header = copy.copy(self.header)
        comment = ("========================================="
                   "Broadband PSF created from {} wavelengths"
                   "=========================================")
        header.add_comment(comment.format(len(self.wavelength)))
        try:
            pyfits.writeto(outfits, data=self.psf,
                           header=header, clobber=clobber)
            print("Output image written in {}".format(outfits))
        except:
            print("Output file {} already existing".format(outfits))

    def _update_header_info(self):
        """
        Ensure the header has the minimum required entries
        """
        if self.header is None:
            hdu = pyfits.PrimaryHDU(data=self.psf)
            self.header = hdu.header

    @staticmethod
    def normalize(data):
        """
        Normalize to unity

        Parameters
        ----------
        data: array_like
            Array data to normalize

        Returns
        -------
        normalized_data: array_like
            Normalized array to unity

        """
        return data / data.sum()

    @staticmethod
    def niemi_effect(data_shape, wavelength, pixel_size,
                     niemi_params=(17.45, -0.237, 45.10, -0.378)):
        """
        Weight function for the wavelength-dependent detector PRF

        Parameters
        ----------
        data_shape: tuple of int
            shape of the RTF of the psf cube
        wavelength: array_like
            wavelength range
        pixel_size: int or float
            output pixel size
        niemi_params: tuple of floats, optional
            Wavelength dependant measured parameters along x and y axes

        Returns
        -------
        niemi_weight: array_like
            Niemi effect to apply to PSF cube

        References
        ----------
        Niemi et al. 2015

        """
        # Niemi et al. wavelength dependence parameters.
        sx_const, sx_power, sy_const, sy_power = niemi_params

        _, psf_size, rtf_size = data_shape

        # Construct wavelength axis
        wavelength = wavelength[:, None, None]

        # Calculate gaussian pixel covariance with scaling factors
        # Equations 8 - 10
        sigma_x = 2 * np.pi * sx_const * wavelength ** sx_power / psf_size
        sigma_y = 2 * np.pi * sy_const * wavelength ** sy_power / psf_size

        # 1 / um^2
        sigma_x_sq = (sigma_x / pixel_size) ** 2
        sigma_y_sq = (sigma_y / pixel_size) ** 2

        # Gaussian weight function.
        y, x = np.indices((psf_size, rtf_size))

        # Separate first and second halves of y values
        y = np.where(y >= psf_size // 2, y - psf_size, y)
        arg = x ** 2 * sigma_x_sq + y ** 2 * sigma_y_sq
        niemi_weight = np.exp(-arg / 2)

        return niemi_weight

    @staticmethod
    def aocs_effect(data_shape, pixel_size,
                    aocs_sigma=0.025, flength=24.75):
        """
        Weight function for the AOCS contribution

        Calculate the gaussian weight function caused by the telescope's
        distorted AOCS errors, and multiply in Fourier space.
        Assumes the sky AOCS gaussian is symmetric with rms size aocs_sigma.

        Parameters
        ----------
        data_shape: tuple of int
            Shape of the broadband data
        pixel_size: float
            Pixel size in microns
        aocs_sigma: float, optional
            AOCS guiding error rms in arcsec
        flength: float, optional
            Default telesope focal length in meters

        Returns
        -------
        aocs_weight: array_like
            AOCS Response to apply to the wavefront

        Notes
        -----
        The telescope distortion is supplied in 2x2 matrix dmatrix.
        This is currently calculated by affineDFT and could be passed into
        this function, but dmatrix will default to an undistorted telescope
        model with the nominal focal length if no distortion matrix is set.

        References
        ----------
        Equations are given in section 2.3, especially equation 7.

        """
        # Height of the output height and width.
        psf_size, rtf_size = data_shape

        # Construct the gaussian covariance matrix for the sky AOCS errors
        # units of arcsec squared in real data this would be provided externally
        # from the AOCS system assume here it's a symmetric diagonal covariance
        # (i.e. a circular smoothing function)

        aocs_covariance_matrix = np.matrix([[aocs_sigma ** 2, 0],
                                            [0, aocs_sigma ** 2]])

        # Now we have to apply the distortion matrix A as in equation 7 of
        # section 2.3. If the distortion matrix has not been provided
        # externally, calculate it assuming no distortion check value of
        # determinant of dmatrix and set no-distortion default values if zero
        # or negative.
        flength_pix = flength * 1e6 / pixel_size
        flength_arcsec = flength_pix * np.pi / 180 / 3600

        dmatrix = np.matrix([[1 / flength_arcsec, 0],
                             [0, 1 / flength_arcsec]])

        # Assuming the input distortion matrix has been supplied as A^T,
        # which is what is calculated inside affineDFT,
        # calculate its inverse to get {A^T}^{-1}
        # (which is the same as {A^{-1}}^T)
        mat_d_matrix = dmatrix.I

        # Then the product A^{-1} C {A^T}^{-1}
        # First transpose {A^T}^{-1} to get A^{-1}
        mat_d_matrix_t = mat_d_matrix.T

        # Multiply by the distortion matrix as in Section 2.3. First create
        # the product C_{G,theta} times the inverse of the distortion matrix.
        # Then multiply to get a new 2x2 matrix, C' as in equation 7.
        aocs_pixcovar = mat_d_matrix_t * aocs_covariance_matrix * mat_d_matrix

        # Apply additional scaling factors from equation 7.
        aocs_pixcovar *= (2 * np.pi / psf_size) ** 2

        # Make sure we deal with the case x or y = 0, when sinc=1.
        y, x = np.indices((psf_size, rtf_size), dtype=float)

        # Separate first and second halves of y values
        y = np.where(y >= psf_size // 2, y - psf_size, y)

        # Matrix product
        yx = np.indices(data_shape).reshape(2, psf_size * rtf_size)
        yval, xval = aocs_pixcovar * yx

        arg = x * xval.reshape(data_shape) + y * yval.reshape(data_shape)
        aocs_weight = np.exp(-arg / 2)

        return aocs_weight

    @staticmethod
    def pixelresponse_effect(data_shape):
        """
        Weight function for the sub-pixel (oversampled) tophat response.

        Calculate the sinc weight function appropriate for the tophat
        oversampled subpixel response.

        Parameters
        ----------
        data_shape: tuple of int
            Shape of the broadband data

        Returns
        -------
        prf_weight: array_like
            Sub-Pixel Response to apply to the wavefront

        References
        ----------
        Equations are given in section 2.3, equation 5 does not depend
        on wavelength does not depend on pixel oversampling.

        """
        psf_size, rtf_size = data_shape

        # Multiply complex array by sinc weight function.
        y, x = np.indices((psf_size, rtf_size), dtype=float)

        # Separate first and second halves of y values
        y = np.where(y >= psf_size // 2, y - psf_size, y)

        xarg = np.pi * x / psf_size
        yarg = np.pi * y / psf_size

        prf_weight = np.sinc(xarg) * np.sinc(yarg)

        return prf_weight

    @staticmethod
    def expand_carray(data, short_range, expanded_range):
        """
        Interpolation on complex data array

        Parameters
        ----------
        data: `ndarray` of complex
            Complex array to expand through interpolation
        short_range: array_like
            Abscissa corresponding to the data array
        expanded_range: array_like
            Abscissa of the interpolated array

        Returns
        -------
        interpolated_data: `ndarray` of complex
            Expanded data array

        """
        real_spl = spi.InterpolatedUnivariateSpline(short_range, data.real)
        imag_spl = spi.InterpolatedUnivariateSpline(short_range, data.imag)

        interpolated_data = (real_spl(expanded_range) +
                             imag_spl(expanded_range) * 1j)

        return interpolated_data


def local_test():
    dpath = '/Users/aboucaud/work/IAS/stage/tperdereau/'
    fitsfile = os.path.join(dpath, 'psf_test_euclid.fits')
    transfile = os.path.join(dpath, 'broadband/throughput/VISThroughBE.dat')
    outfile = os.path.join(dpath, 'psf_broadband_test.fits')

    psf = PSFcube.from_fits(fitsfile, shape=(11, 512, 512),
                            pixel_size=2, wlunit=u.um)
    sed = FlatSED()
    thru = TelescopeThroughput.from_file(transfile, cols=[0, 4], wlunit=u.nm)
    psf.to_broadband(sed=sed, thruput=thru, wl_step=10)
    psf.write_effective_psf(outfile)
    psf.downsample_detector(4)
    psf.write_effective_psf(outfile.replace('.fits', '_dwspl.fits'))


if __name__ == '__main__':
    local_test()
