# Broadband PSF suite

Convert a set of monochromatic PSFs to an effective PSF for a given instrument and stellar SED.

## Preparation

The code runs on single or multiple files.

* The monochromatic PSFs need to be saved as extensions of a single image FITS file.
* The stellar SEDs can be stored as an ASCII file or a FITS table.
* The configuration file must be properly filled with information on the input data

```
[psf]
size = 512          # PSF image size (int)
pixel_size = 2      # PSF pixel size in microns (float)
n_wavelength = 11   # Number of monochromatic PSF in the FITS file (int)
wlunit = um         # Wavelength unit in the FITS header given in astropy.units

[thruput]
filename = /Users/aboucaud/work/Euclid/data/filters/VISThroughBE.dat
wl_col = 0          # Index of the wl column (int)
thruput_col = 4     # Index of the throughput column (int)
wlunit = nm         # Wavelength unit of thruput file given in astropy.units

[sed]
type = ascii        # Type of sed input file (ascii|fits)
wlunit = angstrom   # Wavelength unit of SED file given in astropy.units

[broadband]
wl_step = 10        # Wavelength step for the interpolation (int|float)
use_niemi = False   # Apply Niemi wavelength-dependent weighting (bool)
use_aocs = False    # Apply AOCS response weighting (bool)
use_prf = False     # Apply point response weighting (bool)
down_fact = 2       # Downsampling factor (int >= 1)

[output]
path = .            # Relative path to output folder
nameformat = {psf}_broadband_{sed}.fits
```

## Usage

### Single file example

```bash
python mk_broadband.py monochromatic_psf_cube.fits stellar_sed.fits -c config_file.txt
```

### Multiple files example

```bash
python mk_broadband.py psf*.fits *sed.txt -c config_file.txt
```

This will loop over the input PSF files **and** the various SEDs.

### Multiprocessing option

```bash
python mk_broadband.py psf*.fits *sed.txt -c config_file.txt -n 4
```

This will run the previous command but use 4 processors in parallel instead.
