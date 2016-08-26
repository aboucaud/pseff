# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).


## [Unreleased]

## [0.2.1] - 2016-08-26
### Added
- Installation section in README
- astropy.units to the imported modules inside pseff (as 'u')

### Fixed
- from_basel constructor uses P.Hudelot convenient FITS file now

## [0.2.0] - 2016-07-28
### Added
- Multiprocessing feature to mk_broadband. This is accessible through the
`-n` or `--nprocs` command-line option.
- Add `CHANGELOG`
- Add new example to `README`

## [0.1.1] - 2016-07-28
### Changed
- Switched to relative imports

### Fixed
- Fix all weight function using a comparison with the euclid-vis-psf code.
- Fix weight function call in `to_broadband` method

## [0.1.0] - 2016-07-27
### Changed
- Change name of code to PSefF
- Code refactored as an open source Python project

### Added
- Add documentation to all methods
- Add `README` file
- Add `setup.py`

## [0.0.1] - 2016-07-16
### Added
- Initial commit of the code. One single file.


[Unreleased]: https://git.ias.u-psud.fr/aboucaud/broadband_psf/compare/v0.2.1...HEAD
[0.2.1]: https://git.ias.u-psud.fr/aboucaud/broadband_psf/compare/v0.2.0...v0.2.1
[0.2.0]: https://git.ias.u-psud.fr/aboucaud/broadband_psf/compare/v0.1.0...v0.2.0
[0.1.0]: https://git.ias.u-psud.fr/aboucaud/broadband_psf/compare/
[0.0.1]: https://git.ias.u-psud.fr/aboucaud/broadband_psf/compare/
