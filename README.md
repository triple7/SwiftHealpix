# swift-healpix

swift-healpix is the swift wrapper for the HEALPIX framework created by NASA for efficient spherical to planar conversions of spacecraft data captures. A brief primer is available [here](https://healpix.jpl.nasa.gov/pdf/intro.pdf)

## Technical overview

In brief, Healpix allows for the subdivision to arbitrarily high resolution of spherical discreet functions, originally used for spherical harmonics in cosmology, but has an explicit use the hierarchical equal area iso latitudinal pixellation of the sphere into arbitrarily high degrees of subdivision. 

The HIPS format [(paper)](https://www.ivoa.net/documents/HiPS/20170406/PR-HIPS-1.0-20170406.pdf) uses Healpix for sky viewers such as [Aladin](https://aladin.u-strasbg.fr)
