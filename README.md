# swift-healpix

swift-healpix is the swift wrapper for the HEALPIX framework created by NASA for efficient spherical to planar conversions of spacecraft data captures. A brief primer is available [here](https://healpix.jpl.nasa.gov/pdf/intro.pdf)

## Technical overview

In brief, Healpix allows for the subdivision to arbitrarily high resolution of spherical discreet functions, originally used for spherical harmonics in cosmology, but has an explicit use the hierarchical equal area iso latitudinal pixellation of the sphere into arbitrarily high degrees of subdivision. 

The HIPS format [(paper)](https://www.ivoa.net/documents/HiPS/20170406/PR-HIPS-1.0-20170406.pdf) uses Healpix for sky viewers such as [Aladin](https://aladin.u-strasbg.fr)


## Understanding the library

First, have a read at [Gorski (2005)](http://iopscience.iop.org/article/10.1086/427976/pdf)

The library includes all invertable mappings from each coordinate system as a set of class functions you can call for your use case.

There are a variety of spherical representations in this package:
1. alpha, delta)
2. (theta, phi)
3.  (x, y, z)

They are all normalised here to (z, a) 

The HEALPIX spherical projection is used to map to (t, u). These are the functions:
. za2Tu
. tu2Za

See Section 4.4 and Figure 5 in the paper, where `(t, u)` is called `x_s, y_s).

A simple affine transformation is used to (f, x, y). These functions are:
. tu2Fxy
. fxy2Tu

where `f = {0 .. 11}` is the base pixel index and `(x, y)` is the position within the base pixel in the (north-east, north-west) direction
and `(0, 0)` in the south corner.

From `(f, x, y)`, the HEALPix pixel index in the "nested" scheme is related via the functions:
. fxy2nest
. nest2fxy

And in the "ring" scheme via the functions:
. fxy2ring
. ring2fxy

## Geometrical transformations

There are 2 transformations:

. (z, a) <-> (t, u) is the HEALPix spherical projection. 
. (t, u) <-> (f, x, y) is a 45 deg rotation and scaling for each of the 12 base pixels, so that HEALPix pixels in (x, y) are unit squares

Pixel index computations are relatively straightforward both in the "nested" and "ring" pixelisation schemes.

## Notation

. theta: colatitude (pi/2 - delta) [0 : pi]
. phi: longitude (alpha) [0 : 2pi]
. t: x axis coordinate in spherical projection [0 : 2pi]
. u: y axis coordinate in spherical projection [-pi/2 : pi/2]
. z: cos(theta) [-1 : 1]
. X: sin(theta) * cos(phi) [-1 :  1]
. Y: sin(theta) * sin(phi) [-1 : 1]
. a : phi [0 : 2pi]
. f: base pixel index [0 : 11]
. x: north-east index in base pixel [0 : nside]
. y: north-west index in base pixel [0, nside]
. p: north-east axis in base pixel [0, 1]
. q: north-west axis in base pixel [0, 1]
. j: pixel-in-ring index polar cap: [1 : 4 i] equatorial belt: [1 : 4 nside]
. i: ring index [1 : 4 nside - 1]

## Arbitrary subdivision 

This package is optimised to use integers larger than MaxInt for arbitrarily large integers for identifying tiles in a given hierarchical map. This was initially intended to provide higher levels of subdivision to open up extremely large images that can be subdivided to fine grain tiles. The effect is the ability to zoom from a point where the earth is a small dot in the sky, down to the molecule of a dust particle on an ant's head without overloading either time or space complexity.
