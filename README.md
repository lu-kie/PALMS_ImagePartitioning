# PALMS_ImagePartitioning
PALMS_ImagePartitioning is a Matlab/C++ toolbox for image partitioning 
using the piecewise affine-linear Mumford-Shah model (also known as "affine-linear Potts model").

## Partitioning of color images

   - Supports segmentation of vector-valued images (e.g. RGB images, multispectral images, feature images)
   - Linear complexity in number of channels
   - Label-free: no label discretization required
   
   ![titleimageA](/docs/titleimageA.png)
   ![titleimageB](/docs/titleimageB.png)
   
   Left: natural image; Right: partitioning using the piecewise affine-linear Mumford-Shah model
   
   - Avoids oversegmentation of images with linear trends (e.g. the sky in a landscape image, illumination gradients)
   
   ![PottsAndPALMS](/docs/piecewiseConstantAndPiecewiseAffineLinear.png)
   
   Left: natural image; Center: classical (piecewise constant) Potts model; Right: piecewise affine-linear Mumford-Shah model
   
   - Scale of the partitioning is controlled by a model parameter
   
   ![parameter](/docs/parameter.png)
   
   Large values of the model parameter give few segments, while small choices lead to more segments
and closeness to the data

### Online demo
An online demo can be found at http://www.ipol.im/pub/art/2020/295/

## Installation
### Compiling
The algorithm depends on a mex script that needs to be compiled before execution. For compilation inside MATLAB, cd into the 'src/cpp' folder and run build.m

Requires the Armadillo and OpenMP library

Tested with Armadillo 8.400 https://launchpad.net/ubuntu/+source/armadillo/1:8.400.0+dfsg-2 and OpenMP 4.0.
On Linux, just use your package manager to install it:

sudo apt-get install libarmadillo-dev

sudo apt-get install libomp-dev

### Running
For a test run on the test image "redMacaw", run demo.m. 
demo.m calls the main function, affineLinearPartitioning.m
Arguments of affineLinearPartitioning.m are:
 - f: input image (double)
 - gamma: boundary penalty (larger choice -> less segments)
 - varargin: optional input parameters

## References
- L. Kiefer, M. Storath, A. Weinmann.
    "An efficient algorithm for the piecewise affine-linear Mumford-Shah model based on a Taylor jet splitting."
    IEEE Transactions on Image Processing, 2020.
- L. Kiefer, M. Storath, A. Weinmann.
    "PALMS image partitioning – a new parallel algorithm for the piecewise affine-linear Mumford-Shah model."
     Image Processing On Line, 2020.
