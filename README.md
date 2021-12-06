## Bayesian Image Reconstruction from Cone Mosaic Signal
This repository contains all the analysis code used in our paper: **LQ Zhang, NP Cottaris, and DH Brainard (2021). An Image Reconstruction Framework for Characterizing Early Vision.** 

For a 10-minute introduction to our project, please see our [V-VSS 2020 Talk](https://youtu.be/d5qI0FNCAv4).  

## Dependencies
This project is written in MATLAB, and is built upon three other repositories: 
- [ISETBio](https://github.com/isetbio/isetbio/wiki), which is an open-source, accurate computational model for the early visual system.
- [Reconstruction Toolbox](https://github.com/isetbio/ISETPipelineToolbox), which is a set of routine for our Bayesian image reconstruction algorithm.
- [CSFGenerator](https://github.com/isetbio/ISETBioCSFGenerator), which is a generic and flexible codebase for simulating Contrast Sensitivity Function (CSF).

## Getting Started
- The best way to set up all the dependencies is to use [ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox), which is a MATLAB package manager developed by other wonderful people in our lab. Once you set everything up, you can `clone` (or download) this repo to `Matlab-User-Path/projects/ISETImagePipeline`. Then, simply type `tbUseProject('ISETImagePipeline')` in MATLAB, and `ToolboxToolbox` will set everything up for you!

- Alternatively, you can also manually download the other three dependencies, just make sure everything is on MATLAB's search path before running the code in here!

We have written a set of tutorial/live scripts that are designed to demonstrate the basic usage of our code. You can find them under `recon_2020`. We have tried to comment these well, and they provide detail on how to reproduce the analysis reported in the paper. We suggust you go through them in the order listed in the table below \[1\].

```
...
├── display.mat             # Parameters for the display we are using in the simulation
├── sparsePrior.mat         # A pre-trained sparse coding prior (the basis function)
└── recon_2020
    ├── imageRecon.mlx      # Basic routine for image reconstruction from cone excitation [2]
    ├── dichromacy.mlx      # Image reconstruction from a dichromatic retinal mosaic [2]
    ├── anomTrichroma.mlx   # Image reconstruction from an anomalous trichromacy retinal mosaic [2]
    ├── lmConeRatio.m       # Manipulate the L/M cone ratio and how it influences image reconstruction [3]
    ├── sConeRatio.m        # Manipulate the S cone ratio, chromatic aberration, and lens/pigment density [3]
    ├── priorEffect.m       # Effect of prior on optimal mosaic design (cone ratio) [3]
    ├── reconPeripheral.m   # Image reconstruction at different visual eccentricity 
    ├── aliasing.mlx        # Simulate the experiment of Williams, 1985: Aliasing in human foveal vision
    ├── constructPrior.m    # Construct sparse-coding based prior of natural color images [4]
...
```

\[1\] Our project is organized in a way that if you are interested in the details of our simulation and reconstruction (i.e., you are reading the Methods section), you should also take a look at [Reconstruction Toolbox](https://github.com/isetbio/ISETPipelineToolbox), which is the actual "library" for most of our analysis.  

\[2\] if you are mainly interested in quickly playing with our code, remember to use a **small** parameter for everything so it runs fast (e.g., 0.25 visual degree mosaic, 32 by 32 image size).  

\[3\] These analysis can be computationally intensive (i.e., couple of days on a laptop). We recommend running these on a server, and take advantage of MATLAB's parallel pool functionality. 

\[4\] Running the prior learning routine also requires a large dataset of natural images. We used [ILSVRC](http://www.image-net.org/challenges/LSVRC) in our analysis. We have stored a learned prior `sparsePrior.mat` in case you don't want to repeat this step.  

We have shared both the RGB and hyperspetral image dataset, the parameters used in the simulation including display and cone mosaic setup, as well as some of the intermediate results such as the learned sparse prior, likelihood function (i.e., render matrix), these are available through here: [https://tinyurl.com/26r92c8y](https://tinyurl.com/26r92c8y)  
 

## Others Code & Script
There are many other scripts and functions in this repo, most of them are either things we used to run our analysis at scale, or other things that we have done but didn't go into this paper (e.g., regression-based reconstruction method). Here are some of them you might find useful:

```
...
├── compute
    ├── reconGPU.m          # Yes, if you have a NVIDIA GPU, you should run the reconstruction on it
    ├── ...                 # It is much faster since most of the computation is matrix related
├── csf                     # Code related to the simulation of contrast sensitivity function (CSF)
    ├── spatialCSF_Cone.m   # CSF simulation for Poisson ideal observer based on cone excitation 
    ├── spatialCSF_Recon.m  # CSF simulation for our image reconstruction based observer
    ├── ...
├── hyperspectral
    ├── ratioHyperspec.m    # Image reconstruction and cone ratio analysis with hyperspectral images
...
```

## Contact
Feel free to contact me if you any questions or comments at   
**lingqiz at sas dot upenn dot edu**
