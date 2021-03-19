## ISETImagePipeline
This repository contains all the analysis code used in our paper: LQ Zhang, NP Cottaris, and DH Brainard (2021). An Image Reconstruction Framework for Characterizing Early Vision.  

For a 10 minutes introduction of our project, here is our [V-VSS 2020 Talk](https://youtu.be/d5qI0FNCAv4).  

## Dependencies
This project is written in MATLAB, and is built upon three other repositories: 
- [ISETBio](https://github.com/isetbio/isetbio), which is an open-source, accurate computational model for the early visual system.
- [Reconstruction Toolbox](https://github.com/isetbio/ISETPipelineToolbox), which is a set of routine for our Bayesian image reconstruction algorithm.
- [CSFGenerator](https://github.com/isetbio/ISETBioCSFGenerator), which is a generic and flexible codebase for simulating Contrast Sensitivity Function (CSF).

## Getting started
- The best way to set up all the dependencies is to use [ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox), which is a MATLAB package manager developed by other wonderful people in our lab. Once you set everything up, you can clone (or download) this repo to `Matlab-User-Path/projects/ISETImagePipeline`. Simply type `tbUseProject('ISETImagePipeline')` in MATLAB, and `ToolboxToolbox` will set everything up for you!

- Alternatively, you can also manually download the other three dependencies, just make sure everything is on MATLAB's search path before running the code in here!

We have written a set of tutorial/live scripts that are designed to demonstrate the basic usage of our code. You can find them under `recon_2020`. These code are well commented and contains (hopefully) the details on how to reproduce the analysis reported in the paper. We suggust you go through them in the order as we listed.

```
...
└── recon_2020
    ├── constructPrior.m    # Construct sparse-coding based prior of natural color images
    ├── imageRecon.mlx      # Basic routine for image reconstruction from cone excitation
    ├── dichromacy.mlx      # Image reconstruction from a dichromatic retinal mosaic
    ├── anomTrichroma.mlx   # Image reconstruction from an anomalous trichromacy retinal mosaic
    ├── aliasing.mlx        # Simulate the experiment of Williams, 1985: Aliasing in human foveal vision
    ├── lmConeRatio.m       # Manipulate the L/M cone ratio and how it influences image reconstruction
    ├── sConeRatio.m        # Manipulate the S cone ratio, chromatic aberration, and lens/pigment density 
...
```

If you are interested in the details of our simulation and reconstruction (i.e., you are reading the Methods section), you should also take a look at [ISETBio](https://github.com/isetbio/isetbio) and [Reconstruction Toolbox](https://github.com/isetbio/ISETPipelineToolbox).

## Additional 
There are many other scripts and functions in this repo, most of them are either things we used to run our analysis at scale, or other things that we have done but didn't go into this paper. Here are some of them you might find useful:

```
...
├── compute
│   ├── reconGPU.m          # Yes, if you have a NVIDIA GPU, you should run the reconstruction on it
│   ├── ...                 # It is much faster since most of the computation is matrix related
├── csf                     # Code related to the simulation of contrast sensitivity function (CSF)
│   ├── spatialCSF_Cone.m   # CSF simulation for Poisson ideal observer based on cone excitation 
│   ├── spatialCSF_Recon.m  # CSF simulation for our image reconstruction based observer
│   ├── ...
├── dichromacy
│   ├── MarkovPrior.m       # Gaussian prior for which we have control over its spatial and chromatic correlation
│   ├── priorEffect.m       # The effect of prior on optimal mosaic design (cone ratio)  
│   ├── labDistance.m       # The Spatial CIELAB loss function
│   ├── ...
├── peripheral
│   ├── computePeripheral.m # Image reconstruction at different visual eccentricity
│   ├── ...
...
```

## Contact
Feel free to contact me if you any questions or comments at   
**lingqiz at sas dot upenn dot edu**
