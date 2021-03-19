## ISETImagePipeline
This repository contains all the analysis code used in our paper: LQ Zhang, NP Cottaris, and DH Brainard (2021). An image reconstruction framework for characterizing early vision.  

For a 10 minutes introduction of our project, here is our [V-VSS 2020 Talk](https://youtu.be/d5qI0FNCAv4).  

Feel free to contact me if you any questions or comments!

## Dependencies
This project is written in MATLAB, and is built upon two other repositories: 
- [ISETBio](https://github.com/isetbio/isetbio), which is an open-source, accurate computational model for the early visual system
- [Reconstruction Toolbox](https://github.com/isetbio/ISETPipelineToolbox), which is a set of routine for our Bayesian image reconstruction algorithm.
- The best way to set up all the dependencies is to use [ToolboxToolbox](https://github.com/ToolboxHub/ToolboxToolbox), which is a MATLAB package manager developed by other wonderful people in our lab. Once you set everything up, you can clone (or download) this repo to `Matlab-User-Path/projects/ISETImagePipeline`. Type `tbUseProject(ISETImagePipeline)` in MATLAB, and `ToolboxToolbox` will set everything up for you!

- Alternatively, you can also manually download the other two dependencies, just make sure everything is on MATLAB's search path before running the code in here!

## Getting started
We have written a set of tutorial/live script that are designed to demonstrate the basic usage of our code. You can find them under `recon_2020`. We suggust you go through them in the order as listed.

```
...
└── recon_2020
    ├── constructPrior.m    Construct sparse-coding based prior of natural color images
    ├── imageRecon.mlx      Basic routine for image reconstruction from cone excitation
    ├── dichromacy.mlx      Image reconstruction from a dichromatic retinal mosaic
    ├── anomTrichroma.mlx   Image reconstruction from an anomalous trichromacy retinal mosaic
    ├── aliasing.mlx        Simulate the experiment of Williams, 1985: Aliasing in human foveal vision
    ├── lmConeRatio.m       Manipulate the L/M cone ratio and how it influences image reconstruction
    ├── sConeRatio.m        Manipulate the S cone ratio, chromatic aberration, and lens/pigment density 
...
```