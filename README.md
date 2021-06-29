# Diffusion-Relaxation Correlation Imaging Data Reconstruction

This GitHub repository contains MATLAB code for reconstructing diffusion-relaxation correlation data using a Monte Carlo inversion algorithm. Part of this toolbox is the complete md-dmri toolbox which also covers additional topics such as alternative reconstruction methods, preprocessing methods, and diffusion gradient optimization. An exhaustive readme covering this toolbox may for example be found at https://github.com/daniel-topgaard/md-dmri.

The data and a corresponding mask need to be saved as nii.gz files. The data on how the experiment was performed is supplied as an xps.mat file ('eXperiment Parameter Structure'). More details on how the xps is put together may be found at https://github.com/daniel-topgaard/md-dmri.

Example Usage:

```matlab
run_inversion_pipeline(fn_data, fn_mask, fn_xps);
```

Individual bootstraps may be collected in a single data struct using jm_bs2m. dtr2r1d parameters may then be calculated using the function jm_m2pars.