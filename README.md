# SyncSSCModel.jl
Radiative models for simulating emission spectra in various black hole jet models.

You can find the [documentation on GitHub](https://astro-group-bristol.github.io/SyncSSCModel.jl/dev/).

## Load the module
Installing the Spectral Fitting package version 0.3.11.
```julia
julia>] registry add https://github.com/astro-group-bristol/AstroRegistry
julia>] add SpectralFitting@0.3.11
```

Installing the SyncSSCModel package
```julia
julia>] add SyncSSCModel
```
or, you can also check out this repository in development using the following command lines
```julia
julia>] dev --local https://github.com/astro-group-bristol/SyncSSCModel.jl 
julia>] instantiate
```

## Simple fit example
Find the example scripts in `examples`. Two files are required for each source. 
* The first file contains the Julia script to do the fitting (e.g. `PicA.jl`). You can enter parameter values according to the source and define which parameters to freeze and which to fit.
* The second file is a `.txt` file (e.g. `PicA.txt`) where you put all your observation data sets in a logarithmic scale. The first column is for the frequency ``log_{10} \\nu`` in Hz, and the second column is for the spectral power flux ``log_{10} \\nu F_{\\nu}`` in cgs units ``erg cm^{-2} s^{-1}``.
