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

## Simple fit examples
Find the example scripts in [examples](examples). Two files are required for each source.
- The first file contains the Julia script to do the fitting (e.g. [PicA.jl](examples/PicA.jl)). You can enter parameter values according to the source and define which parameters to freeze and which to fit.
- The second file is a `.txt` file (e.g. [PicA.txt](examples/PicA.txt)) where you put all your observation data sets in a logarithmic scale. The first column is for the frequency $\log_{10} \nu$ in Hz, and the second column is for the spectral power flux $\log_{10} \nu F_{\nu}$ in cgs units $\text{erg cm}^{-2} \text{s}^{-1}$.
- There are example fits to Pictor A ([PicA.jl](examples/PicA.jl)) and PKS 0637-752 ([PKS0637-752.jl](examples/PKS0637-752.jl)).

## Pluto interactive examples
You can use [Pluto](https://plutojl.org/) to interactively run the examples.

```julia
import Pkg; Pkg.add("Pluto")
import Pluto; Pluto.run()
```

Then open the file [examples/interactivePicA.jl](examples/interactivePicA.jl) or [examples/interactivePKS.jl](examples/interactivePKS.jl) in Pluto, in your web browser.

*Note*: the examples requires the [PlutoUI.jl](https://github.com/JuliaPluto/PlutoUI.jl) module.