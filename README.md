![icon](assets/logo.png)


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gatocor.github.io/FlowCytometry.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gatocor.github.io/FlowCytometry.jl/dev/)
# FlowCytometry.jl

The package wants to offer a complete set of tools for the analysis of cytometry dataset in Julia. 

The package can be used to work with any kind of cytometry data:

 - Flow cytometry
 - Spectral cytometry
 - Mass cytometry

The package provides with functionalities for:

 - Loading .ftc files.
 - Data structures for the manipulation and analysis of cytometry data.
 - Gating in a manual or automatic manner.
 - Quality control analysis.
 - Compensation tools.
 - Dimensionality reduction.
 - Cluster discovery.
 - Comparative analysis for perturbation discovery (healthy-ill,...).
 - Loading and saving of processed data in HDF5 format, compatible to be read in many other platforms.

## Installation

Until we register the package, the package can be installed from github as:

```
using Pkg; pkg.add("https://github.com/gatocor/FlowCytometry.jl")
```

