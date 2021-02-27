# SAIMFitKit.jl
Julia implementation of the Scanning Angle Interference Microscopy (SAIM) fitting algorithms

## Description
Fitting algorithm reconstructs sample height from a stack of scanning angle interference images. More details on the fiting algorithms can be found in the following papers:

Colville MJ, Park S, Zipfel WR, and Paszek MJ. High-speed device synchronization in optical microscopy with an open-source hardware control platform. Scientific Reports. 9: 1-13. 

Paszek MJ, Dufort C, Rubashkin MG, Davidson MW, Thorn KS, Liphardt JT, and Weaver VM.  2012.  Scanning angle interference microscopy reveals cell dynamics at the nano-scale.  Nature Methods 9: 825-7

Existing Functionality
----------------------

`fit = fit_SAIM(file_path, file_name, optic, angles, init_params, lower_bounds, upper_bounds; kwargs...)`:

* `file_path`: string that provides the path to the image file
* `file_name`: string that provides the name of the image file (.tif extension should not be included in name)`
* `optic`: composite type of imaging parameters ('SAIMOptic')
* `angles`: the incidence angles in degrees for the acquisition sequence 
* `lower_bounds`: lower bounds for the fit parameters (e.g. A, B, H)
* `upper_bounds`: lower bounds for the fit parameters (e.g. A, B, H)
* `kwargs`: addition parameters for fitting, such as `glb` (run global grid search when true), 'step' (step size for grid search), `show` (display plot of fit heights when true), 'color' (number of excitation wavelengths; typically one color)
* `fit`: composite type of results (`SAIMFitResult`)

Fit SAIM data using non-linear least squares (LsqFit.jl) using either a local or global grid search.

----
