# SAIMFitKit.jl
Julia implementation of the Scanning Angle Interference Microscopy (SAIM) fitting algorithms

## Description
Fitting algorithm reconstructs sample height from a stack of scanning angle interference images. More details on the fiting algorithms can be found in the following papers:

Colville MJ, Park S, Zipfel WR, and Paszek MJ. High-speed device synchronization in optical microscopy with an open-source hardware control platform. Scientific Reports. 9: 1-13. 

Paszek MJ, Dufort C, Rubashkin MG, Davidson MW, Thorn KS, Liphardt JT, and Weaver VM.  2012.  Scanning angle interference microscopy reveals cell dynamics at the nano-scale.  Nature Methods 9: 825-7

Basic Usage
-----------

The method `fit_SAIM()` is the main data fitting routine. The current implementation supports analysis of SAIM image sequences acquired with one, two, or three (in development) laser excitation wavelengths. Currently only analysis of tif image stacks are supported.  Basic usage for analysis of a typical experiment with a single excitation wavelength is provided below:

```julia
using SAIMFitKit

# Optical parameters are stored in a composite structure
opt = SAIMOptics()
opt.nB = 1.33;	                #Refractive index of the ambient buffer / cytoplasm
opt.dOx = 1924.9;	              #Thickness of the SiO2 layer in units of nm
opt.λ_Ex_1 = 642.0;		        	#Excitation laser wavelength #1 in units of nm
opt.nOx_1 = 1.4719;		          #Refractive index of SiO2 at the excitation wavelength
opt.nSi_1 = 3.8660 + 0.017933im;    #The complex refractive index of Si at the excitation wavelength
```
The function fit_SAIM applies the per observation function `A*Pex(H, angle[i])+B` to the full dataset given by the various channels at each pixel in a SAIM image stack.  Each channel in the SAIM image stack corresponds to a microscope acquisition at an excitation laser incidence angle[i].  A, B, and H are the three fit parameters for single wavelength experiment.

```julia
# Fit parameters
p0 = [0.5, 1.0, 75.0]		      #Initial guesses for parameters [A, B, H]
lb = [0.0, 0.0, 0.0]			    #Lower bounds for [A, B, H]
ub = [10000., 10000., 200.]		#Upper bounds for [A, B, H]

# Image information
angles = range(5.0, length=32, stop=43.75)	    #Incidence angle in degrees for each image frame	
f = "21TR_5c"                                   #Name of image stack; .tif extension should not be included in the image name
path = "[user home folder]\\.julia\\dev\\SAIMFitKit\\example\\TestImages" #Path to image file location

# Conduct fit with global grid search (glb=true) with step size of 40 nm (step=40.0) and displaying the plot of the fit heights (show=true)
fit = fit_SAIM(path, f, opt, angles, p0, lb, ub, glb=true, step=40.0, show=true)

Existing Functionality
----------------------

`fit = fit_SAIM(file_path, file_name, optic, angles, init_params, lower_bounds, upper_bounds; kwargs...)`:

* `file_path`: string that provides the path to the image file location
* `file_name`: string that provides the name of the image file (.tif extension should not be included in name)`
* `optic`: composite type of imaging parameters ('SAIMOptic')
* `angles`: the incidence angles in degrees for the acquisition sequence 
* `lower_bounds`: lower bounds for the fit parameters (e.g. A, B, H)
* `upper_bounds`: lower bounds for the fit parameters (e.g. A, B, H)
* `kwargs`: addition parameters for fitting, such as `glb` (run global grid search when true), `step` (step size for grid search), `show` (display plot of fit heights when true), `color` (number of excitation wavelengths; typically one color)
* `fit`: composite type of results (`SAIMFitResult`)

Fit SAIM data using non-linear least squares (see LsqFit.jl) with either a local or global grid search. Results are automatically saved to an output folder. 

----

`plot_fit(image_file, JLD_file, row, col)`:

* `image_file`: string that provides the full path to the image file
* `JLD_name`: Output JLD file generated by `fit_SAIM`
* `row`: row number of image pixel to view fit result 
* `col`: column number of image pixel to view fit result 

Plot raw image intensity data and a previously generate fit result at an image pixel given by (row,col) 
