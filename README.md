# SAIMFitKit.jl
Julia implementation of the Scanning Angle Interference Microscopy (SAIM) fitting algorithms

## Description
Fitting algorithm reconstructs sample height from a stack of scanning angle interference images. More details on the fiting algorithms can be found in the following papers:

Colville MJ, Park S, Zipfel WR, and Paszek MJ. High-speed device synchronization in optical microscopy with an open-source hardware control platform. Scientific Reports. 9: 1-13. 

Paszek MJ, Dufort C, Rubashkin MG, Davidson MW, Thorn KS, Liphardt JT, and Weaver VM.  2012.  Scanning angle interference microscopy reveals cell dynamics at the nano-scale.  Nature Methods 9: 825-7

Existing Functionality
----------------------

`fit = curve_fit(model, [jacobian], x, y, [w,] p0; kwargs...)`:

* `model`: function that takes two arguments (x, params)
* `jacobian`: (optional) function that returns the Jacobian matrix of `model`
* `x`: the independent variable
* `y`: the dependent variable that constrains `model`
* `w`: (optional) weight applied to the residual; can be a vector (of `length(x)` size or empty) or matrix (inverse covariance matrix)
* `p0`: initial guess of the model parameters
* `kwargs`: tuning parameters for fitting, passed to `levenberg_marquardt`, such as `maxIter`, `show_trace` or `lower` and `upper` bounds
* `fit`: composite type of results (`LsqFitResult`)


This performs a fit using a non-linear iteration to minimize the (weighted) residual between the model and the dependent variable data (`y`). The weight (`w`) can be neglected (as per the example) to perform an unweighted fit. An unweighted fit is the numerical equivalent of `w=1` for each point  (although unweighted error estimates are handled differently from weighted error estimates even when the weights are uniform).

----

`sigma = stderror(fit; atol, rtol)`:

* `fit`: result of curve_fit (a `LsqFitResult` type)
* `atol`: absolute tolerance for negativity check
* `rtol`: relative tolerance for negativity check

This returns the error or uncertainty of each parameter fit to the model and already scaled by the associated degrees of freedom.  Please note, this is a LOCAL quantity calculated from the jacobian of the model evaluated at the best fit point and NOT the result of a parameter exploration.

If no weights are provided for the fits, the variance is estimated from the mean squared error of the fits. If weights are provided, the weights are assumed to be the inverse of the variances or of the covariance matrix, and errors are estimated based on these and the jacobian, assuming a linearization of the model around the minimum squared error point.
