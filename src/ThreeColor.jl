# ---------------------------- fit_3c_local ---------------------------------------
# Function to conduct local non-linear least squares curve fitting for each
# pixel in a two-color SAIM image stack
# Images sequences should be saved as a tif stack
# Sequences should be organized as (λ1₁...λ1ₙ λ2₁...λ2ₙ) where n indicates the
# nth and last angle; λ1ᵢ and λ2ᵢ are the ith images acquired with the first and
# and second excitation wavelengths, respectively
# INPUTS:
# name: A string that stores base file name of tif stack
# file_ext: file extension for iumage stack; Only ".tif" has been tested so far
# optic1: Optical struct containing parameters for excitation wavelength #1
# optic2: Optical struct containing parameters for excitation wavelength #2
# optic3: Optical struct containing parameters for excitation wavelength #3
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
# init_params: Intial parameters, [A1ₒ, B1ₒ, A2ₒ, B2ₒ, A3ₒ, B3ₒ]
# lower_bounds: Lower bounds for [A1, B1, A2, B2, A3, B3, H]
# upper_bounds: Lower bounds for [A1, B1, A2, B2, A3, B3, H]
# disp: When true, a heatmap for the fit heights will be displayed
# OUTPUTS: 1) Heatmap with fitted heights; 2) 16-bit image image with intensities
# proportional to best-fit pixel height -  Pixel intensity = 100*height in nm;
# 3) JLD file with results for each pixel - Fields are "A1," "B1," "A2," "B2," "A3", "B3",
# "H," and "errors"; errors are the standard errors for fit [A1, B1, A2, B2, H]
function fit_3c_local(file_path::String, file_name::String, optic1::SAIMOptics, optic2::SAIMOptics, 
	optic3::SAIMOptics, angles::AbstractArray, init_params::AbstractArray, lower_bounds::AbstractArray,
	upper_bounds::AbstractArray, disp::Bool=false)

    #Add code
end

# ---------------------------- fit_3c_global ----------------------------------
# Function to conduct global non-linear least squares curve fitting for each
# pixel in a three-color SAIM image stack
# Executes a simple grid search, stepping through initial guesses for H; search
# is conducted along H from lower_bounds to upper_bounds
# Images sequences should be saved as a tif stack
# Sequences should be organized as (λ1₁...λ1ₙ λ2₁...λ2ₙ) where n indicates the
# nth and last angle; λ1ᵢ and λ2ᵢ are the ith images acquired with the first and
# and second excitation wavelengths, respectively
# INPUTS:
# name: A string that stores base file name of tif stack
# file_ext: file extension for iumage stack; Only ".tif" has been tested so far
# optic1: Optical struct containing parameters for excitation wavelength #1
# optic2: Optical struct containing parameters for excitation wavelength #2
# optic3: Optical struct containing parameters for excitation wavelength #3
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
# init_params: Intial parameters, [A1ₒ, B1ₒ, A2ₒ, B2ₒ, A3ₒ, B3ₒ]
# lower_bounds: Lower bounds for [A1, B1, A2, B2, A3, B3, H]
# upper_bounds: Lower bounds for [A1, B1, A2, B2, A3, B3, H]
# step: Step size for initial guesses for H in grid search
# disp: When true, a heatmap for the fit heights will be displayed
# OUTPUTS: 1) Heatmap with fitted heights; 2) 16-bit image image with intensities
# proportional to best-fit pixel height -  Pixel intensity = 100*height in nm;
# 3) JLD file with results for each pixel - Fields are "A1," "B1," "A2," "B2," "A3", "B3",
# "H," and "errors"; errors are the standard errors for fit [A1, B1, A2, B2, H]
function fit_3c_global(file_path::String, file_name::String, optic1::SAIMOptics, optic2::SAIMOptics, 
	optic3::SAIMOptics, angles::AbstractArray, init_params::AbstractArray, lower_bounds::AbstractArray,
	upper_bounds::AbstractArray, step::Float64=40., disp::Bool=false)

    #Add code
end

# ------------------------ Calculate_Constants_3c -----------------------------
# Function that returns constants for calculating Fresnel coefficients (real
# and imag components) and phase constants in two-color SAIM fitting routines
# INPUTS:
# optic1: Optical struct containing parameters for excitation wavelength #1
# optic2: Optical struct containing parameters for excitation wavelength #2
# optic3: Optical struct containing parameters for excitation wavelength #3
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
function calculate_constants_3c(optic1, optic2, optic3, angles)

    #Add code
end




#Optical model for three-color experiments: y = A*Pex(H,theta)+B
#p = [A1, B1, A2, B2, A3, B3, H]
#x[i,:] holds the optical model constants for each excitation incidence angle / wavelength combination
function model_3c(x, p)
    #Add code
end


#Generate the analytical Jacobian for model_3c
#p = [A1, B1, A2, B2, A3, B3, H]
#x[i,:] holds the optical model constants for each excitation incidence angle / wavelength combination
function jacobian_model_3c(x,p)
    #Add code
end