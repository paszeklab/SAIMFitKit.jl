# ---------------------------- fit_2c_local ---------------------------------------
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
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
# init_params: Intial parameters, [A1ₒ, B1ₒ, A2ₒ, B2ₒ, Hₒ]
# lower_bounds: Lower bounds for [A1, B1, A2, B2, H]
# upper_bounds: Lower bounds for [A1, B1, A2, B2, H]
# disp: When true, a heatmap for the fit heights will be displayed
# OUTPUTS: 1) Heatmap with fitted heights; 2) 16-bit image image with intensities
# proportional to best-fit pixel height -  Pixel intensity = 100*height in nm;
# 3) JLD file with results for each pixel - Fields are "A1," "B1," "A2," "B2,"
# "H," and "errors"; errors are the standard errors for fit [A1, B1, A2, B2, H]
function fit_2c_local(file_path::String, file_name::String, optic1::SAIMOptics, optic2::SAIMOptics, 
	angles::AbstractArray, init_params::AbstractArray, lower_bounds::AbstractArray,
	upper_bounds::AbstractArray, disp::Bool=false)

	#Calculate the optical model constants for each image frame angle
	constants = calculate_constants_2c(optic1, optic2, angles)

	#------------- OPEN IMAGE STACK AND INITIALIZE CONTAINERS -------------
	file = file_path*"\\"*file_name*".tif"
	#Open an image
	img1_HWC = load(file);								#image stored as height X width X channels
	img_CHW = permutedims(img1_HWC, (3, 1, 2))  		#image stored as channels X height X width; CHW indexing more memory friendly
	channels = channelview(float64.(img_CHW))*65535; 	#Each channel corresponds to a specific excitation angle
	rows = size(channels,2)  							#Image size is (number channels X rows X cols)
	cols = size(channels,3)								#Image size is (number channels X rows X cols)

	#Initialize container to store the fitted heights and fitting results
	fit_params = Array{Float64}(undef, rows, cols, 5)
	fit_errors = fill(-1.0, rows, cols, 5)

	#-------------- CONDUCT THE FITTING ----------------------
	println("Conducting fit for ", file_name)  #Threads.@threads
	counter = 0
	@time Threads.@threads for i = 1:rows
		for j = 1:cols
			ydata = channels[:,i,j]

			result = curve_fit(model_2c, jacobian_model_2c, constants, ydata,
				init_params, lower=lower_bounds, upper=upper_bounds)

			fit_params[i,j,:] = coef(result)

			#Calculate the standard errors for the fitted parameters
			try
           		fit_errors[i,j,:] = stderror(result)
       		catch e
           		counter = counter + 1
       		end
		end
	end
	println("  ", (rows*cols - counter), " out of ", (rows*cols), " pixels fit within bounds")

	#-------------- GENERATE THE OUTPUTS ----------------------
	println("Generating outputs for ", file_name)

	#Make an output directory if needed
	dir = file_path*"\\output\\"*file_name

	if !isdir(file_path*"\\output")
		mkdir(file_path*"\\output")
	end
	
	if !isdir(dir)
		mkdir(dir)
	end

	#Save height information in a 16-bit image; Pixel intensity = 100*height in nm
	save(File(format"PNG", dir*"\\"*file_name*"_H.png"), n0f16.((fit_params[:,:,5]/1000.)*(100000. / 65535.)))

	#Save the fit information as a JLD file (includes fit parameters, errors, and metadata)
	save(File(format"JLD", dir*"\\"*file_name*"_results.jld"),"A1", fit_params[:,:,1], "B1",
		fit_params[:,:,2], "A2", fit_params[:,:,3], "B2", fit_params[:,:,4],
		"H",fit_params[:,:,5], "errors", fit_errors, "constants", constants,
		"angles", angles, "ip", init_params, "lb", lower_bounds, "ub",
		upper_bounds, "optic1", optic1, "optic2", optic2,
		"type", "two-color global")

	#Save the fit parameters and height standard err as CSV files
	if !isdir(dir*"\\CSV")
		mkdir(dir*"\\CSV")
	end
	CSV.write(dir*"\\CSV\\"*file_name*"_A1.csv",  DataFrame(fit_params[:,:,1]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_B1.csv",  DataFrame(fit_params[:,:,2]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_A2.csv",  DataFrame(fit_params[:,:,3]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_B2.csv",  DataFrame(fit_params[:,:,4]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H.csv",  DataFrame(fit_params[:,:,5]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H_sterr.csv",  DataFrame(fit_errors[:,:,5]), writeheader=false)

	#Construct a height map and save
	hm = heatmap(fit_params[:,:,5], yflip=true, colorbar_title="height (nm)",
		color=:RdYlBu_10, aspect_ratio=:equal, xaxis = (nothing, false),
		 yaxis = (nothing, false))
	ylims!((1,rows))
	xlims!((1,cols))
	display(hm)
	savefig(hm, dir*"\\"*file_name*"_Hmap.pdf")
	println("Completed ", file_name)
	println()
	if disp
		display(hm)
	end
end

# ---------------------------- fit_2c_global ----------------------------------
# Function to conduct global non-linear least squares curve fitting for each
# pixel in a two-color SAIM image stack
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
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
# init_params: Intial parameters, [A1ₒ, B1ₒ, A2ₒ, B2ₒ]
# lower_bounds: Lower bounds for [A1, B1, A2, B2, H]
# upper_bounds: Lower bounds for [A1, B1, A2, B2, H]
# step: Step size for initial guesses for H in grid search
# disp: When true, a heatmap for the fit heights will be displayed
# OUTPUTS: 1) Heatmap with fitted heights; 2) 16-bit image image with intensities
# proportional to best-fit pixel height -  Pixel intensity = 100*height in nm;
# 3) JLD file with results for each pixel - Fields are "A1," "B1," "A2," "B2,"
# "H," and "errors"; errors are the standard errors for fit [A1, B1, A2, B2, H]
function fit_2c_global(file_path::String, file_name::String, optic1::SAIMOptics, optic2::SAIMOptics, 
	angles::AbstractArray, init_params::AbstractArray, lower_bounds::AbstractArray,
	upper_bounds::AbstractArray, step::Float64=40., disp::Bool=false)

	#Calculate the optical model constants for each image frame angle
	constants = calculate_constants_2c(optic1, optic2, angles)

	#------------- OPEN IMAGE STACK AND INITIALIZE CONTAINERS -------------
	file = file_path*"\\"*file_name*".tif"
	#Open an image
	img1_HWC = load(file);								#image stored as height X width X channels
	img_CHW = permutedims(img1_HWC, (3, 1, 2))  		#image stored as channels X height X width; CHW indexing more memory friendly
	channels = channelview(float64.(img_CHW))*65535; 	#Each channel corresponds to a specific excitation angle
	rows = size(channels,2)  							#Image size is (number channels X rows X cols)
	cols = size(channels,3)								#Image size is (number channels X rows X cols)

	#Initialize container to store the fitted heights and fitting results
	fit_params = fill(0.0, rows, cols, 5)
	fit_errors = fill(-1.0, rows, cols, 5)

	#-------------- CONDUCT THE FITTING ----------------------
	println("Conducting fit for ", file_name)  #Threads.@threads
	guesses = collect( lower_bounds[5]:step:upper_bounds[5] )		#initial guesses for heights
	params = Array{Float64}(undef, 5)
	params[1] = init_params[1]
	params[2] = init_params[2]
	params[3] = init_params[3]
	params[4] = init_params[4]

	@time Threads.@threads for i = 1:rows
		for j = 1:cols
			ydata = channels[:,i,j]
			err_lo = 1.0e10				#Set arbitrarily high
			#Global fitting routine - grid search
			for Hₒ in guesses
				params[5] = Hₒ
				result = curve_fit(model_2c, jacobian_model_2c, constants, ydata,
					params, lower=lower_bounds, upper=upper_bounds)
				#Calculate the standard errors for the fitted parameters
				if result.converged
					try
						err = stderror(result)
						if err[5] < err_lo
							fit_params[i,j,:] = coef(result)
							fit_errors[i,j,:] = err
							err_lo = err[5]
						end
					catch e
					end
				end
			end
		end
	end

	#-------------- GENERATE THE OUTPUTS ----------------------
	println("Generating outputs for ", file_name)

	#Make an output directory if needed
	dir = file_path*"\\output\\"*file_name

	if !isdir(file_path*"\\output")
		mkdir(file_path*"\\output")
	end
	
	if !isdir(dir)
		mkdir(dir)
	end

	#Save height information in a 16-bit image; Pixel intensity = 100*height in nm
	save(File(format"PNG", dir*"\\"*file_name*"_H.png"), n0f16.((fit_params[:,:,5]/1000.)*(100000. / 65535.)))

	#Save the fit information as a JLD file (includes fit parameters, errors, and metadata)
	save(File(format"JLD", dir*"\\"*file_name*"_results.jld"),"A1", fit_params[:,:,1], "B1",
		fit_params[:,:,2], "A2", fit_params[:,:,3], "B2", fit_params[:,:,4],
		"H",fit_params[:,:,5], "errors", fit_errors, "constants", constants,
		"angles", angles, "ip", init_params, "lb", lower_bounds, "ub",
		upper_bounds, "step", step, "optic1", optic1, "optic2", optic2,
		"type", "two-color global")

	#Save the fit parameters and height standard err as CSV files
	if !isdir(dir*"\\CSV")
		mkdir(dir*"\\CSV")
	end
	CSV.write(dir*"\\CSV\\"*file_name*"_A1.csv",  DataFrame(fit_params[:,:,1]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_B1.csv",  DataFrame(fit_params[:,:,2]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_A2.csv",  DataFrame(fit_params[:,:,3]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_B2.csv",  DataFrame(fit_params[:,:,4]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H.csv",  DataFrame(fit_params[:,:,5]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H_sterr.csv",  DataFrame(fit_errors[:,:,5]), writeheader=false)

	#Construct a height map and save
	hm = heatmap(fit_params[:,:,5], yflip=true, colorbar_title="height (nm)",
		color=:RdYlBu_10, aspect_ratio=:equal, xaxis = (nothing, false),
		 yaxis = (nothing, false))
	ylims!((1,rows))
	xlims!((1,cols))
	savefig(hm, dir*"\\"*file_name*"_Hmap.png")
	savefig(hm, dir*"\\"*file_name*"_Hmap.pdf")
	println("  Completed ", file_name)
	println()
	if disp
		display(hm)
	end
end

# ------------------------ calculate_constants_2c -----------------------------
# Function that returns constants for calculating Fresnel coefficients (real
# and imag components) and phase constants in two-color SAIM fitting routines
# INPUTS:
# optic1: Optical struct containing parameters for excitation wavelength #1
# optic2: Optical struct containing parameters for excitation wavelength #2
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
function calculate_constants_2c(optic1, optic2, angles)

	#Initialize containers for holding the constants
	num_angles = length(angles)
	constants = Array{Float64}(undef, 2*num_angles, 3)

	#Calculate the constants for the first laser wavelength
	wavelength = optic1.lambda_Ex	#The Excitation wavelength
	dOx = optic1.d_Oxide		#The thickness of the SiO2 layer in units of nm
	nB = optic1.n_Buffer		#The refractive index of the ambient media / cytoplasm
	nOx = optic1.n_Oxide		#The refractive index of SiO2
	nSi = optic1.n_Silicon		#The efractive index of Si

	inds = axes(angles, 1)
	for i = inds
		#Calculate the incidence angle of excitation light in each material
		θB = angles[i]						#Incidence angle in buffer
		θOx = asin(nB*sin(θB)/nOx) 			#Incidence angle in silicon oxide
		θSi = asin(nOx*sin(θOx)/nSi)		#Incidence angle in silicon

		#Calculate the Fresnesl coefficient (perpendicular components) for the layer system (Si02, Si)
		kOx = 2*pi*nOx/wavelength
		lOx = kOx*dOx*cos(θOx)
		pSi = nSi*cos(θSi)
		pOx = nOx*cos(θOx)
		pB = nB*cos(θB)

		m11TE = cos(lOx)
		m12TE = -(im/pOx)*sin(lOx)
		m21TE = -im*pOx*sin(lOx)
		m22TE = m11TE

		rTE = ((m11TE + m12TE*pSi)*pB - (m21TE + m22TE*pSi)) / ((m11TE + m12TE*pSi)*pB + (m21TE + m22TE*pSi));

		constants[i,1] = real(rTE)		#Real component of rTE
		constants[i,2] = imag(rTE)		#Imaginary compnent of rTE
		constants[i,3] = (4*pi*nB/wavelength)*cos(θB) 	#phase shift/height
	end

	#Calculate the constants for the second laser wavelength
	wavelength = optic2.lambda_Ex	#The Excitation wavelength
	dOx = optic2.d_Oxide		#The thickness of the SiO2 layer in units of nm
	nB = optic2.n_Buffer		#The refractive index of the ambient media / cytoplasm
	nOx = optic2.n_Oxide		#The refractive index of SiO2
	nSi = optic2.n_Silicon		#The efractive index of Si

	#inds = axes(angles, 1)
	for i = inds
		#Calculate the incidence angle of excitation light in each material
		θB = angles[i]						#Incidence angle in buffer
		θOx = asin(nB*sin(θB)/nOx) 			#Incidence angle in silicon oxide
		θSi = asin(nOx*sin(θOx)/nSi)		#Incidence angle in silicon

		#Calculate the Fresnesl coefficient (perpendicular components) for the layer system (Si02, Si)
		kOx = 2*pi*nOx/wavelength
		lOx = kOx*dOx*cos(θOx)
		pSi = nSi*cos(θSi)
		pOx = nOx*cos(θOx)
		pB = nB*cos(θB)

		m11TE = cos(lOx)
		m12TE = -(im/pOx)*sin(lOx)
		m21TE = -im*pOx*sin(lOx)
		m22TE = m11TE

		rTE = ((m11TE + m12TE*pSi)*pB - (m21TE + m22TE*pSi)) / ((m11TE + m12TE*pSi)*pB + (m21TE + m22TE*pSi));

		constants[i+num_angles,1] = real(rTE)		#Real component of rTE
		constants[i+num_angles,2] = imag(rTE)		#Imaginary compnent of rTE
		constants[i+num_angles,3] = (4*pi*nB/wavelength)*cos(θB) 	#phase shift/height
	end

	return constants
end

#Optical model for two-color experiments: y = A*Pex(H,theta)+B
#p = [A1, B1, A2, B2, H]
#x[i,:] holds the optical model constants for each excitation incidence angle / wavelength combination
function model_2c(x, p)

	y= Vector{Float64}(undef, length(x[:,1]))
	num_angles::Int16 = length(x[:,1])/2

	#First excitation wavelength
	for i = 1:num_angles
		c = x[i,1]			#Real part of rTE
		d = x[i,2]			#Imaginary part of rTE
		ϕₕ = x[i,3]			#ϕₕ*H is the phase shift

		Pex = 1 + 2*c*cos(ϕₕ*p[5]) - 2*d*sin(ϕₕ*p[5]) + c*c + d*d
        y[i] = p[1]*Pex + p[2]
	end

	#Second excitation wavelength
	for i = (num_angles+1):(2*num_angles)
		c = x[i,1]			#Real part of rTE
		d = x[i,2]			#Imaginary part of rTE
		ϕₕ = x[i,3]		#ϕₕ*H is the phase shift

		Pex = 1 + 2*c*cos(ϕₕ*p[5]) - 2*d*sin(ϕₕ*p[5]) + c*c + d*d
	    y[i] = p[3]*Pex + p[4]
    end
    y
end

#Generate the analytical Jacobian for model2c
#p = [A1, B1, A2, B2, H]
#x[i,:] holds the optical model constants for each excitation incidence angle / wavelength combination
function jacobian_model_2c(x,p)

	J = Array{Float64}(undef, length(x[:,1]), length(p))
	num_angles::Int16 = length(x[:,1])/2

	#First excitation wavelength
	for i = 1:num_angles
		c = x[i,1]			#Real part of rTE
		d = x[i,2]			#Imaginary part of rTE
		ϕₕ = x[i,3]			#ϕₕ*H is the phase shift

        J[i,1] = 1 + 2*c*cos(ϕₕ*p[5]) - 2*d*sin(ϕₕ*p[5]) + c*c + d*d 	  	 #dmodel/dA1
        J[i,2] = 1.0										 				#dmodel/dB1
		J[i,3] = 0.0														#dmodel/dA2
		J[i,4] = 0.0														#dmodel/dB2
		J[i,5] = -2*p[1]*ϕₕ*(c*sin(ϕₕ*p[5]) + d*cos(ϕₕ*p[5])) 				 #dmodel/dH
    end

	#Second excitation wavelength
	for i = (num_angles+1):(2*num_angles)
		c = x[i,1]			#Real part of rTE
		d = x[i,2]			#Imaginary part of rTE
		ϕₕ = x[i,3]			#ϕₕ*H is the phase shift

		J[i,1] = 0.0 	  	 												#dmodel/dA1
		J[i,2] = 0.0										 				#dmodel/dB1
		J[i,3] = 1 + 2*c*cos(ϕₕ*p[5]) - 2*d*sin(ϕₕ*p[5]) + c*c + d*d		 #dmodel/dA2
		J[i,4] = 1.0														#dmodel/dB2
		J[i,5] = -2*p[3]*ϕₕ*(c*sin(ϕₕ*p[5]) + d*cos(ϕₕ*p[5])) 				 #dmodel/dH
	end
    J
end