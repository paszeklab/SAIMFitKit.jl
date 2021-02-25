function test_working()
	println("yes")
end

# ---------------------------- SAIM_fit ---------------------------------------
# Function to conduct non-linear least squares curve fitting for each pixel
# in a one-color SAIM image stack
# Images sequences should be saved as a tif stack
# Sequences should be organized as (λ1₁...λ1ₙ) where n indicates the
# nth and last angle; λ1ᵢ is the ith image acquired with the excitation
# wavelength
# INPUTS:
# name: A string that stores base file name of tif stack
# file_ext: file extension for iumage stack; Only ".tif" has been tested so far
# optic: Optical struct containing parameters for the excitation wavelength
# angles: 1D array containing incidence angles in radians for the acquisition
# 	sequence
# init_params: Intial parameters, [Aₒ, Bₒ, Hₒ]
# lower_bounds: Lower bounds for [A, B, H]
# upper_bounds: Lower bounds for [A, B, H]
# disp: When true, a heatmap for the fit heights will be displayed
# OUTPUTS: 1) Heatmap with fitted heights; 2) 16-bit image image with intensities
# proportional to best-fit pixel height -  Pixel intensity = 100*height in nm;
# 3) JLD file with results for each pixel - Fields are "A," "B,"  "H,"
# and "errors"; errors are the standard errors for fit [A, B, H]
function fit_SAIM(file_path::String, file_name::String, optic::SAIMOptics, angles::AbstractArray, 
	init_params::AbstractArray, lower_bounds::AbstractArray, upper_bounds::AbstractArray, disp::Bool=false)

	#Calculate the optical model constants for each image frame angle
	constants = calculate_constants(optic, angles)

	#------------- OPEN IMAGE STACK AND INITIALIZE CONTAINERS -------------
	file = file_path*"\\"*file_name*".tif"
	#Open an image
	img1_HWC = load(file);						#image stored as height X width X channels
	img_CHW = permutedims(img1_HWC, (3, 1, 2))  		#image stored as channels X height X width; CHW indexing more memory friendly
	channels = channelview(float64.(img_CHW))*65535; 	#Each channel corresponds to a specific excitation angle
	rows = size(channels,2)  							#Image size is (number channels X rows X cols)
	cols = size(channels,3)								#Image size is (number channels X rows X cols)

	#Initialize container to store the fitted heights and fitting results
	fit_params = Array{Float64}(undef, rows, cols, 3)
	fit_errors = fill(-1.0, rows, cols, 3)

	#-------------- CONDUCT THE FITTING ----------------------
	println("Conducting fit for ", file_name)  #Threads.@threads
	counter = 0
	@time Threads.@threads for i = 1:rows
		for j = 1:cols
			ydata = channels[:,i,j]

			result = curve_fit(model_1c, jacobian_model_1c, constants, ydata,
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
	save(File(format"PNG", dir*"\\"*file_name*"_H.png"), n0f16.((fit_params[:,:,3]/1000.)*(100000. / 65535.)))

	#Save the fit information as a JLD file (includes fit parameters, errors, and metadata)
	save(File(format"JLD", dir*"\\"*file_name*"_results.jld"),"A", fit_params[:,:,1], "B",
		fit_params[:,:,2], "H",fit_params[:,:,3], "errors", fit_errors,
		"constants", constants, "angles", angles, "ip", init_params,
		"lb", lower_bounds, "ub", upper_bounds, "optic",
		optic, "type", "one-color local")

	#Save the fit parameters and height standard err as CSV files
	if !isdir(dir*"\\CSV")
		mkdir(dir*"\\CSV")
	end
	CSV.write(dir*"\\CSV\\"*file_name*"_A.csv",  DataFrame(fit_params[:,:,1]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_B.csv",  DataFrame(fit_params[:,:,2]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H.csv",  DataFrame(fit_params[:,:,3]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H_sterr.csv",  DataFrame(fit_errors[:,:,3]), writeheader=false)

	#Construct a height map and save
	hm = heatmap(fit_params[:,:,3], yflip=true, colorbar_title="height (nm)",
		color=:RdYlBu_10, aspect_ratio=:equal, xaxis = (nothing, false),
		 yaxis = (nothing, false))
	ylims!((1,rows))
	xlims!((1,cols))

	savefig(hm, dir*"\\"*file_name*"_Hmap.pdf")
	println("  Completed ", file_name)
	println()
	if disp
		display(hm)
	end
end

# ---------------------------- fit_SAIM_g ----------------------------------
# Function to conduct global non-linear least squares curve fitting for each
# pixel in a two-color SAIM image stack
# Executes a simple grid search, stepping through initial guesses for H; search
# is conducted along H from lower_bounds to upper_bounds
# Images sequences should be saved as a tif stack
# Sequences should be organized as (λ1₁...λ1ₙ) where n indicates the
# nth and last angle; λ1ᵢ is the ith image acquired with the excitation
# wavelength
# INPUTS:
# name: A string that stores base file name of tif stack
# file_ext: file extension for iumage stack; Only ".tif" has been tested so far
# optic: Optical struct containing parameters for the excitation wavelength
# angles: 1D array containing incidence angles in radians for the acquisition
# sequence
# init_params: Intial parameters, [Aₒ, Bₒ]
# lower_bounds: Lower bounds for [A, B, H]
# upper_bounds: Lower bounds for [A, B, H]
# step: Step size for initial guesses for H in grid search
# disp: When true, a heatmap for the fit heights will be displayed
# OUTPUTS: 1) Heatmap with fitted heights; 2) 16-bit image image with intensities
# proportional to best-fit pixel height -  Pixel intensity = 100*height in nm;
# 3) JLD file with results for each pixel - Fields are "A," "B,"  "H,"
# and "errors"; errors are the standard errors for fit [A, B, H]
function fit_SAIM_g(file_path::String, file_name::String, optic::SAIMOptics, angles::AbstractArray, 
	init_params::AbstractArray, lower_bounds::AbstractArray, upper_bounds::AbstractArray, 
	step::Float64=40., disp::Bool=false)


	#Calculate the optical model constants for each image frame angle
	constants = calculate_constants(optic, angles)

	#------------- OPEN IMAGE STACK AND INITIALIZE CONTAINERS -------------
	file = file_path*"\\"*file_name*".tif"
	#Open an image
	img1_HWC = load(file);						#image stored as height X width X channels
	img_CHW = permutedims(img1_HWC, (3, 1, 2))  		#image stored as channels X height X width; CHW indexing more memory friendly
	channels = channelview(float64.(img_CHW))*65535; 	#Each channel corresponds to a specific excitation angle
	rows = size(channels,2)  							#Image size is (number channels X rows X cols)
	cols = size(channels,3)								#Image size is (number channels X rows X cols)

	#Initialize container to store the fitted heights and fitting results
	fit_params = fill(0.0, rows, cols, 3)
	fit_errors = fill(-1.0, rows, cols, 3)

	#-------------- CONDUCT THE FITTING ----------------------
	println("Conducting fit for ", file_name)  #Threads.@threads
	guesses = collect( lower_bounds[3]:step:upper_bounds[3] )		#initial guesses for heights
	params = Array{Float64}(undef, 3)
	params[1] = init_params[1]
	params[2] = init_params[2]

	@time Threads.@threads for i = 1:rows
		for j = 1:cols
			ydata = channels[:,i,j]
			err_lo = 1.0e10				#Set arbitrarily high
			#Global fitting routine - grid search
			for Hₒ in guesses
				params[3] = Hₒ
				result = curve_fit(model_1c, jacobian_model_1c, constants, ydata,
					params, lower=lower_bounds, upper=upper_bounds)
				#Calculate the standard errors for the fitted parameters
				if result.converged
					try
						err = stderror(result)
						if err[3] < err_lo
							fit_params[i,j,:] = coef(result)
							fit_errors[i,j,:] = err
							err_lo = err[3]
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
	save(File(format"PNG", dir*"\\"*file_name*"_H.png"), n0f16.((fit_params[:,:,3]/1000.)*(100000. / 65535.)))

	#Save the fit information as a JLD file (includes fit parameters, errors, and metadata)
	save(File(format"JLD", dir*"\\"*file_name*"_results.jld"),"A", fit_params[:,:,1], "B",
		fit_params[:,:,2], "H",fit_params[:,:,3], "errors", fit_errors,
		"constants", constants, "angles", angles, "ip", init_params,
		"lb", lower_bounds, "ub", upper_bounds, "step", step, "optic",
		optic, "type", "one-color global")

	#Save the fit parameters and height standard err as CSV files
	if !isdir(dir*"\\CSV")
		mkdir(dir*"\\CSV")
	end
	CSV.write(dir*"\\CSV\\"*file_name*"_A.csv",  DataFrame(fit_params[:,:,1]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_B.csv",  DataFrame(fit_params[:,:,2]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H.csv",  DataFrame(fit_params[:,:,3]), writeheader=false)
	CSV.write(dir*"\\CSV\\"*file_name*"_H_sterr.csv",  DataFrame(fit_errors[:,:,3]), writeheader=false)

	#Construct a height map and save
	hm = heatmap(fit_params[:,:,3], yflip=true, colorbar_title="height (nm)",
		color=:RdYlBu_10, aspect_ratio=:equal, xaxis = (nothing, false),
		 yaxis = (nothing, false))
	ylims!((1,rows))
	xlims!((1,cols))
	savefig(hm, dir*"\\"*file_name*"_Hmap.pdf")
	println("  Completed ", file_name)
	println()
	if disp
		display(hm)
	end

end

# ------------------------ Calculate_Constant -----------------------------
# Function that returns constants for calculating Fresnel coefficients (real
# and imag components) and phase constants in one-color SAIM fitting routines
# INPUTS:
# optic: Optical struct containing parameters for the excitation wavelength
# angles: 1D array containing incidence angles in radians for the acquisition
# sequence
function calculate_constants(optic, angles)

	#Initialize containers for holding the constants
	constants = Array{Float64}(undef, length(angles), 3)

	wavelength = optic.lambda_Ex	#The Excitation wavelength
	dOx = optic.d_Oxide		#The thickness of the SiO2 layer in units of nm
	nB = optic.n_Buffer		#The refractive index of the ambient media / cytoplasm
	nOx = optic.n_Oxide		#The refractive index of SiO2
	nSi = optic.n_Silicon		#The efractive index of Si

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

	return constants
end

#Optical model for one-color experiments: y = A*Pex(H,theta)+B
#p = [A, B, H]
#x[i,:] holds the optical model constants for each excitation incidence angle
function model_1c(x, p)
    y= Vector{Float64}(undef, length(x[:,1]))
    for i in eachindex(y)
		c = x[i,1]			#Real part of rTE
		d = x[i,2]			#Imaginary part of rTE
		ϕₕ = x[i,3]		#ϕₕ*H is the phase shift
		Pex = 1 + 2*c*cos(ϕₕ*p[3]) - 2*d*sin(ϕₕ*p[3]) + c*c + d*d
        y[i] = p[1]*Pex + p[2]
    end
    y
end

#Generate the analytical Jacobian for model
#p = [A, B, H]
#x[i,:] holds the optical model constants for each excitation incidence angle
function jacobian_model_1c(x,p)
    J = Array{Float64}(undef, length(x[:,1]), length(p))
    for i = 1:length(x[:,1])
		c = x[i,1]			#Real part of rTE
		d = x[i,2]			#Imaginary part of rTE
		ϕₕ = x[i,3]			#ϕₕ*H is the phase shift

        J[i,1] = 1 + 2*c*cos(ϕₕ*p[3]) - 2*d*sin(ϕₕ*p[3]) + c*c + d*d 	  	 #dmodel/dA
        J[i,2] = 1										 		#dmodel/dB
		J[i,3] = -2*p[1]*ϕₕ*(c*sin(ϕₕ*p[3]) + d*cos(ϕₕ*p[3])) 	#dmodel/dH
    end
    J
end