#Function to show example fits for a given pixel in a SAIM image sequence
function plot_fit(image_file, JLD_file, row, col)

	println(image_file)
	println(JLD_file)

	#Load the constants from the JLD_file
	constants = load(JLD_file, "constants")
	angles = load(JLD_file, "angles")
	angles_deg = angles*360/(2*pi)		#convert from radians to degrees

	#Load the image sequence and extract the channel intensities at pixel (row, col)
	img_HWC = load(image_file);
    channels = channelview(float64.(img_HWC))*65535; 	#Each channel corresponds to a specific excitation angle
    ydata = channels[row,col,:]

    #Load the fit results
    A = load(JLD_file, "A")
    B = load(JLD_file, "B")
	H = load(JLD_file, "H")
    constants = load(JLD_file, "constants")
    p = [A[row,col], B[row,col], H[row,col]]
    
	#Plot the example fit
	yfit = model_1c(constants, p) #model(constants, p)
    colors = palette(:tab10)
	plt = plot(angles_deg, ydata, label="data", seriestype = :scatter, markercolor = colors[1])
	plot!(angles_deg, yfit, label="fit", linecolor = colors[4], linewidth = 1.5)
    xlabel!("Incidence angle (degrees)")
	ylabel!("Intensity (A.U.)")
	
    savefig(plt, "example_fit.png")
	savefig(plt, "example_fit.pdf")
	savefig(plt, "example_fit.svg")
    display(plt)
end

#=calculate_constants_1c(optic, angles)

	#Initialize containers for holding the constants
	constants = Array{Float64}(undef, length(angles), 3)

	wavelength = optic.λ_Ex_1	#The Excitation wavelength
	dOx = optic.dOx 			#The thickness of the SiO2 layer in units of nm
	nB = optic.nB				#The refractive index of the ambient media / cytoplasm
	nOx = optic.nOx_1			#The refractive index of SiO2
	nSi = optic.nSi_1			#The efractive index of Si
=#

#Plot a SAIM curve
function plot_curve_1c(optic::SAIMOptics, heights::AbstractArray, angle_range::AbstractArray)

	#-------------- ANGLES  ---------------------
	angles_deg = range(angle_range[1], length=100, stop=angle_range[2])	#Incendece angle in degree for each image frame
	angles_rad = angles_deg*2*pi/360 					#Incendece angle in radians for each image frame
	num_angles = length(angles_rad)

	#-------------- SIMULATION----------------------
	colors = palette(:tab10)

	constants = calculate_constants_1c(optic, angles_rad) 
	println(heights)
	
	#Optical model for one-color experiments: y = A*Pex(H,theta)+B
	#p = [A, B, H]
	#x[i,:] holds the optical model constants for each excitation incidence angle
	#p = [A, B, H]
	plt = plot()
	for h in heights

		p = [1.0, 0.0, h]
		Pex = model_1c(constants, p)
		#plot!(angles_deg, Pex, linecolor = colors[1], linewidth = 2, label = "H = 94 nm", legend=:outertopright)
		plot!(angles_deg, Pex, linewidth = 2, label = "H = "*string(h)*" nm", legend=:outertopright)
	end
		#p = [1.0, 0.0, 365.0]
		#Pex = model_1c(constants, p)
		#plot!(angles_deg, Pex, linecolor = colors[4], linewidth = 2, label = "H = 365 nm")

	xlabel!("Incidence angle (degrees)")
	ylabel!("Intensity (A.U.)")
	display(plt)
	savefig(plt, "94_365_642nm.png")
	savefig(plt, "94_365_642nm.pdf")
	savefig(plt, "94_365_642nm.svg")

end
