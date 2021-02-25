#Plot a SAIM curve
function plot_curve()

	#-------------- OPTICAL PARAMETERS----------------------
	nB = 1.33;			#The refractive index of the ambient media / cytoplasm
	dOx = 1924.9;		#The thickness of the SiO2 layer in units of nm
	λ_Ex_1 = 642.0; 		#The wavelength of excitation laser #1 in units of nm
	nOx_1 = 1.463;		#The refractive index of SiO2 at excitation wavelength #1
	nSi_1 = 4.3638; 	#The complex refractive index of Si at excitation wavelength #1; #4.3 + 0.073im;

	opt1 = Optical(nB, nOx_1, nSi_1, dOx, λ_Ex_1)

	#-------------- ANGLES  ---------------------
	angles_deg = range(5.0, length=100, stop=43.75)	#Incendece angle in degree for each image frame
	angles_rad = angles_deg*DTR						#Incendece angle in radians for each image frame
	num_angles = length(angles_rad)

	#-------------- SIMULATION----------------------
	colors = palette(:tab10)

	constants = Calculate_Constants1c(opt1, angles_rad)

	Pex = Pex_from_Constants(94.0, constants)
	plt = plot(angles_deg, Pex, linecolor = colors[1],
		linewidth = 2, label = "H = 94 nm", legend=:outertopright)

	Pex = Pex_from_Constants(365.0, constants)
	plot!(angles_deg, Pex, linecolor = colors[4],
		linewidth = 2, label = "H = 365 nm")

	xlabel!("Incidence angle (degrees)")
	ylabel!("Intensity (A.U.)")
	display(plt)
	savefig(plt, "94_365_642nm.png")
	savefig(plt, "94_365_642nm.pdf")
	savefig(plt, "94_365_642nm.svg")

end

#Function to show example fits for a given pixel in a SAIM image sequence
function plot_fit(image_file, JLD_file, row, col)

	println(image_file)
	println(JLD_file)

	#Load the constants from the JLD_file
	constants = load(JLD_file, "constants")
	angles = load(JLD_file, "angles")
	angles = RTD*angles		#convert from radians to degrees

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
	yfit = model(constants, p)
    colors = palette(:tab10)
	plt = plot(angles, ydata, label="data", seriestype = :scatter, markercolor = colors[1])
	plot!(angles, yfit, label="fit", linecolor = colors[4], linewidth = 1.5)
    xlabel!("Incidence angle (degrees)")
	ylabel!("Intensity (A.U.)")
	
    savefig(plt, "example_fit.png")
	savefig(plt, "example_fit.pdf")
	savefig(plt, "example_fit.svg")
    display(plt)
end