#Store details of the optical system 
struct SAIMOptics
	n_Buffer::Float64		        #The refractive index of the ambient media / cytoplasm
	n_Oxide::Float64		        #The refractive index of SiO2
	n_Silicon::Complex{Float64}	    #The complex refractive index of Si
	d_Oxide::Float64		        #The thickness of the SiO2 layer in units of nm
	lambda_Ex::Float64		        #The wavelength of the microscope excitation laser #1 in units of nm
end


