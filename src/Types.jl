#Store details of the optical system 
mutable struct SAIMOptics
	#Required parameters for standard SAIM experiments with single excitation
	nB::Float64		        	#The refractive index of the ambient buffer / cytoplasm
	dOx::Float64		        #The thickness of the SiO2 layer in units of nm
	λ_Ex_1::Float64		        #The wavelength of excitation light #1 in units of nm
	nOx_1::Float64		        #The refractive index of SiO2 for excitation #1
	nSi_1::Complex{Float64}	    #The complex refractive index of Si for excitation #1
	
	#Additional parameters for two-color and three-color SAIM; these parameters are
	#not used for single excitation experiments
	λ_Ex_2::Float64		        #The wavelength of excitation light #2 in units of nm
	nOx_2::Float64		        #The refractive index of SiO2 for excitation #2
	nSi_2::Complex{Float64}	    #The complex refractive index of Si for excitation #2
	λ_Ex_3::Float64		        #The wavelength of excitation light #3 in units of nm
	nOx_3::Float64		        #The refractive index of SiO2 for excitation #3
	nSi_3::Complex{Float64}	    #The complex refractive index of Si for excitation #3

	#Constructor for the structure with all variaable initialized to zero
	function SAIMOptics()
		this = new(0.0, 0.0, 0.0, 0.0, 0.0+0.0im, 0.0, 0.0, 0.0+0.0im, 0.0, 0.0, 0.0+0.0im)
	end
end


