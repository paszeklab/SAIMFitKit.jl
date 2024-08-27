using SAIMFitKit

#Construct a SAIM curve

opt = SAIMOptics()
opt.nB = 1.33;	                    #The refractive index of the ambient buffer / cytoplasm
opt.dOx = 1924.9;	                #The thickness of the SiO2 layer in units of nm
opt.λ_Ex_1 = 642.0;		        	#The wavelength of the microscope excitation laser #1 in units of nm
opt.nOx_1 = 1.4719;		            #The refractive index of SiO2
opt.nSi_1 = 3.8660 + 0.017933im;    #The complex refractive index of Si

x_range = [0.0, 45.0]
heights = [94.0, 200.0]
plot_curve_1c(opt, heights, x_range)