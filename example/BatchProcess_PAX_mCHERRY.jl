using SAIMFitKit

#Example for processing an example paxillin-mCherry data set

#-------------- OPTICAL PARAMETERS----------------------

opt = SAIMOptics()
opt.nB = 1.33;	                    #The refractive index of the ambient buffer / cytoplasm
opt.dOx = 2040;#1924.9	                #The thickness of the SiO2 layer in units of nm
opt.λ_Ex_1 = 561.0; #642		        	#The wavelength of the microscope excitation laser #1 in units of nm
opt.nOx_1 = 1.4719;		            #The refractive index of SiO2
opt.nSi_1 = 3.8660 + 0.017933im;    #The complex refractive index of Si

#-------------- FIT PARAMETERS----------------------
p0 = [100, 500.0, 80.0]		    #Initial guesses for parameters [A, B, H]
lb = [0.0, 0.0, 0.0]			#Lower bounds for [A, B, H]
ub = [10000., 10000., 200.]		#Upper bounds for [A, B, H]

#-------------- IMAGE STACKS ---------------------
angles_deg = range(5.0, length=32, stop=40)	    #Incidece angle in degrees for each image frame	
angles = angles_deg*2*pi/360  	                    #Incidece angle in radians for each image frame	
files = ["PAX_1"]                             #Names of image stacks; multiple experiments can be listed to run in batch	
path = "/Users/matthewpaszek/Documents/GitHub/SAIMFitKit.jl/example/TestImages" #Path to image files

#-------------- CONDUCT FITS ---------------------
for f in files
    fit_SAIM(path, f, opt, angles, p0, lb, ub, glb=false, step=20.0, show=true)
end
