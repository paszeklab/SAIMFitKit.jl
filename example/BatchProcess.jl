using SAIMFitKit

#-------------- OPTICAL PARAMETERS----------------------

opt = SAIMOptics()
opt.nB = 1.33;	                    #The refractive index of the ambient buffer / cytoplasm
opt.dOx = 1924.9;	                #The thickness of the SiO2 layer in units of nm
opt.λ_Ex_1 = 642.0;		        	#The wavelength of the microscope excitation laser #1 in units of nm
opt.nOx_1 = 1.4719;		            #The refractive index of SiO2
opt.nSi_1 = 3.8660 + 0.017933im;    #The complex refractive index of Si

#-------------- FIT PARAMETERS----------------------
p0 = [0.5, 1.0, 75.0]		    #Initial guesses for parameters [A, B, H]
lb = [0.0, 0.0, 0.0]			#Lower bounds for [A, B, H]
ub = [10000., 10000., 200.]		#Upper bounds for [A, B, H]
g = true					    #Conduct global grid search when true
step = 40.0						#Step size for global grid search
d = true                        #Display generated heatmaps when true

#-------------- IMAGE STACKS ---------------------
angles = range(5.0, length=32, stop=43.75)	    #Incidece angle in degrees for each image frame	
files = ["21TR_5c"] #, "21TR_5c"]               #Names of image stacks; multiple experiments can be listed to run in batch	
path = "C:\\Users\\matth\\Documents\\Julia Scripts\\SAIMFitKit\\example\\TestImages" #Path to image files

#-------------- CONDUCT FITS ---------------------
for f in files
    fit_SAIM(path, f, opt, angles, p0, lb, ub, glb=g, disp=d)
end

