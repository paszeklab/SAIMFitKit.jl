using SAIMFitKit

#-------------- OPTICAL PARAMETERS----------------------
nB = 1.33;			#The refractive index of the ambient media / cytoplasm
dOx = 1924.9;		#The thickness of the SiO2 layer in units of nm

λ_Ex_1 = 642.0; 	#The wavelength of excitation laser #2 in units of nm
nOx_1 = 1.4719;		#The refractive index of SiO2 at excitation wavelength #2
nSi_1 = 3.8660 + 0.017933im; 	#The complex refractive index of Si at excitation wavelength #2; #4.3 + 0.073im;k = 0.017933

opt1 = SAIMOptics(nB, nOx_1, nSi_1, dOx, λ_Ex_1)

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
    fit_SAIM(path, f, opt1, angles, p0, lb, ub, glb=g, disp=d)
end

