using SAIMFitKit

#-------------- OPTICAL PARAMETERS----------------------
nB = 1.33;			#The refractive index of the ambient media / cytoplasm
dOx = 1900.1;		#The thickness of the SiO2 layer in units of nm

λ_Ex_1 = 488.0; 		        #The wavelength of excitation laser #1 in units of nm
nOx_1 = 1.4775;		            #The refractive index of SiO2 at excitation wavelength #1
nSi_1 = 4.3707 + 0.080068im; 	#The complex refractive index of Si at excitation wavelength #1; #4.3 + 0.073im; 0.080068

λ_Ex_2 = 561.0;		        	#The wavelength of excitation laser #3 in units of nm
nOx_2 = 1.4595;		            #The refractive index of SiO2 at excitation wavelength #3
nSi_2 = 4.0506 + 0.033698im;    #The complex refractive index of Si at excitation wavelength #3; 

λ_Ex_3 = 642.0; 		        #The wavelength of excitation laser #2 in units of nm
nOx_3 = 1.4719;		            #The refractive index of SiO2 at excitation wavelength #2
nSi_3 = 3.8660 + 0.017933im; 	#The complex refractive index of Si at excitation wavelength #2; #4.3 + 0.073im;k = 0.017933

opt1 = SAIMOptics(nB, nOx_1, nSi_1, dOx, λ_Ex_1)
opt2 = SAIMOptics(nB, nOx_2, nSi_2, dOx, λ_Ex_2)
opt3 = SAIMOptics(nB, nOx_3, nSi_3, dOx, λ_Ex_3)

#-------------- FIT PARAMETERS----------------------
p0 = [0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 75.0]		            #Initial guesses for parameters [A1, B1, A2, B2, H]
lb = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]					#Lower bounds for [A1, B1, A2, B2, H]
ub = [10000., 10000., 10000., 10000., 10000., 10000., 450.]		#Upper bounds for [A1, B1, A2, B2, H]
grid = false						                #Conduct global grid search when true
step = 40.0						                    #Step size for global grid search
show = true                                        #Display generated heatmaps when true

#-------------- IMAGE STACKS ---------------------
angles = range(5.0, length=32, stop=43.75)	#Incidece angle in degrees for each image frame	
files = ["488_560_640"]				#Names of image stacks; multiple experiments can be listed to run in batch
path = "D:\\Google Drive\\Admin\\2017-2023 Cornell\\Research\\3_Raw_Data\\210304_SP_SAIM_Three_Color\\D6_MemGlow_488_590_640\\Test_4_2" #Path to image files

#-------------- CONDUCT FITS ---------------------
for f in files
    fit_SAIM(path, f, opt1, opt2, opt3, angles, p0, lb, ub, disp=show, glb=grid)
end
