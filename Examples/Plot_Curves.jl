using SAIMFitKit
using Plots

#Example illustrates how to contruct and plot a one-color SAIM curve or curves
#for a single emitter height or an array of heights (94, and 200 nm in this example)
#Uses the utility function plot_curve_1c
#Construct a SAIM curve

opt = SAIMOptics()
opt.nB = 1.33;	                    #The refractive index of the ambient buffer / cytoplasm
opt.dOx = 2049.1;	                #The thickness of the SiO2 layer in units of nm
opt.λ_Ex_1 = 488.0;		        	#The wavelength of the microscope excitation laser #1 in units of nm
opt.nOx_1 = 1.4719;		            #The refractive index of SiO2
opt.nSi_1 = 3.8660 + 0.017933im;    #The complex refractive index of Si

x_range = [0.0, 35.0]
heights = [0.0, 20.0, 40.0, 60.0]
plt1 = plot_curve_1c(opt, heights, x_range)
savefig(plt1, "SAIM_Curve.png")
savefig(plt1, "SAIM_Curve.pdf")
savefig(plt1, "SAIM_Curve.svg")
display(plt1)

ifile = "/Users/matthewpaszek/Documents/GitHub/SAIMFitKit.jl/Examples/Beads/240913_Green_1.tif"
jfile = "/Users/matthewpaszek/Documents/GitHub/SAIMFitKit.jl/Examples/Beads/output/240913_Green_1/240913_Green_1_results.jld"
plt2 = plot_fit(ifile, jfile, 40, 40)

savefig(plt2, "example_fit1.png")
savefig(plt2, "example_fit1.pdf")
savefig(plt2, "example_fit1.svg")
display(plt2)