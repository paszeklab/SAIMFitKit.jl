module SAIMFitKit

    export  fit_SAIM,
            fit_SAIM_g,
            fit_SAIM2c,
            fit_SAIM2c_g,
            test_working,
            #Types
            SAIMOptics

    # using statement for external packages -
    using Plots
    using LsqFit
    using Images, FileIO
    using JLD
    using Statistics
    using CSV
    using DataFrames

    # include package codes -
    include("Types.jl")
    include("OneColor.jl")
    include("TwoColor.jl")
    include("Utilities.jl")

end # module
