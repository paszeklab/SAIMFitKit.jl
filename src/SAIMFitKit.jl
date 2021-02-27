module SAIMFitKit

    export  fit_SAIM,
            plot_fit,
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
    include("DataFit.jl")
    include("OneColor.jl")
    include("TwoColor.jl")
    include("ThreeColor.jl")
    include("Utilities.jl")

end # module
