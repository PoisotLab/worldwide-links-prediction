# Activate project environment
import Pkg; Pkg.activate(".")

Pkg.instantiate()

# Load required packages
import CSV
using DataFrames
using Distributions
using Plots
using ProgressMeter
using Statistics
using Turing

using ArchGDAL
using Mangal
using SimpleSDMLayers

# Load scripts
include("src/1_import_mangal_data.jl")
include("src/2_get_measures.jl")
include("src/3_make_figures.jl")
