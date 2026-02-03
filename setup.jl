#!/usr/bin/env julia
using Pkg

println("Setting up project environment...")
# 激活或创建项目
Pkg.activate(".")

# 添加所有依赖
Pkg.add([
    "DataFrames",
    "CSV", 
    "Plots",
    "Statistics",
    "Printf",
    "LinearAlgebra",
    "DelimitedFiles",   
    "StatsBase",
    "UnROOT",
    "Distributions",
    "Random",
    "JLD2",
    "FileIO",
    "JSON",
    "IntervalSets",
    "OrderedCollections",
    "EmpiricalDistributions",
    "ValueShapes",
    "BAT" ,
    "DensityInterface",
    "StatsPlots",
    "AutoDiffOperators",
    "Cuba",
    "AdvancedHMC",
    "ForwardDiff",
    "Random123",
    "Optim",
    "HDF5"
])

# 生成/更新 Manifest.toml
Pkg.resolve()

println("Environment setup complete!")
println("Manifest.toml has been generated/updated.")

