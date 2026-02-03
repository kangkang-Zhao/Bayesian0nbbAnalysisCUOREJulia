module LineshapeModule
export LineshapeModule
 
using UnROOT
using JSON
using LinearAlgebra
using StatsBase
using Distributions
using Random
using Printf

using Statistics
using Printf
using LinearAlgebra
using Distributions
using Random
using JLD2
using FileIO
using Plots
using StatsBase: Histogram

using ..readHistModule
 
"""
Implemention of LineshapeGauss3 in Julia.
    GetRandom、GetResponsePDF、SetPolyScaling ...
"""
# -----------------------
abstract type AbstractLineshape end

struct ParametersGauss3
    Qvalue::Float64
    Sigma::Float64
    Variance::Float64
    EnergyRatio::Float64  # ratio of peak amplitude: left-gauss / main-peak
    EnergyRatio2::Float64 # ratio of peak amplitude: right-gauss / main-peak
    Ratio::Float64        # ratio of energy left-shoulder / main-peak
    Ratio2::Float64       # ratio of energy right-shoulder / main-peak
end

struct BaselineParameters
    BkgSigma::Float64
    CalSigma::Float64
    BkgVariance::Float64
    CalVariance::Float64
end
# -----------------------
# Scaling helpers
# -----------------------
@enum ScalableParameter kQvalue kSigma
# Constants definition
const M_PI = π
const E208Tlpeak = 2614.511 


mutable struct LineshapeGauss3 <: AbstractLineshape
    fInputFileName::String
    fInputBkgFileName::String
    fInputCalFileName::String

    fParChannel::Dict{Int, ParametersGauss3}
    fParBaselineChannel::Dict{Int, BaselineParameters}

    fPolyScalingMap::Dict{ScalableParameter, Vector{Real}}
    fHistoScalingMap::Dict{ScalableParameter, Histogram}
    fTmpPolyScalingMap::Dict{ScalableParameter, Vector{Real}}
    fTmpPolyScalingInUse::Dict{ScalableParameter, Bool}

    fCovarianceMatrixPar::Dict{ScalableParameter, Matrix{Real}}
    fInverseCovarianceMatrixPar::Dict{ScalableParameter, Matrix{Real}}

    fResoModel::Int
    fDs::Int
    fExposure::Dict{Int, Float64}
end


# arguments with default values have to be put afer those without default values
function LineshapeGauss3(
    # peakshape file
    inputfilename::String = "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape/calibration/PeakShape_ds3819.root", 
    ds::Int = 3819)

    fParChannel = Dict{Int, ParametersGauss3}()
    fParBaselineChannel = Dict{Int, BaselineParameters}()
    fPolyScalingMap = Dict{ScalableParameter, Vector{Real}}()
    fHistoScalingMap = Dict{ScalableParameter, Histogram}()
    fTmpPolyScalingMap = Dict{ScalableParameter, Vector{Real}}()
    fTmpPolyScalingInUse = Dict{ScalableParameter, Bool}()
    fCovarianceMatrixPar = Dict{ScalableParameter, Matrix{Real}}()
    fInverseCovarianceMatrixPar = Dict{ScalableParameter, Matrix{Real}}()
    fExposure = Dict{Int, Float16}()

    # set initial default values
    fTmpPolyScalingInUse[kQvalue] = false
    fTmpPolyScalingInUse[kSigma] = false
    treename::String = "fitParamTree"
    LSfile = ROOTFile(inputfilename)
    if !haskey(LSfile, treename)
        error("Tree '$treename' not found in file")
    end
    LStree = LazyTree(LSfile, treename, 
        ["Dataset", 
        "Channel", 
        "Q_value", 
        "Sigma", 
        "SubPeakEnergyRatio", 
        "SubPeakEnergyRatio2", 
        "SubPeakRatio", 
        "SubPeakRatio2", 
        "Exposure"])

    ReallocDataset = LStree.Dataset
    ReallocChannel = LStree.Channel
    ReallocQvalue = LStree.Q_value
    ReallocSigma = LStree.Sigma
    ReallocEnergyRatio = LStree.SubPeakEnergyRatio
    ReallocEnergyRatio2 = LStree.SubPeakEnergyRatio2
    ReallocRatio = LStree.SubPeakRatio
    ReallocRatio2 = LStree.SubPeakRatio2
    ReallocExposure = LStree.Exposure
    Variance = ReallocSigma .^ 2
    
    for i in 1:lastindex(ReallocDataset)
        ch = ReallocChannel[i]
        if ReallocDataset[i] != ds
            continue
        end
        if ds == 3519 && ch == 755 
            @warn "Channel 755 in DS3519 is inactive, please check!"
            @warn "WARNING WARNING WARNING: manually removing ch 755 for ds 3519."
            @warn "If you still see this message, ask yourself if you should fix the code!!!"
            @warn "For further info, contact Alice, Gio, or Guido." 
            continue
        end
        if !haskey(fParChannel, ch)
            if ReallocExposure[i] <= 0.0
                @warn "Channel $(ch) in DS$(ds) has non-positive exposure $(ReallocExposure[i]), skipping."
                continue
            end

            par = ParametersGauss3(
                ReallocQvalue[i], 
                ReallocSigma[i], 
                Variance[i], 
                ReallocEnergyRatio[i], 
                ReallocEnergyRatio2[i], 
                ReallocRatio[i], 
                ReallocRatio2[i])

            fParChannel[ch] = par
            fExposure[ch] = ReallocExposure[i]
        else
            @warn "Duplicate entry for channel $(ch) in DS$(ds), skipping-position0."
        end
    end
    close(LSfile)
    fResoModel = 1

    println("LineshapeGauss3 initialized for DS$(ds) with $(length(fParChannel)) active channels.")
    # return LineshapeGauss3(
    return LineshapeGauss3(
        inputfilename, 
        "", 
        "", 
        fParChannel, 
        fParBaselineChannel, 
        fPolyScalingMap, 
        fHistoScalingMap, 
        fTmpPolyScalingMap, 
        fTmpPolyScalingInUse, 
        fCovarianceMatrixPar, 
        fInverseCovarianceMatrixPar, 
        fResoModel, 
        ds, 
        fExposure)
end


# Lineshape constructor with InputFileName, InputBkgSigma, InputCalSigma and Ds
function LineshapeGauss3(
    inputfilename::String,# = "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape/calibration/PeakShape_ds3819.root", 
    inputbkgsigma::String,# = "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/baseline_sigmas/bkg/ds3819_baseline.root", 
    inputcalsigma::String,# = "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/baseline_sigmas/calib/ds3819_baseline.root", 
    ds::Int)
    # Call the first constructor
    instance = LineshapeGauss3(inputfilename, ds)
    # Then read bkg/cal files and fill fParBaselineChannel
    instance.fInputBkgFileName = inputbkgsigma
    instance.fInputCalFileName = inputcalsigma
    bkgSigmaFile = ROOTFile(inputbkgsigma)
    calibSigmaFile = ROOTFile(inputcalsigma)
    treename_bkg::String = "info"
    if !haskey(bkgSigmaFile, treename_bkg)
        error("Tree '$treename_bkg' not found in file")
    end
    bkgTree = LazyTree(bkgSigmaFile, treename_bkg, ["channel", "dataset", "sigma"])
    ReallocDataset_bkg = bkgTree.dataset
    ReallocChannel_bkg = bkgTree.channel
    ReallocBkgSigma = bkgTree.sigma
    BkgVariance = ReallocBkgSigma .^ 2
    # read calibSigmaFile
    treename_calib::String = "info"
    if !haskey(calibSigmaFile, treename_calib)
        error("Tree '$treename_calib' not found in file")
    end
    calibTree = LazyTree(calibSigmaFile, treename_calib, ["channel", "dataset", "sigma"])
    ReallocDataset_calib = calibTree.dataset
    ReallocChannel_calib = calibTree.channel
    ReallocCalibSigma = calibTree.sigma
    CalibVariance = ReallocCalibSigma .^ 2

    if length(ReallocDataset_bkg) != length(ReallocDataset_calib)
        error("Baseline sigma files for bkg and calib have different lengths!")
    end 
    

    for i in 1:lastindex(ReallocDataset_bkg)
        ch = ReallocChannel_bkg[i]
        if !haskey(instance.fParChannel, ch)
            continue
        end
        if ReallocDataset_bkg[i] != ds
            continue
        end
        if ds == 3519 && ch == 755 
            @warn "WARNING WARNING WARNING: manually removing ch 755 for ds 3519."
            @warn "Baseline parameters found for inactive channel $(ch) in DS$(ds), skipping."
            continue
        end

        if ReallocBkgSigma[i] <= 0.0 || ReallocCalibSigma[i] <= 0.0
            @warn "Channel $(ch) in DS$(ds) has non-positive baseline sigma (bkg: $(ReallocBkgSigma[i]), cal: $(ReallocCalibSigma[i])), skipping."
            continue
        end
        sigma_peak::Float64 = instance.fParChannel[ch].Sigma
        if sigma_peak < ReallocCalibSigma[i]
            
            baselinepar = BaselineParameters(
                ReallocBkgSigma[i], 
                sigma_peak, 
                BkgVariance[i], 
                sigma_peak^2
                )
            @warn "For DsCh: $ds -- $ch  2615 sigma < baseline sigma: \n  $sigma_peak  < $(ReallocCalibSigma[i]), will set cal baseline = 2615 sigma"
        else
            baselinepar = BaselineParameters(
                ReallocBkgSigma[i], 
                ReallocCalibSigma[i], 
                BkgVariance[i], 
                CalibVariance[i])
        end
        instance.fParBaselineChannel[ch] = baselinepar
     
    end
    if length(instance.fParChannel) != length(instance.fParBaselineChannel)
        @warn "Warning: number of channels with fit parameters ($(length(instance.fParChannel))) does not match number with baseline parameters ($(length(instance.fParBaselineChannel)))."
    end

    close(bkgSigmaFile)
    close(calibSigmaFile)
    # Set the lineshape resolution model to 2 (baseline model)
    instance.fResoModel = 2
    return instance
end

# -----------------------
# Append ! to names of functions is a naming convention used for functions that mutate an argument.  
# When the function is supposed to be used in a loop, it is better to avoid unnecessary allocations.
# -----------------------

"""
    SetPolyScaling!(ls::LineshapeGauss3, name::Symbol, coeff::Vector{Float64})
"""
function SetPolyScaling!(ls::LineshapeGauss3, name::ScalableParameter, coeff::AbstractVector{<:Real})
    if isempty(coeff)
        @printf("WARNING: setting a polynomial scaling with no coefficients.\n")
        return false
    end   
    if haskey(ls.fPolyScalingMap, name)
        @printf("WARNING: overwriting polynomial scaling coefficients.\n")
    end   
    ls.fPolyScalingMap[name] = coeff
    # ??? Why set false here? fTmpPolyScaling and fPolyScaling are different options.
    ls.fTmpPolyScalingInUse[name] = false
    return true
end

"""
    SetPolyScaling!(ls::LineshapeGauss3, name::Symbol, inputrootfile::String, order::Int)
    To be implemented: if needed
从 ROOT 文件读取多项式拟合系数和协方差矩阵。
相当于 C++ 版的 LineshapeGauss3::SetPolyScaling。
"""
function SetPolyScaling!(ls::LineshapeGauss3, name::ScalableParameter, inputrootfile::String, order::Int)
    @warn "Function SetPolyScaling with input root file is to be implemented."
    return 0
end

function set_poly_scaling(ls::LineshapeGauss3, name::ScalableParameter, json_file::String, order::Int)

    # ----------------------------
    if name != kQvalue && name != kSigma
        error("Wrong parameter name passed. Scaling does not exist for it.")
    end
    # ----------------------------
    data = JSON.parsefile(json_file)
    key = name == kSigma ? "Sigma" : "Q"

    if !haskey(data, key)
        error("JSON does not contain scaling information for $key")
    end

    params = Float64.(data[key]["params"])

    if isempty(params)
        @warn "Setting a polynomial scaling with no coefficients. This may cause problems."
        return false
    end

    if length(params) < order + 1
        error("Requested polynomial order $order exceeds available parameters.")
    end

    coeff = params[1:order+1]

    # ----------------------------
    if haskey(ls.fPolyScalingMap, name)
        @warn "Overwriting polynomial scaling coefficients. Make sure this is REALLY what you want."
    end

    ls.fPolyScalingMap[name] = coeff
    ls.fTmpPolyScalingInUse[name] = false

    # ----------------------------
    cov_raw = data[key]["cov"]
    cov = reduce(vcat, (Float64.(row)' for row in cov_raw))

    ls.fCovarianceMatrixPar[name] = cov

    return true
end


# Overload: set histo scaling without hist -> set identity
function SetHistoScaling(ls::LineshapeGauss3, name::ScalableParameter)
    ls.fPolyScalingMap[name] = [0.0]
    ls.fTmpPolyScalingInUse[name] = false
    return true
end

function SetHistoScaling(ls::LineshapeGauss3, name::ScalableParameter, hist::Histogram)
    ls.fHistoScalingMap[name] = hist
    # take mode as single coefficient as in C++
    idx = argmax(hist.weights)
    mode = (hist.edges[1][idx] + hist.edges[1][idx+1]) / 2
    ls.fPolyScalingMap[name] = [mode]
    ls.fTmpPolyScalingInUse[name] = false
    return true
end


function GetPolyScaling!(ls::LineshapeGauss3, name::ScalableParameter)
    if !haskey(ls.fPolyScalingMap, name)
        @warn "GetHistoScaling called for non-existing histogram scaling for $(name)."
        return nothing
    end
    return ls.fPolyScalingMap[name]
end

function SetTmpPolyScaling(
    ls::LineshapeGauss3,
    name::ScalableParameter,
    # coeff::AbstractVector{<:Real}
    coeff::Vector{<:Real}
)
    if !haskey(ls.fPolyScalingMap, name)
        error("SetTmpPolyScaling: base poly not set for $name")
    end
    if length(coeff) != length(ls.fPolyScalingMap[name])
        error("SetTmpPolyScalingChecked: size mismatch")
    end
    ls.fTmpPolyScalingInUse[name] = true
    ls.fTmpPolyScalingMap[name] = copy(coeff)
end

function ResetTmpPolyScaling(ls::LineshapeGauss3, name::ScalableParameter)
    ls.fTmpPolyScalingInUse[name] = false
end

# -----------------------
# Response function (3-gauss) and PDF / random
# -----------------------
function ResponseFunction(E::Real, p::AbstractVector{<:Real})
    r1 = p[1]
    r2 = p[2]
    qValue = p[3]
    σ = p[4]
    ER1 = p[5]
    ER2 = p[6]
    # Using promote_type to ensure consistent type for calculations
    T = promote_type(typeof(E), eltype(p))

    invdenom = one(T) / (2 * σ^2)
    norm = (one(T) + r1 + r2) * sqrt(2 * T(pi)) * σ

    val =
        (r1 * exp(-((E - ER1*qValue)^2) * invdenom) +
         r2 * exp(-((E - ER2*qValue)^2) * invdenom) +
         exp(-((E - qValue)^2) * invdenom)) / norm

    return val
end

function BuildScalingPoly(ls::LineshapeGauss3, name::ScalableParameter, energy::Real)
    if ls.fTmpPolyScalingInUse[name] && haskey(ls.fTmpPolyScalingMap, name)
        myScaling = ls.fTmpPolyScalingMap[name]
    elseif haskey(ls.fPolyScalingMap, name)
        myScaling = ls.fPolyScalingMap[name]
    else
        return 0.0
    end
    
    output = 0.0
    maxOrder = length(myScaling)
    # reverse loop
    for order in maxOrder:-1:1
        output *= energy
        output += myScaling[order]
    end
    
    return output #maxbin along with energy
end

function GetQValueScaling(ls::LineshapeGauss3, value::Float64, energy::Real, scalePoly::Real)
    return value * energy / E208Tlpeak + scalePoly
end

    
function GetSigmaValueScaling(ls::LineshapeGauss3, chan::Int, value::Real, resoModel::Int, scalePoly::Real)
    if resoModel == 1
        return scalePoly * value
    elseif resoModel == 2
        if !haskey(ls.fParBaselineChannel, chan)
            error("Channel $chan not found in baseline parameters")
        end
        
        bpar = ls.fParBaselineChannel[chan]
        bkgVar = bpar.BkgVariance
        calVar = bpar.CalVariance
        
        return sqrt(bkgVar + scalePoly * (value^2 - calVar))
    else
        error("Unknown resolution model: $resoModel")
    end
end

function GetParameterPolyScaling(
    ls::LineshapeGauss3,
    name::ScalableParameter,
    chan::Int,
    value::Real,
    energy::Real,
    resoModel::Int
)
    if !haskey(ls.fPolyScalingMap, name) && !(ls.fTmpPolyScalingInUse[name] && haskey(ls.fTmpPolyScalingMap, name))
        if name == kSigma
            return value
        elseif name == kQvalue
            return value * energy / E208Tlpeak
        end
    end
    
    scalePoly = BuildScalingPoly(ls, name, energy)
    
    if name == kQvalue
        return GetQValueScaling(ls, value, energy, scalePoly)
    elseif name == kSigma
        return GetSigmaValueScaling(ls, chan, value, resoModel, scalePoly)
    end
    
    return scalePoly
end


function GetResponsePDF(
    ls::LineshapeGauss3,
    channel::Int,
    energy::Real,
    expectedEnergy::Real
)
    if !haskey(ls.fParChannel, channel)
        @printf(
            "ERROR: a call to GetResponsePDF for a non-active ch-ds was done. Channel = %d, Dataset = %d\n",
            channel, ls.fDs
        )
        return zero(energy)
    end

    params = ls.fParChannel[channel]

    T = promote_type(typeof(energy), typeof(expectedEnergy))

    p = Vector{T}(undef, 6)

    p[1] = params.Ratio
    p[2] = params.Ratio2
    p[3] = GetParameterPolyScaling(
        ls, kQvalue, channel, params.Qvalue, expectedEnergy, ls.fResoModel
    )
    p[4] = GetParameterPolyScaling(
        ls, kSigma, channel, params.Sigma, expectedEnergy, ls.fResoModel
    )
    p[5] = params.EnergyRatio
    p[6] = params.EnergyRatio2

    return ResponseFunction(energy, p)
end

# Sampling from the response (maps to TF1::GetRandom)
function GetRandom(ls::LineshapeGauss3, channel::Int32, expectedEnergy::Float64, min::Float64, max::Float64)
    if !haskey(ls.fParChannel, channel)
        @warn "GetRandom on inactive channel $channel"
        return 0.0
    end
    ppars = ls.fParChannel[channel]
    # evaluate scaled Q and S
    Q = GetParameterPolyScaling(ls, kQvalue, channel, ppars.Qvalue, expectedEnergy, ls.fResoModel)
    S = GetParameterPolyScaling(ls, kSigma, channel, ppars.Sigma, expectedEnergy, ls.fResoModel)
    # components: three Gaussians centered at ER*Q, ER2*Q, Q with same sigma σ, weights proportional to (r1, r2, 1)
    weights = [ppars.Ratio, ppars.Ratio2, 1.0]
    total = sum(weights)
    probs = weights ./ total
    # choose one random component according to weights
    comp = rand(Categorical(probs)) # gets random index value 1,2,3
    means = [ppars.EnergyRatio * Q, ppars.EnergyRatio2 * Q, Q]
    # sample from truncated normal in [min,max]
    μ = means[comp]
    σ = S
    # rejection sampling simple
    for iter in 1:10000
        x = rand(Normal(μ,σ))
        if x>=min && x<=max
            return x
        end
    end
    @warn "GetRandom: rejection sampling failed after many tries; returning μ"
    return μ
end


function GetExposure(ls::LineshapeGauss3)
    return ls.fExposure  
end

function GetTotExposure(ls::LineshapeGauss3)
    exposure = 0.0
    for exp_val in values(ls.fExposure)
        exposure += exp_val
    end 
    return exposure
end

function GetDs(ls::LineshapeGauss3)
    return ls.fDs
end

function GetIsActiveChannel(ls::LineshapeGauss3, channel::Int32)
    return haskey(ls.fParChannel, channel)
end

function GetPolyScaling(ls::LineshapeGauss3, name::ScalableParameter)
    return get(ls.fPolyScalingMap, name, Real[])
end
# -----------------------
# Export a small API
# -----------------------
export LineshapeGauss3, ScalableParameter, kQvalue, kSigma
export GetResponsePDF, GetRandom, SetPolyScaling, GetPolyScaling
export SetTmpPolyScaling, ResetTmpPolyScaling
export GetExposure, GetTotExposure, GetDs, GetIsActiveChannel
export set_poly_scaling
export SetHistoScaling
end # module

