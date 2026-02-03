# Dataset.jl - Dataset module for handling experimental datasets
module DatasetModule
export DatasetModule

using Printf
using LinearAlgebra
using Statistics
using DelimitedFiles
using UnROOT
using Distributions
using Random
using CSV
using DataFrames
using StatsBase

# import Lineshape module
using ..readHistModule
using ..PhysicalConstants
using ..LineshapeModule

# ======================================================
# AbstractDataset（corresponds to C++: class Dataset）
# ======================================================
abstract type AbstractDataset end

# Imutable Event struct
struct Event
    run::UInt32
    multiplicity::Int64
    type::Int32  # only used for toy-MC. 0=bkg, 1=0nbb, 2=60Co
    energy::Vector{Float64}
    channel::Vector{Int32}
end

# Event constructors
Event() = Event(0, -1, -1, Float64[], Int32[])

function Event_init(energy::Vector{Float64}, channel::Vector{Int32}, mult::Int64=-1, run::UInt32=0, event_type::Int32=-1)
    # mult == length(energy) || error("mult != energy.size()")
    if mult != length(energy)
        error("FATAL Error in Event constructor: mult != energy.size()")
    end
    if length(energy) != length(channel)
        error("FATAL Error in Event constructor: energy.size() != channel.size()")
    end
    return Event(run, mult, event_type, energy, channel)
end
# Event initializer functions
function Event_init(energy::Float64, channel::Int32, event_type::Int32=-1)
    return Event(0, -1, event_type, [energy], [channel])
end

function get_energy(ev::Event, index::Int=1)
    if index > ev.multiplicity
        error("Error in get_energy: $index > multiplicity ($(ev.multiplicity))")
    end
    return ev.energy[index]
end

function get_channel(ev::Event, index::Int=1)
    if index > ev.multiplicity
        error("Error in get_channel: $index > multiplicity ($(ev.multiplicity))")
    end
    return ev.channel[index]
end


mutable struct Dataset <: AbstractDataset
    # basic info
    ds::Int
    filename::String
    lineshape::LineshapeModule.LineshapeGauss3
    
    # optins/parameters
    cut_efficiency::Float64
    cut_efficiency_err::Float64

    containment_efficiency::Float64
    containment_efficiency_err::Float64
    delta_t::Float64

    q_scaling::Float64
    q_scaling_err::Float64
    sigma_scaling::Float64
    sigma_scaling_err::Float64


    emin::Float64
    emax::Float64
    
    n_events::Int 
    # options
    exposure_override::Bool
    th1_efficiency::Bool
    account_pulser_xtalk::Bool

    file_efficiency::String
    file_triggereff::String
    file_exposure::String
    livetime_filename::String
    
    eff_prior::Histogram  # placeholder for TH1D*

    histo_ls_scaling::Bool
    bias_ls_scaling::Histogram  
    reso_ls_scaling::Histogram  
    
    # data containers
    exposure_channel::Dict{Int, Float64}
    total_efficiency_channel::Dict{Int, Float64}
    total_efficiency_err_channel::Dict{Int, Float64}
    trigger_efficiency_channel::Dict{Int, Float64}
    events_channel::Dict{Int, Vector{Event}}
    events::Vector{Event}
    
    # other
    co_reduction::Float64
          
end 

# Dataset initializer/constructor
# Arguments with defalult values are for DS3819 Super Reduced Background
function Dataset(
    fLineshape::LineshapeModule.LineshapeGauss3,
    fDs::Int=3819,
    fFilename::String = "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/super_reduced/background/unblinded/SuperReduced_Background_ds3819.root",
    
    fCut_efficiency::Float64=0.947,
    fCut_efficiency_err::Float64=0.0075,
    fContainment_efficiency::Float64=0.88345,
    fContainment_efficiency_err::Float64=0.00085,
    fDelta_t::Float64=237.0,

    fQ_scaling::Float64=1.0,
    fQ_scaling_err::Float64=0.0,
    fSigma_scaling::Float64=1.0,
    fSigma_scaling_err::Float64=0.0,

    fEmin::Float64=2465.0,
    fEmax::Float64=2575.0,

    fTH1_efficiency::Bool=false,
    # fHisto_ls_scaling::Bool=false,
    fExposure_override::Bool=false,
    fAccount_pulser_xtalk::Bool=false,

    fFile_efficiency::String="/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/CombinedEfficiency/CombinedEfficiencies.root",
    fFile_triggereff::String="",
    fFile_exposure::String="/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/Exposures/Exposures_ds3819.txt",
    fLivetime_filename::String="",

    fHisto_ls_scaling::Bool = false,
    fFile_HistoScaling::String=""
)
    fN_events::Int = -1
    fExposure_channel = Dict{Int, Float64}()
    fTotal_efficiency_channel = Dict{Int, Float64}()
    fTotal_efficiency_err_channel = Dict{Int, Float64}()
    fTrigger_efficiency_channel = Dict{Int, Float64}()
    fEvents_channel = Dict{Int, Vector{Event}}()
    fEvents = Vector{Event}()

    if fTH1_efficiency
        if fFile_efficiency == ""
            error("fTH1_efficiency is set to true but no efficiency file provided")
        elseif isfile(fFile_efficiency)
            println("File exists.")
        else
            error("fTH1_efficiency is set to true but efficiency file does not exist: $fFile_efficiency")
        end
        # read TH1D* from fFile_TH1_efficiency
        histname = "CombinedEff_DS$(fDs)"
        println("Reading histogram '$histname' from $fFile_efficiency ...")
        edges, vals = read_root_hist(fFile_efficiency, histname)
        
        # println(pwd())
        # outfile = "$histname.png"
        # println("Plotting and saving to $outfile ...")

        h = Histogram(edges, vals)
        fEff_prior = normalize_histogram(h)

    else
        fEff_prior = Histogram(Float64[1,2], Float64[1])
    end 

    fCo_reduction = exp(-fDelta_t / PhysicalConstants.Tau60Cobalt)

    if fHisto_ls_scaling
        if fFile_HistoScaling == ""
            error("fHisto_ls_scaling is set to true but no LS scaling file provided")
        elseif isfile(fFile_HistoScaling)
            println("LS scaling file exists.")
        else
            error("fHisto_ls_scaling is set to true but LS scaling file does not exist: $fFile_HistoScaling")
        end
        println("Reading LS scaling histograms from $fFile_HistoScaling ...")
        percentile::Float64 = 0.997
        println("The histogram priors may be a bit odd out at the tails, so truncate them at percentile of $percentile")


        # read TH1D* from fFile_HistoScaling
        histname_bias = "bias_Qbb_ds$(fDs)"
        histname_reso = "reso_Qbb_ds$(fDs)"
        println("Reading histogram '$histname_bias' ...")
        edges_bias, vals_bias = read_root_hist(fFile_HistoScaling, histname_bias)
        histo_bias = truncate_histogram(Histogram(edges_bias, vals_bias),percentile)
        fBias_ls_scaling = normalize_histogram(histo_bias)
        println("Reading histogram '$histname_reso' ...")
        edges_reso, vals_reso = read_root_hist(fFile_HistoScaling, histname_reso)
        histo_reso = truncate_histogram(Histogram(edges_reso, vals_reso),percentile)
        fReso_ls_scaling = normalize_histogram(histo_reso)
    else
        println("Histogram LS scaling not used.")
        fBias_ls_scaling = Histogram(Float64[1,2], Float64[1])
        fReso_ls_scaling = Histogram(Float64[1,2], Float64[1])
    end

    dataset = Dataset(   
        fDs, 
        fFilename,
        fLineshape,

        fCut_efficiency,
        fCut_efficiency_err,
        fContainment_efficiency,
        fContainment_efficiency_err,
        fDelta_t,

        fQ_scaling,
        fQ_scaling_err,
        fSigma_scaling,
        fSigma_scaling_err,
        fEmin,
        fEmax,
        
        fN_events,

        fExposure_override,
        fTH1_efficiency,
        fAccount_pulser_xtalk,

        fFile_efficiency,
        fFile_triggereff,
        fFile_exposure,
        fLivetime_filename,

        fEff_prior,

        fHisto_ls_scaling,
        fBias_ls_scaling,
        fReso_ls_scaling,

        fExposure_channel,
        fTotal_efficiency_channel,
        fTotal_efficiency_err_channel,
        fTrigger_efficiency_channel,
        fEvents_channel,
        fEvents,
        fCo_reduction)
    read_data(dataset)
    set_all_efficiencies(dataset)
    # set_all_exposures(dataset, 1.0)
    dataset.exposure_channel = LineshapeModule.GetExposure(dataset.lineshape)
    set_all_active_exposures_from_file(dataset)
    # when take into account the Pulser Cross Talk
    if fAccount_pulser_xtalk
        set_livetime_fraction_from_file(dataset)
    end

    @info "Dataset initialized: DS=$(dataset.ds), Events=$(dataset.n_events), Exposure sum=$(sum(values(dataset.exposure_channel))) kg·yr"
    return dataset

end



# Load data from ROOT file into Dataset
function read_data(dataset::Dataset)
    ds = dataset.ds
    Emin = dataset.emin
    Emax = dataset.emax
    println("Opening file $(dataset.filename)")
    
    file = ROOTFile(dataset.filename)
    # check the existence of the certain Tree
    treename = "qtree_bkg"
    if !haskey(file, treename)
        # print all the content in the file
        println("Keys in root file:")
        for key in keys(file)
            println("  ", key)
        end
        error("Tree '$treename' not found in file")
    end

    readtree = LazyTree(file, treename, 
        ["Run",
        "Dataset",
        "Channel",
        "Energy",
        "Multiplicity",
        "GoodForAnalysisV",
        "AnalysisBaseCutV",
        "PCACutV"])

    ReallocRun              = readtree.Run
    ReallocDataset          = readtree.Dataset
    ReallocChannel          = readtree.Channel
    ReallocEnergy           = readtree.Energy
    ReallocMulti            = readtree.Multiplicity
    ReallocGoodForAnalysisV = readtree.GoodForAnalysisV
    ReallocAnalysisBaseCutV = readtree.AnalysisBaseCutV
    ReallocPCACutV          = readtree.PCACutV

    nEntries = lastindex(ReallocRun)
    println("Number of entries in tree: $nEntries")
    for i in 1:nEntries
        # Skip non-matching datasets
        if ReallocDataset[i] != ds
            continue
        end

        # check if event quality is good
        if ReallocMulti[i] != 1 || !ReallocGoodForAnalysisV[i][1] || !ReallocAnalysisBaseCutV[i][1] || !ReallocPCACutV[i][1]
            continue
        end
        
        energy_val = ReallocEnergy[i]
        channel_val = ReallocChannel[i]
        # Chenk energy window
        if energy_val > Emin && energy_val < Emax
            # check if channel is active (non-zero exposure) in lineshape
            if !LineshapeModule.GetIsActiveChannel(dataset.lineshape, channel_val)
                continue
            end
            event = Event_init([energy_val], [channel_val], 1, ReallocRun[i], Int32(0))
            push!(dataset.events, event)
            
            if !haskey(dataset.events_channel, channel_val)
                dataset.events_channel[channel_val] = [event]
            else
                push!(dataset.events_channel[channel_val], event)
            end
        end
    end
    dataset.n_events = length(dataset.events)
    close(file)
end


function set_all_efficiencies(dataset::Dataset)
    efficiency_flag = false
    # Set trigger efficiency by-channel from input ASCII file
    if !isempty(dataset.file_triggereff)
        @warn "ASCII file supplied for trigger efficiency by-channel. Read program not implemented yet."
        # try
        
        #     # read ASCII file, to be implemented!!!
        #     data = readdlm(dataset.file_triggereff, skipstart=0)
        #     for row in eachrow(data)
        #         if length(row) >= 2
        #             chan = Int(row[1])
        #             effi = Float64(row[2])
                    
        #             if effi > 1.0 || effi < 0.0 || chan < 1 || chan > 988
        #                 error("Channel $chan has efficiency $effi")
        #             end
                    
        #             dataset.trigger_efficiency_channel[chan] = effi
        #             prinln("Channel $chan trigger efficiency set to $effi")
        #         end
        #     end
            
        # catch e
        #     error("Could not read input ASCII efficiency file $(dataset.file_efficiency): $e")
        # end
    else
        @warn "No input ASCII for trigger efficiency by-channel. Setting all to 1."
        efficiency_flag = true
    end
    
    # compute total efficiency and its error
    for ch in 1:988
        if efficiency_flag
            dataset.trigger_efficiency_channel[ch] = 1.0
        end
        dataset.total_efficiency_channel[ch] = dataset.containment_efficiency *
                                              dataset.cut_efficiency *
                                              dataset.trigger_efficiency_channel[ch]
        
        dataset.total_efficiency_err_channel[ch] = dataset.trigger_efficiency_channel[ch] *
                                                  sqrt((dataset.containment_efficiency_err * dataset.cut_efficiency)^2 +
                                                       (dataset.containment_efficiency * dataset.cut_efficiency_err)^2)
    end
end


# Get the energy and channel of an event from a Dataset
function get_energy(dataset::Dataset, i_event::Int, index::Int=1)
    return get_energy(dataset.events[i_event], index)
end

function get_channel(dataset::Dataset, i_event::Int, index::Int=1)
    return get_channel(dataset.events[i_event], index)
end

# Set all exposures
function set_all_exposures(dataset::Dataset, t::Float64)
    dataset.exposure_override = true
    empty!(dataset.exposure_channel)
    for ch in 1:988
        dataset.exposure_channel[ch] = t
    end
end

function set_all_active_exposures(dataset::Dataset, t::Float64)
    dataset.exposure_override = true
    for ch in keys(dataset.exposure_channel)
        dataset.exposure_channel[ch] = t
    end
end

function set_all_active_exposures_from_file(dataset::Dataset)
    if !dataset.exposure_override
        @warn "You are calling a function to override the exposure with option to override disabled"
    end
    
    try
        # use CSV.jl for large throughput files
        data = CSV.read(dataset.file_exposure, 
                    DataFrame;
                    header=false,
                    delim='\t',
                    ignorerepeated=true,
                    types=[Int, Float64])
    
        # rename columns
        rename!(data, [:Index, :Exposure])
    
        for i in 1:nrow(data)
            chan = Int(data[i, :Index])
            exposure = Float64(data[i, :Exposure])
        
            if chan < 1 || chan > 988
                error("You are trying to set the exposure of NON bolometric channels: $chan")
            end  
            if !haskey(dataset.exposure_channel, chan) && exposure > 0.0
                error("Channel $chan previously was not in the exposure map, you are trying to add it")
            end
            dataset.exposure_channel[chan] = exposure
        end
        
        # Remove channels with zero exposure
        channels_to_remove = Int[]
        for (ch, exp) in dataset.exposure_channel
            if exp == 0.0
                push!(channels_to_remove, ch)
            end
        end
        for ch in channels_to_remove
            delete!(dataset.exposure_channel, ch)
        end

    catch e
        error("Could not read input ASCII exposure file $(dataset.file_exposure): $e")
    end
end


function set_livetime_fraction_from_file(dataset::Dataset)
    @warn "set_livetime_fraction_from_file function is not implemented yet."
    # try
    #     data = readdlm(dataset.livetime_filename, skipstart=1)  
    #     for row in eachrow(data)
    #         if length(row) >= 4
    #             chan = Int(row[1])
    #             frac = Float64(row[2])
    #             invalid = occursin("1.0", string(row[4]))  
                
    #             if haskey(dataset.exposure_channel, chan)
    #                 if invalid
    #                     dataset.exposure_channel[chan] = 0.0
    #                 else
    #                     dataset.exposure_channel[chan] *= frac
    #                 end
    #             end
    #         end
    #     end
        
    # catch e
    #     error("Could not read livetime fraction file $(dataset.livetime_filename): $e")
    # end
end

# 获取总曝光
function get_total_exposure(dataset::Dataset)
    exposure = 0.0
    for exp_val in values(dataset.exposure_channel)
        exposure += exp_val
    end 
    if dataset.exposure_override
        @warn "Exposures were written manually with a setter. The result might be inaccurate. DEBUG USE ONLY!!!"
    end
    return exposure
end

# used for toy-MC and sensitivity studies, to be implemented later
function compute_expected_bkg(
    dataset::Dataset,
    min_energy::Float64,
    max_energy::Float64;
    bkg_rate::Float64=0.0,
    bkg_rate_err::Float64=0.0,
    roi_width::Float64=15.0,
    blinded::Bool=true
)
    tot_range = max_energy - min_energy 
    if bkg_rate != 0.0
        bkg = bkg_rate * tot_range * get_total_exposure(dataset)
        if bkg > sqrt(bkg_rate * tot_range * get_total_exposure(dataset))
            bkg_err = bkg_rate_err * tot_range * get_total_exposure(dataset)
        else
            bkg_err = sqrt(bkg)
        end
        return bkg, bkg_err
    end
    
    bkg = 0.0
    bkg_err = 0.0
    
    # Retrieve the energy of 60Co and Qbb
    co60 = PhysicalConstants.ESumPeak60Cobalt
    qbb = PhysicalConstants.Qvalue0nbb130Te
    
    allowed_regions = Tuple{Float64, Float64}[]
    if roi_width > co60 - min_energy || roi_width > max_energy - qbb
        error("ROI width is way too large")
    end
    if roi_width < qbb - co60
        push!(allowed_regions, (min_energy, co60 - 0.5 * roi_width))
        push!(allowed_regions, (co60 + 0.5 * roi_width, qbb - 0.5 * roi_width))
        push!(allowed_regions, (qbb + 0.5 * roi_width, max_energy))
    else
        push!(allowed_regions, (min_energy, co60 - 0.5 * roi_width))
        push!(allowed_regions, (qbb + 0.5 * roi_width, max_energy))
    end
    
    # Loop over events
    for event in dataset.events
        if event.multiplicity > 1
            continue
        end
        
        energy_val = get_energy(event, 1)
        for (region_min, region_max) in allowed_regions
            if energy_val > region_min && energy_val < region_max
                bkg += 1.0
                break
            end
        end
    end
    
    bkg_err = sqrt(bkg)
    
    used_range = 0.0
    for (region_min, region_max) in allowed_regions
        used_range += (region_max - region_min)
    end
    
    bkg *= tot_range / used_range
    bkg_err *= tot_range / used_range
    
    return bkg, bkg_err
end


export Dataset, Event
export read_data, get_energy, get_channel, get_total_exposure, compute_expected_bkg
export set_all_efficiencies, set_all_exposures, set_all_active_exposures, set_all_active_exposures_from_file, set_livetime_fraction_from_file

end #module
