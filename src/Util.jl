module UtilsModule
export UtilsModule


using Distributions
using Statistics
using StatsBase
using IntervalSets
using OrderedCollections
using EmpiricalDistributions
using LinearAlgebra
using ValueShapes


using ..LineshapeModule
using ..DatasetModule
using ..PhysicalConstants
using ..HistoPriorModule
using ..readHistModule

# sum over all groups' ds-ch of Effi*Exposure
function active_exposure_sum(datasets::Vector{DatasetModule.Dataset})
    total_expo::Float64 = 0.0
    senst_expo::Float64 = 0.0
    for dataset in datasets
        num::Int64 = 0
        for (ch, expo) in dataset.exposure_channel
            eff::Float64 = get(dataset.total_efficiency_channel, ch, 0)    
            # expo::Float64 = get(dataset.exposure_channel,ch, 0) 
            if (eff==0.0 || expo == 0.0)
                num += 1
            end
            total_expo += expo
            senst_expo += eff * expo
        end
        println("Dataset $(dataset.ds): $num channels with zero efficiency or zero exposure in ROI.")
    end
    return total_expo, senst_expo
end

function total_exposure_sum(datasets::Vector{DatasetModule.Dataset})
    total_expo::Float64 = 0.0
    senst_expo::Float64 = 0.0
    for dataset in datasets
        num::Int64 = 0
        for ch in (1:1:988)
            eff::Float64 = get(dataset.total_efficiency_channel, ch, 0)    
            expo::Float64 = get(dataset.exposure_channel,ch, 0) 
            if (eff==0.0 || expo == 0.0)
                num += 1
            end
            total_expo += expo
            senst_expo += eff * expo
        end
        println("Dataset $(dataset.ds): $num channels with zero efficiency or zero exposure in total.")

    end
    return total_expo, senst_expo
end

#------------------------------------------------
function signal_prior(datasets::Vector{DatasetModule.Dataset})
    sensExposure = active_exposure_sum(datasets)[2]
    # compute f_norm
    # Abundance of 130Te is considered constant here
    # i.e., included in fNorm directly
    # if it is to be floated, then it has to be included
    f_norm = PhysicalConstants.N_A * 1000. * PhysicalConstants.Abundance_130Te * sensExposure / PhysicalConstants.mass_TeO2

    nSgnCandidates::Float64 = 0.
    nBkgCandidates::Float64 = 0.
    for dataset in datasets
        for ev in dataset.events
            eventEnergy = ev.energy[1]
            if eventEnergy > PhysicalConstants.Emin0nbbPeak20keV &&
               eventEnergy < PhysicalConstants.Emax0nbbPeak20keV
                nSgnCandidates += 1.
            else
                nBkgCandidates += 1.
            end
        end
    end
    # width in keV of signal region
    sgnRange = PhysicalConstants.Emax0nbbPeak20keV - PhysicalConstants.Emin0nbbPeak20keV
    # Emax and Emin can be put outside if all datasets have the same
    fitRange = datasets[1].emax - datasets[1].emin  # assuming all datasets have the same fit range
    
    errNSgnCandidates = sqrt( nSgnCandidates + nBkgCandidates * sgnRange / (fitRange-sgnRange) )
    nSgnCandidates -= nBkgCandidates * sgnRange / (fitRange-sgnRange)
    
    minR = nSgnCandidates - 15.0 * errNSgnCandidates
    if minR < 0.0
        minR = 0.0
    end
    maxR = nSgnCandidates + 15.0 * errNSgnCandidates
    if maxR < 3.0
        maxR = 15.0
    end
    minR /= f_norm
    maxR /= f_norm

    distrS = minR .. maxR

    return distrS
end

#------------------------------------------------
# Add BI cts/kev/kg/yr parameters, shared among all datasets
#****************  B I  ***************************************
function BI_prior(datasets::Vector{DatasetModule.Dataset})
    nEvents::Float64 = 0.
    for dataset in datasets
        for ev in dataset.events
            eventEnergy = ev.energy[1]
            if eventEnergy < PhysicalConstants.Emin0nbbPeak40keV ||
                eventEnergy > PhysicalConstants.Emax0nbbPeak40keV
                nEvents += 1.
            end
        end
    end
    # width in keV of signal region
    exclusionRange = PhysicalConstants.Emax0nbbPeak40keV - PhysicalConstants.Emin0nbbPeak40keV
    fitRange = datasets[1].emax - datasets[1].emin  # assuming all datasets have the same fit range
    # # of bkg events = N_called outside 4 sigma region around Qbb times Fraction(= fitRange / (fitRange - exclusionRange))
    nEvents *= fitRange / (fitRange - exclusionRange)
    # HARD-CODED AVERAGE VALUE FOR TOTAL DETECTOR EFFICIENCY!!!
    nEvents *= 0.95; 
    # err_bkg_events = sqrt( N_called ) * fitRange / (fitRange - exclusionRange )
    errNEvents = sqrt( nEvents * fitRange / (fitRange - exclusionRange) )
    
    minNEvents = nEvents - 15. * errNEvents
    if minNEvents < 0.0 
        minNEvents = 0.0
    end
    maxNEvents = nEvents + 15. * errNEvents
    totExposure = active_exposure_sum(datasets)[1]  
    minBI = minNEvents / (fitRange*totExposure)
    maxBI = maxNEvents / (fitRange*totExposure)
    distrBI = minBI .. maxBI
    return distrBI
end


# end of BI prior
#*******************************************************  E N D   O F   B I

#------------------------------------------------
# Add Linear Bkg parameters, shared among all datasets
#********* L I N E A R   B K G ********************************************** 
function BISlope_prior(datasets::Vector{DatasetModule.Dataset})
    Emax::Float64 = 0.0
    Emin::Float64 = 0.0
    Index::Int32 = 0
    for dataset in datasets
        if Index == 0
            Emax = dataset.emax
            Emin = dataset.emin
            Index += 1
            continue
        end
        if dataset.emax > Emax
            Emax = dataset.emax
        end
        if dataset.emin > Emin
            Emin = dataset.emin
        end
        Index += 1
    end
    # fMidRange = 0.5 * ( datasets[1].Emax + datasets[1].Emin )
    dE = 0.5 * (Emax - Emin)
    minSlope = -1.0 / dE
    maxSlope = 1.0 / dE
    distrBISlope = minSlope .. maxSlope
    return distrBISlope
end

# end of Linear Bkg prior
#********* E N D   O F   L I N E A R   B K G ********************************************** 

#------------------------------------------------
function Co60_prior(datasets::Vector{DatasetModule.Dataset})
#set one 60Co parameter shared
    # Floating 60Co mean value FIXME--hardcoded

    A::Float64 = 0.0 # [60Co events / kg / yr]
    errA::Float64 = 0.0
    B::Float64 = 0.0 # [bkg  events / kg / yr]
    errB::Float64 = 0.0
    totExp::Float64 = 0.0

    minCoW::Float64   = 2497.6
    maxCoW::Float64   = 2517.6
    max0nbbW::Float64 = 2537.6

    # Loop over datasets
    for dataset in datasets
        cutEff::Float64 = dataset.cut_efficiency
        CoReduction::Float64 = exp(-dataset.delta_t / PhysicalConstants.Tau60Cobalt)
        N60Co::Float64  = 0.0
        Nbkg::Float64   = 0.0
        energy::Float64 = 0.0
        # Loop over all channels
        for (ch, events) in dataset.events_channel
            for ev in events
                energy = ev.energy[1]
                if energy > minCoW && energy < maxCoW
                    N60Co += 1.0
                end
                if energy < minCoW || energy > max0nbbW
                    Nbkg += 1.0
                end
            end
        end
        A    +=      N60Co  / cutEff / CoReduction
        errA += sqrt(N60Co) / cutEff / CoReduction
        B    +=      Nbkg   / cutEff / CoReduction
        errB += sqrt(Nbkg)  / cutEff / CoReduction
    end# End of loop over datasets
    totExp = active_exposure_sum(datasets)[1]

    A /= totExp
    errA /= totExp
    B /= totExp
    errB /= totExp

    dE = maxCoW - minCoW
    fitRange = datasets[1].emax - datasets[1].emin
    DE = fitRange - ( max0nbbW - minCoW )
    R = A - B * dE / DE
    errR = sqrt(errA^2 + (errB * dE / DE)^2)
    minR = R - 10. * errR
    maxR = R + 10. * errR
    if( minR < 0. )
        minR = 0.
    end
        
    distrCo60 = minR .. maxR
    return distrCo60

end
#------------------------------------------------


function Co60Mean_prior()
    distrCo60Mean = 2497. .. 2517.
    return distrCo60Mean

end

# --------------------------------------------------------------
# Add a nuisance parameter for the cut efficiency, if required

function cutEff_prior(datasets::Vector{DatasetModule.Dataset})
    Eff_priors = OrderedDict()
    # Check if efficiency uncertainty is provided
    flag = true
    # Loop over datasets
    for ds in datasets
        if ds.cut_efficiency_err == 0.0
            flag = false
        end
    end
    if !flag
        for ds in datasets
            if !ds.th1_efficiency
                error("Efficiency set as nuisance parameter but its uncertainty is not provided. Abort.")
            end
        end
    end

    # Loop over datasets
    for ds in datasets

        # --------------------------------------------------
        # Case 1: Gaussian efficiency prior
        if !ds.th1_efficiency

            eff    = ds.cut_efficiency
            effErr = ds.cut_efficiency_err

            minEff = eff - 5 * effErr
            maxEff = eff + 5 * effErr

            minEff = max(minEff, 0.0)
            maxEff = min(maxEff, 1.0)

            name      = "EffCut_ds$(ds.ds)"
            Eff_priors[Symbol(name)] = Truncated(Normal(eff, effErr), minEff, maxEff)
            
        # --------------------------------------------------
        # Case 2: Histogram-based efficiency prior
        else
            histo = truncate_histogram(ds.eff_prior,0.999)
            # d_histo = UvBinnedDist(histo)
            d_histo = HistoPriorModule.HistogramAsUvDistribution(histo)
            name    = "EffCut_ds$(ds.ds)"
            Eff_priors[Symbol(name)] = d_histo
        end
    end
    return Eff_priors
    
end
# --------------------------------------------------------------
# Add a nuisance parameter for the Monte Carlo efficiency
# (common to all datasets)
function MCEff_prior(datasets::Vector{DatasetModule.Dataset})
    # Check if MC efficiency and uncertainty are provided and consistent
    flag = true
    for ds in datasets
        if ds.containment_efficiency != PhysicalConstants.MCEfficiency ||
           ds.containment_efficiency_err  != PhysicalConstants.MCEfficiencyErr ||
           ds.containment_efficiency_err  == 0.0
            flag = false
        end
    end
    if !flag
        error(
            "Monte Carlo efficiency set as nuisance parameter but in the current formulation (PRL19 analysis):\n" *
            "This is a parameter common to all datasets (both value and error).\n" *
            "This error can mean that:\n" *
            "Option 1: the MC efficiency value or the error differs from PhysicalConstants0nbb\n" *
            "Option 2: the uncertainty on the MC efficiency is not provided (null). Abort."
        )
    end
    minMCEff = PhysicalConstants.MCEfficiency - 5.0*PhysicalConstants.MCEfficiencyErr
    maxMCEff = PhysicalConstants.MCEfficiency + 5.0*PhysicalConstants.MCEfficiencyErr

    println(
        "This is the range for the MC efficiency (added as NP to the Model for the ML fit): ",
        minMCEff,
        " ",
        maxMCEff
    )
    return Truncated(Normal(PhysicalConstants.MCEfficiency, PhysicalConstants.MCEfficiencyErr), minMCEff, maxMCEff)
end
# ------------------------------------
# Add a parameter for Qbb
function Qββ_prior()
    minQββ = PhysicalConstants.Qββ130Te - 5.0*PhysicalConstants.Qββ130TeErr
    maxQββ = PhysicalConstants.Qββ130Te + 5.0*PhysicalConstants.Qββ130TeErr
    return Truncated(Normal(PhysicalConstants.Qββ130Te, PhysicalConstants.Qββ130TeErr), minQββ, maxQββ)
end
# ------------------------------------------------------
# Add a parameter for the isotopic fraction
function IsoFrac_prior()

    minAb = PhysicalConstants.NaturalIsotopicAbuTe130 - 5.0*PhysicalConstants.NaturalIsotopicAbuTe130Err
    maxAb = PhysicalConstants.NaturalIsotopicAbuTe130 + 5.0*PhysicalConstants.NaturalIsotopicAbuTe130Err
    
    return Truncated(Normal(PhysicalConstants.NaturalIsotopicAbuTe130, PhysicalConstants.NaturalIsotopicAbuTe130Err), minAb, maxAb)
   
end
# ------------------------------------------------------
# Add a parameter for the PSA efficiency systematic
function SystEff_prior()
    @warn "You need to take care the prior building of systematical efficiency!"
    SystEffMin = 0.9
    SystEffMax = 1.2

    return SystEffMin .. SystEffMax
end


function trim_histo_zerobins(hist::Histogram)
    bin_start::Int64 = 1
    bin_end::Int64 = length(hist.weights)
    println("Original bin range: start = '$bin_start', end = '$bin_end'.")
    for i in 1:1:length(hist.weights)
        if hist.weights[i]==0.0
            continue
        end
        bin_start = i
        break
    end
    for i in length(hist.weights):-1:1
        if hist.weights[i]==0.0
            continue
        end
        bin_end = i
        break
    end
    println("Trimmed bin range: start = '$bin_start', end = '$bin_end'.")
    trimmed_hist = Histogram(hist.edges[1][bin_start:bin_end+1],hist.weights[bin_start:bin_end])
    return trimmed_hist
end

# --------------------------------------------------------------
# Lineshape scaling nuisance parameters
function LineshapeScaling_prior(datasets::Vector{DatasetModule.Dataset}, fHistoLineshape::Bool)

    @warn "Lineshape scaling still under development. Run at your own risk."
    # Loop over datasets
    bias_priors = OrderedDict()
    reso_priors = OrderedDict()
    Qval_priors = OrderedDict()
    Sigma_priors = OrderedDict()
    for ds in datasets

        # ======================================================
        # Case 1: histogram-based lineshape priors
        if fHistoLineshape

            # ---------------- Bias ----------------
            hist_bias = ds.bias_ls_scaling
            hist_bias = trim_histo_zerobins(hist_bias)
            hist_bias = truncate_histogram(hist_bias,0.999)
            bias_name = "bias_ds$(ds.ds)"
            # bias_priors[Symbol(bias_name)] = UvBinnedDist(hist_bias)
            bias_priors[Symbol(bias_name)] = HistogramAsUvDistribution(hist_bias)

            # ---------------- Resolution ----------------
            hist_reso = ds.reso_ls_scaling
            hist_reso = trim_histo_zerobins(hist_reso)
            hist_reso = truncate_histogram(hist_reso,0.999)
            reso_name = "reso_ds$(ds.ds)"
            # reso_priors[Symbol(reso_name)] = UvBinnedDist(hist_reso)
            reso_priors[Symbol(reso_name)] = HistogramAsUvDistribution(hist_reso)

        # ======================================================
        # Case 2: analytic lineshape with covariance matrices
        else
            ls = ds.lineshape
            # Q-value scaling
            pQval = ls.fPolyScalingMap[kQvalue]
            mQval = ls.fCovarianceMatrixPar[kQvalue]
            # Sigma scaling
            pSigma = ls.fPolyScalingMap[kSigma]
            mSigma = ls.fCovarianceMatrixPar[kSigma]

            mvQval = MvNormal(pQval, Symmetric(mQval))
            mvSigma = MvNormal(pSigma, Symmetric(mSigma))

            Qval_name = "Qval_ds$(ds.ds)"
            Qval_priors[Symbol(Qval_name)] = mvQval

            Sigma_name = "Sigma_ds$(ds.ds)"
            Sigma_priors[Symbol(Sigma_name)] = mvSigma
            
        end
    end
    return bias_priors, reso_priors, Qval_priors, Sigma_priors
end

export active_exposure_sum, total_exposure_sum,
       signal_prior, BI_prior, BISlope_prior,
       Co60_prior, Co60Mean_prior,
       cutEff_prior, MCEff_prior,
       Qββ_prior, IsoFrac_prior,
       SystEff_prior,
       LineshapeScaling_prior
end
# --------------------------------------------------------------
