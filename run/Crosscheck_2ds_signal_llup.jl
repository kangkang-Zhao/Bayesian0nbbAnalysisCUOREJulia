include("../readHist.jl")
include("../PhysicalConstants.jl")
include("../Lineshape.jl")
include("../Dataset.jl")
include("../HistogramAsUvDistribution.jl")
include("../Util.jl")

using .LineshapeModule
using .DatasetModule
using .PhysicalConstants
using .readHistModule
# using .HistoPriorModule
using .UtilsModule
LS1 = LineshapeGauss3(
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape/calibration/PeakShape_ds3819.root",
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/baseline_sigmas/bkg/ds3819_baseline.root",
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/baseline_sigmas/calib/ds3819_baseline.root",
    3819)
println(LS1.fResoModel)
# set_poly_scaling(LS1, kQvalue, "/Users/zhaokangkang/gssiwork/julia-dev/code-test/test_json/lineshape_scaling_ds3021.json",2)
# set_poly_scaling(LS1, kSigma, "/Users/zhaokangkang/gssiwork/julia-dev/code-test/test_json/lineshape_scaling_ds3021.json",2)
SetHistoScaling(LS1, kQvalue)
SetHistoScaling(LS1, kSigma)
println(LS1.fResoModel)
DS1 = DatasetModule.Dataset(LS1,
    3819,
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/super_reduced/background/unblinded/SuperReduced_Background_ds3819.root",
    0.947,
    0.0075,
    0.88345,
    0.00085,
    237.0,
    1.0,
    0.0,
    1.0,
    0.0,
    2465.0,
    # 2520.0,
    2575.0,
    true,
    false,
    false,
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/CombinedEfficiency/CombinedEfficiencies.root",
    "",
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/Exposures/Exposures_ds3819.txt",
    "",
    true,
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape_scaling_output/CombinedReso_model2_pol2bias.root")

LS2 = LineshapeGauss3(
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape/calibration/PeakShape_ds3820.root",
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/baseline_sigmas/bkg/ds3820_baseline.root",
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/baseline_sigmas/calib/ds3820_baseline.root",
    3820)
println(LS2.fResoModel)
# set_poly_scaling(LS2, kQvalue, "/Users/zhaokangkang/gssiwork/julia-dev/code-test/test_json/lineshape_scaling_ds3021.json",2)
# set_poly_scaling(LS2, kSigma, "/Users/zhaokangkang/gssiwork/julia-dev/code-test/test_json/lineshape_scaling_ds3021.json",2)
SetHistoScaling(LS2, kQvalue)
SetHistoScaling(LS2, kSigma)
println(LS2.fResoModel)
DS2 = DatasetModule.Dataset(LS2,
    3820,
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/super_reduced/background/unblinded/SuperReduced_Background_ds3820.root",
    0.9155,
    0.0095,
    0.88345,
    0.00085,
    301.0,
    1.0,
    0.0,
    1.0,
    0.0,
    2465.0,
    2575.0,
    true,
    false,
    false,
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/CombinedEfficiency/CombinedEfficiencies.root",
    "",
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/Exposures/Exposures_ds3820.txt",
    "",
    true,
    "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape_scaling_output/CombinedReso_model2_pol2bias.root"
)

datasets = Vector{DatasetModule.Dataset}()
# using push! function for individual elements or append! for multiple elements
push!(datasets, DS1)
push!(datasets, DS2)


# using Distributions, IntervalSets
using OrderedCollections
datasets_group = [[DS1], [DS2]]
BI_priors = OrderedDict()
numm2::Int32 = 1
for datasets in datasets_group
    BI_priors[Symbol("BI_$(numm2)")] = BI_prior(datasets)
    global numm2 += 1
end

priors = OrderedDict()
priors[Symbol("Gamma_0ν")] = signal_prior(datasets)
merge!(priors, BI_priors)
priors[Symbol("Co60")] = Co60_prior(datasets)
merge!(priors,cutEff_prior(datasets))
priors[Symbol("Qββ")] = Qββ_prior()
prior_tmp = LineshapeScaling_prior(datasets, true)
for i in 1:1:length(prior_tmp)
    merge!(priors, prior_tmp[i])
end

using BAT
using Distributions
using Random
using StatsPlots
using ValueShapes
using LinearAlgebra
using Statistics
using DensityInterface

prior_my = distprod(; priors...)
# sensExposure = active_exposure_sum(datasets)[2]
f_norm = PhysicalConstants.N_A * 1000. * PhysicalConstants.Abundance_130Te / PhysicalConstants.mass_TeO2
# f_norm = PhysicalConstants.N_A * 1000. * PhysicalConstants.Abundance_130Te * sensExposure / PhysicalConstants.mass_TeO2
likelihood = DensityInterface.logfuncdensity(
    function (param::NamedTuple)
        total_ll = 0.0

        par_signal = param[:Gamma_0ν]
        par_co60 = param[Symbol("Co60")]
        par_Qbb = param[:Qββ]
        numm::Int32 = 1
        par_BI = []
        for datasets in datasets_group
            push!(par_BI, param[Symbol("BI_$(numm)")])
            numm += 1
            for ds in datasets
                par_bi = par_BI[numm-1]
                par_bias = param[Symbol("bias_ds$(ds.ds)")]*1.0
                par_reso = param[Symbol("reso_ds$(ds.ds)")]*1.0
                par_effcut = param[Symbol("EffCut_ds$(ds.ds)")]*1.0
                ls = ds.lineshape
                SetTmpPolyScaling(ls, kQvalue, [par_bias])
                SetTmpPolyScaling(ls, kSigma, [par_reso])
                for (ch, expo) in ds.exposure_channel
                    # eff = par_effcut * ds.containment_efficiency * ds.trigger_efficiency_channel[ch]
                    s = par_signal * f_norm * par_effcut * ds.containment_efficiency * ds.trigger_efficiency_channel[ch] * ds.exposure_channel[ch]
                    b = par_bi * (ds.emax - ds.emin) * par_effcut * ds.exposure_channel[ch]
                    # c = par_co60 * par_effcut * exp(-(ds.delta_t-datasets[1].delta_t) / PhysicalConstants.Tau60Cobalt) * ds.trigger_efficiency_channel[ch] * ds.exposure_channel[ch]
                    c = par_co60 * par_effcut * exp(-(ds.delta_t) / PhysicalConstants.Tau60Cobalt) * ds.trigger_efficiency_channel[ch] * ds.exposure_channel[ch]
                    λ = s + b + c
                    total_ll -= λ
                    if !haskey(ds.events_channel, ch)
                        continue
                    else
                        events = ds.events_channel[ch]
                    end

                     
                    Qco = ls.fParChannel[ch].Qvalue * PhysicalConstants.ESumPeak60Cobalt / PhysicalConstants.E208Tlpeak + par_bias
                    Qbb = ls.fParChannel[ch].Qvalue * par_Qbb / PhysicalConstants.E208Tlpeak + par_bias

                    bpar = ls.fParBaselineChannel[ch]
                    bkgVar = bpar.BkgVariance
                    calVar = bpar.CalVariance
                    reso =  sqrt(bkgVar + par_reso * (ls.fParChannel[ch].Sigma^2 - calVar))
                    invdenom = 1.0/(2.0*reso^2)
                    norm = (1.0 + ls.fParChannel[ch].Ratio + ls.fParChannel[ch].Ratio2) * sqrt(2*pi) * reso
                    for ev in events
                        energy = get_energy(ev)
                        # total_ll += log(
                        #     s * GetResponsePDF(ls, ch, energy, par_Qbb*1.0)
                        #     + c * GetResponsePDF(ls, ch, energy, PhysicalConstants.ESumPeak60Cobalt)
                        #     + b/(ds.emax-ds.emin)
                        # )
                                           
                        total_ll += log(
                                        s * (ls.fParChannel[ch].Ratio * exp(-((energy-ls.fParChannel[ch].EnergyRatio*Qbb)^2)*invdenom) + ls.fParChannel[ch].Ratio2 * exp(-((energy-ls.fParChannel[ch].EnergyRatio2*Qbb)^2)*invdenom) + exp(-((energy-Qbb)^2)*invdenom)) / norm
                                        + c * (ls.fParChannel[ch].Ratio * exp(-((energy-ls.fParChannel[ch].EnergyRatio*Qco)^2)*invdenom) + ls.fParChannel[ch].Ratio2 * exp(-((energy-ls.fParChannel[ch].EnergyRatio2*Qco)^2)*invdenom) + exp(-((energy-Qco)^2)*invdenom)) / norm
                                        + b / (ds.emax-ds.emin)
                                       )
                    end
                end
            end
        end
        return total_ll
    end,
)
posterior = PosteriorDensity(likelihood, prior_my)

using CPUTime
algorithm = MCMCSampling(
    mcalg = MetropolisHastings(),
    # mcalg = HamiltonianMC(),
    nsteps = 10^5,
    nchains = 5,
    burnin = MCMCMultiCycleBurnin(nsteps_per_cycle = 10000, max_ncycles = 50),  # 简化burnin
    # convergence = BrooksGelmanConvergence(threshold = 1.3)  #loose the convergence criteria
)
CPUtic()
@time @CPUtime samples = bat_sample(posterior, algorithm)
cpu_time = CPUtoc()
println("Run completed in $cpu_time seconds")

using Optim
#Comparison of Truth and Best Fit
samples_mode = mode(samples.result)
samples_mode isa NamedTuple

findmode_result = bat_findmode(
    posterior,
    OptimAlg(optalg = Optim.NelderMead(), init = ExplicitInit([samples_mode]))
)
# 分析结果
println("\nBAT拟合参数值:")
fit_par_values = findmode_result.result
println(fit_par_values)
# 可视化结果
using Plots

for key in keys(priors)
    println("Key: ", typeof(key))
    # 创建变量名
    # var_name = Symbol("p_$i")
    # 创建并赋值
    p = plot(size=(800,500), layout=(1,1), labelfontsize=12, tickfontsize=10, legendfontsize=7)
    pars = []
    par = key
    p = plot!(samples.result, par, subplot=1, label = "Posterior")
    idx = nothing
    for samp in samples.result
        v = samp.v
        weight = samp.weight
        for w = 1:1:weight
            if (idx == nothing)
                append!(pars, v[par])
            else
                append!(pars, v[par][idx])
            end
        end
    end
    x = range(minimum(pars), stop = maximum(pars), length = 100)
    y = pdf(prior_my[par], x)
    plot!(x, y, label = "prior", color = "black")
    savefig(p, "./result/$(par)_withSignal_Feb2_llup.pdf")

end

import HDF5
bat_write("./result/SignalIncluded_Feb2_llup.h5", samples.result)

