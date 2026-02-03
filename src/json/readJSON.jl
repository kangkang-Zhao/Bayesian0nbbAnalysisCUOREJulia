using JSON
using LinearAlgebra

j = JSON.parsefile("lineshape_scaling_ds3021.json")

pQ   = Vector{Float64}(j["Q"]["params"])
errQ = Vector{Float64}(j["Q"]["errors"])
cov_raw = j["Q"]["cov"]
covQ = reduce(vcat, (Float64.(row)' for row in cov_raw))

chi2 = j["Q"]["chi2"]
ndf  = j["Q"]["ndf"]

@show pQ
@show covQ
@show chi2, ndf

pSigma   = Vector{Float64}(j["Sigma"]["params"])
errSigma = Vector{Float64}(j["Sigma"]["errors"])
cov_rawSigma = j["Sigma"]["cov"]
covSigma = reduce(vcat, (Float64.(row)' for row in cov_rawSigma))

chi2Sigma = j["Sigma"]["chi2"]
ndfSigma  = j["Sigma"]["ndf"]

@show pSigma
@show covSigma
@show chi2Sigma, ndfSigma

