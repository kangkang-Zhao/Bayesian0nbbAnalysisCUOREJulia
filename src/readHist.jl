module readHistModule
export readHistModule

using DelimitedFiles
using Plots
using StatsBase: Histogram

export read_root_hist, get_percentile_bounds, truncate_histogram, normalize_histogram
"""
    read_root_hist(file::String, histname::String; exe="./read_hist")

Calls the C++ program `read_hist` to read a histogram from the root file.
return (bin_edges, contents)：
- bin_edges: Vector{Float64}, n+1 bin edges
- contents: Vector{Float64}, n bin content
"""
function read_root_hist(file::String, histname::String; exe="/Users/zhaokangkang/gssiwork/julia-dev/code-test/cpp/readHist")
    # Calls the C++ program and captures the output
    output = read(`$exe $file $histname`, String)

    # 按行分割，去掉注释行
    lines = filter(l -> !startswith(l, "#"), split(output, '\n'))

    lows, highs, vals = Float64[], Float64[], Float64[]
    for line in lines
        isempty(line) && continue
        parts = split(line, ",")
        push!(lows, parse(Float64, parts[1]))
        push!(highs, parse(Float64, parts[2]))
        push!(vals, parse(Float64, parts[3]))
    end

    # 构造 bin_edges (从 lows[1] 到 highs[end])
    edges = vcat(lows, highs[end])
    return edges, vals
    """
    # 构造 StatsBase.Histogram
    h = Histogram(edges, vals)
    return h
    """
end


"""
    get_percentile_bounds(h::Histogram, percentile::Float64; side::Int=0)

get the boundaries (minval, maxval) of a histogram with the given percetile

- percentile: range(0, 1), you want to preserve, 0.997 corresponds to 99.7%
- side = -1: remove from the right side (preserve the left side of percentile, 0~percetile)
- side =  0: remove from the both sides (kepp 1/2-percentile/2 ~ 1/2+percentile/2)
- side =  1: remove from the left side (preserve the right side of percentile, 1-percetile, 1)
"""
function get_percentile_bounds(h::Histogram, percentile::Float64; side::Int=0)
    """
    # 展开 histogram 为数据点 (根据 bin 内容重复)
    values = Float64[]
    println(length(h.weights))
    println(length(h.edges[1]))
    edges = h.edges[1]
    for i in 1:length(h.weights)
        n = h.weights[i]
        append!(values, ((edges[i]+edges[i+1])/2, n))
    end

    if isempty(values)
        return (0.0, 0.0)
    end
    """

    if percentile < 0.0 || percentile > 1.0
        @error "percentile must be in the range [0, 1]. Clamping to valid range."
    end
    # calculate the cumulative distribution function (CDF)    
    edges = h.edges[1]             # bin edges
    weights = h.weights            # bin contents
    cdf = cumsum(weights) ./ sum(weights)


    if side == -1
        qs = [0.0, percentile]
    elseif side == 0
        qs = [(1 - percentile)/2, (1 + percentile)/2]
    elseif side == 1
        qs = [1 - percentile, 1.0]
    else
        @warn "side must be -1, 0, or 1. Defaulting to 0."
        qs = [(1 - percentile)/2, (1 + percentile)/2]
    end

    #bounds = quantile(values, qs)
    #return (bounds[1], bounds[2])

    function quantile_from_cdf(q)
        idx = findfirst(>=(q), cdf)
        return edges[idx]   # return the left edge of the bin
    end

    return (quantile_from_cdf(qs[1]), quantile_from_cdf(qs[2]))
end


"""
    truncate_histogram(h::Histogram, percentile::Float64; side::Int=0)

truncate histogram, keep only the data in the range of given `percentile`
return the new truncated histogram。
"""
function truncate_histogram(h::Histogram, percentile::Float64; side::Int=0)
    if percentile < 0.0 || percentile > 1.0
        @error "percentile must be in the range [0, 1]. Clamping to valid range."
    end
    if percentile == 1.0
        return h
    end
    minval, maxval = get_percentile_bounds(h, percentile; side=side)

    # find the index of the boundaries in the original histogram
    lows = h.edges[1][1:end-1]
    highs = h.edges[1][2:end]
    mask = (highs .> minval) .& (lows .< maxval)

    # new boundaries
    new_edges = vcat(max(minval, lows[findfirst(mask)]), highs[mask])
    new_weights = h.weights[mask]

    return Histogram(new_edges, new_weights)
end

function normalize_histogram(h::Histogram)
    total = sum(h.weights)
    if total == 0
        return h
    end
    return Histogram(h.edges, h.weights ./ total)
end



function main()
    histname = "reso_Qbb_ds3819"
    filename = "/Users/zhaokangkang/gssiwork/julia-dev/data/taup23_2ds/lineshape_scaling_output/CombinedReso.root"
    println("Reading histogram '$histname' from $filename ...")
    edges, vals = read_root_hist(filename, histname)
    
    outfile = "$histname.png"
    println("Plotting and saving to $outfile ...")

    h = Histogram(edges, vals)
    # Draw
    plot(h, xlabel="x", ylabel="counts", title="Histogram from ROOT")
    savefig(outfile)
    println("Done.")
end

# Run the main function if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end