using BAT, HDF5, Plots, Distributions

# 读取两个文件
samples1 = bat_read("./result/SignalIncluded_Feb2_llup.h5").result
samples2 = bat_read("./result/SignalIncluded_Jan30_llup_HMC.h5").result  # 请替换为您的第二个文件名

# 获取所有参数的键（假设两个文件有相同的参数结构）
keys_list = collect(keys(samples1.v[2]))
num_params = length(keys_list)

# 计算布局：每行最多3个图
ncols = min(3, num_params)
nrows = ceil(Int, num_params / ncols)

# 创建包含所有子图的大画布
p = plot(size=(3*400*ncols, 3*300*nrows), 
         layout=(nrows, ncols), 
         labelfontsize=12, 
         tickfontsize=10, 
         legendfontsize=9)

# 为每个参数绘制对比图
for (idx, key) in enumerate(keys_list)
    # 获取当前子图的行列位置
    row = ceil(Int, idx / ncols)
    col = idx - (row-1)*ncols
    
    # 在两个样本中绘制同一参数的分布
    plot!(p, samples1, key, 
          subplot=idx, 
          label="MH", 
          linewidth=2,
          title=string(key),
          titlefontsize=11)
    
    plot!(p, samples2, key, 
          subplot=idx, 
          label="HMC", 
          linewidth=2,
          linestyle=:dash)
            
    '''
    if haskey(prior_my, par)
        x_range = range(minimum(vcat(pars, [bkg_only.v[i][par] for i in 1:min(100, length(bkg_only.v))])),
                           maximum(vcat(pars, [bkg_only.v[i][par] for i in 1:min(100, length(bkg_only.v))])),
                           length=100)
        y = pdf(prior_my[par], x_range)
        plot!(p, x_range, y,
                  subplot=idx,
                  label="Prior",
                  linewidth=3,
                  color=:black,
                  linestyle=:dash)
        
    end
    '''
end

# 保存为单个PDF文件
savefig(p, "MH_HMC_comparison.pdf")

# 也可以保存为PNG以便预览
savefig(p, "all_parameters_comparison.png")


