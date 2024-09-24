#----------------------------------------------------------------------------------------------------------
#生成团簇网络
#----------------------------------------------------------------------------------------------------------
module Cluster
using Graphs, Distributions, JLD2, SimpleWeightedGraphs, DataFrames, CSV, ProgressBars, ProgressMeter

#网络信息的structure
struct netWork 
    x::Array{Float64}
    y::Array{Float64}
    A::Matrix{Float64}
    N::Int64
    Ncs::Array{Int64}
    n::Array{Int64}
    lmargin::Array{Float64}
    sizeV::Array{Int64}
end

#----------------------------------------------------------------------------------------------------------
#按概率函数f随机产生一个团簇大小    
function random_cluster_size(N, f_size; sizeLimit=50)
    sizeV = Int64[]
    Nsum::Int64 = 0
    while Nsum != N
        r1 = rand() * sizeLimit
        r2 = rand() * f_size(0)
        if r2 <= f_size(r1)
            size0 = convert(Int64, ceil(r1))
            push!(sizeV, size0)
            Nsum += size0
            if Nsum > N
                Nsum -= size0
                pop!(sizeV)
            end
        end
    end
    sizeV
end

#----------------------------------------------------------------------------------------------------------
#随机产生一个团簇
function cluster_generate(size, distance, sigma, region)
    x0, y0 = 0.0, 0.0 #设定初始点为原点
    x, y = Float64[x0], Float64[y0] #坐标数组
    Nlattice = 20 #一个RyR分格数
    cell = 30.0 / (2.0 * Nlattice + 1.0) #格的大小 
    NTrans = convert(Int64, ceil(region / cell)) # 格点数

    occupy1 = zeros(Int64, NTrans * 2, NTrans * 2)
    for j in NTrans-Nlattice:NTrans+Nlattice, k in NTrans-Nlattice:NTrans+Nlattice
        occupy1[j, k] = 1
    end
    i = 1
    while i < size
        phi = rand() * 2 * pi #随机方向
        r = rand(Normal(distance, sigma)) #距离，正规分布
        x0 += r * cos(phi) #平移x
        y0 += r * sin(phi) #平移y
        Mx, My = convert.(Int64, cld.(x0, cell)) + NTrans, convert.(Int64, cld.(y0, cell)) + NTrans #所在格子
        Mxl, Mxh, Myl, Myh = Mx - Nlattice, Mx + Nlattice, My - Nlattice, My + Nlattice #外扩xx个格子
        #println(Mxl,Mxh,Myl,Myh)

        if Mxl < 1 || Myl < 1 || Mxh > NTrans * 2 || Myh > NTrans * 2 #超出区域
            continue
        elseif sum(occupy1[Mxl:Mxh, Myl:Myh]) >= 1 #2Nlattice+1格中任意一点已经被占据
            continue
        else
            for j in Mxl:Mxh
                for k in Myl:Myh
                    occupy1[j, k] = 1
                end
            end
            push!(x, x0) #记录x
            push!(y, y0) #记录y
            i += 1
        end
    end
    return x, y
end

#----------------------------------------------------------------------------------------------------------
#判断一个团簇是否与已有的点重合且在范围内
function seperate(size, NC, Nx, Ny, occupy, regionl)
    N = min(max(convert(Int64, round(size * sqrt(NC) / 100)), 5),10)
    #N=2
    for i in 1:size #循环所有的点
        if Nx[i] > regionl || Ny[i] > regionl || Nx[i] < 1 || Ny[i] < 0  #在(0，0),(region,region)区域内
            return false
        end
        Nxl, Nxh, Nyl, Nyh = Nx[i] - N, Nx[i] + N, Ny[i] - N, Ny[i] + N #考虑该点所在及外扩N个格
        if Nxl < 1
            Nxl = 1
        end #防止i<=0
        if Nyl < 1
            Nyl = 1
        end
        if Nxh > regionl
            Nxh = regionl
        end #防止i<=0
        if Nyh > regionl
            Nyh = regionl
        end
        if sum(occupy[Nxl:Nxh, Nyl:Nyh]) >= 1 #(2N+1)^2格中任意一点已经被占据  
            return false
        end
    end
    true
end
#随机分布团簇    
function random_cluster_distribution(; sizeV, clusterDistance, distance, sigma)
    #sizeV团簇大小数组;clusterDistance团簇间距离; distance单个团簇参数; sigma单个团簇参数;lattice格点大小
    NC = length(sizeV) #团簇个数
    region = clusterDistance * sqrt(NC) #分布范围
    lattice=10
    regionl = convert(Int64, ceil(region / lattice)) # 格点数
    x, y, Ncs, n = Float64[], Float64[], Float64[], Int64[] #保存所有受体位置的数组
    occupy = zeros(Int64, regionl, regionl) #记录占据状况的数组
    i = 1
    label = 0
    x_cluster, y_cluster = Float64[], Float64[]
    println("Generate and distribut clusters")
    p = Progress(NC, dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:black)
    for i in 1:NC
        x_cluster, y_cluster = cluster_generate(sizeV[i], distance, sigma, region) #产生团簇
        label = 1
        while label >= 1
            if label > 100
                x_cluster, y_cluster = cluster_generate(sizeV[i], distance, sigma, region) #产生团簇
                #println(label)
            end
            label += 1
            xc = rand() * region #随机产生团簇中心x
            yc = rand() * region #随机产生团簇中心y
            x_c, y_c = x_cluster .+ xc, y_cluster .+ yc #平移团簇至随机产生的中心位置
            Nx, Ny = convert.(Int64, cld.(x_c, lattice)), convert.(Int64, cld.(y_c, lattice)) #判断团簇内各点所在格点

            if seperate(sizeV[i], NC, Nx, Ny, occupy, regionl) #判断一个团簇是否与已有的点重合且在范围内
                for ioccupy in 1:sizeV[i]
                    for xoc in -3:3, yoc in -3:3
                        if regionl > Nx[ioccupy] + xoc > 0 && regionl > Ny[ioccupy] + yoc > 0
                            occupy[Nx[ioccupy]+xoc, Ny[ioccupy]+yoc] = 1
                        end
                    end
                end #无重合则将团簇所在格点标记为已被占据
                append!(x, x_cluster .+ xc) #记录x
                append!(y, y_cluster .+ yc) #记录y
                append!(Ncs, fill(sizeV[i], sizeV[i], 1))
                append!(n, fill(i, sizeV[i], 1))
                label = 0
            end
        end
        next!(p)
    end



    lmargin = Float64[]
    N = 8
    for i in eachindex(x)
        Nx, Ny = convert.(Int64, cld.(x[i], lattice)), convert.(Int64, cld.(y[i], lattice)) #判断团簇内各点所在格点
        Nmargin = 0
        for xmargin in -N:N, ymargin in -N:N
            if regionl > xmargin + Nx > 0 && regionl > Ny + ymargin > 0
                Nmargin += occupy[Nx+xmargin, Ny+ymargin]
            end
        end
        push!(lmargin, (1.0 - Nmargin / (2 * N + 1)^2) * 2.)

    end
    return x, y, Ncs, n, lmargin

end

#实验数据
function data(lattice)
    data = CSV.read("data.csv", DataFrame)

    Cluster = []
    sizeV = Int64[]
    i0 = 1
    N = size(data)[1]
    for i in 2:N
        if typeof(data.Nc[i]) == Int64
            push!(Cluster, [data.x[i0:i-1], data.y[i0:i-1]])
            push!(sizeV, i - i0)
            i0 = i
        elseif i == N
            push!(Cluster, [data.x[i0:i], data.y[i0:i]])
            push!(sizeV, N - i0 + 1)
        end
    end
    permvec = sortperm(sizeV, rev=true)
    sizeV = sizeV[permvec]
    Cluster = Cluster[permvec]
    #println(sizeV,sum(sizeV))

    x, y = Float64[], Float64[] #坐标数组 

    for i in eachindex(sizeV)
        append!(x, Cluster[i][1])
        append!(y, Cluster[i][2])
    end


    Ncs, n = Float64[], Int64[]
    for i in eachindex(sizeV)
        append!(Ncs, fill(sizeV[i], sizeV[i], 1))
        append!(n, fill(i, sizeV[i], 1))
    end

    region = convert(Int64, round(max(maximum(x), maximum(y))))
    occupy = zeros(Int64, region, region) #记录占据状况的数组
    lmargin = Int64[]
    for i in eachindex(x)
        Nx, Ny = convert.(Int64, cld.(x[i], lattice)), convert.(Int64, cld.(y[i], lattice)) #判断团簇内各点所在格点
        Nmargin = 0
        for xmargin in -2:2, ymargin in -2:2
            if region > xmargin + Nx > 0 && region > Ny + ymargin > 0
                Nmargin += occupy[Nx+xmargin, Ny+ymargin]
            end
        end
        if Nmargin > 20
            push!(lmargin, 0)
        else
            push!(lmargin, 1)
        end

    end

    #print(lmargin)

    return x, y, Ncs, n, lmargin, sizeV, N
end

#----------------------------------------------------------------------------------------------------------
#转化几何分布为几何网络
function euclidean2net(x::Array{Float64}, y::Array{Float64})
    N = length(x) #节点数
    g = SimpleGraph(N, 0) #产生一个N个节点的空图

    for i in 1:N, j in i:N
        if i != j
            r = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) #计算两点距离
            ran = rand() #产生随机数
            if ran <= exp(-r / 30) * 100 # (0.777 * exp(-r / 4.3) + 0.217 * exp(-r / 15.5)  +1. * 0.0019 * exp(-r / 300.0)) * 200 #按几率判断，r0Ca为有链接的标度
                add_edge!(g, i, j) #添加链接
            end
        end
    end
    g
end
#转化为加权网络（weighted network）
function euclidean2net_w(x::Array{Float64}, y::Array{Float64})
    N = length(x) #节点数
    iv = Int64[]
    jv = Int64[]
    wv = Float64[]
    print("Transfer to a geometric network...")
    for i in 1:N, j in i:N
        if i != j
            r = sqrt((x[i] - x[j])^2 + (y[i] - y[j])^2) #计算两点距离
            w::Float64 = exp(-r / 60.0) / 1e-5   #tauRR变大，tau分布变大
            push!(iv, i)
            push!(jv, j)
            push!(wv, w)
        else
        end
    end
    net = SimpleWeightedGraph(iv, jv, wv)
    println("done")
    print("Save network...")
    return net
end

#----------------------------------------------------------------------------------------------------------
#产生网络，团簇网络和去掉rogue RyR的网络。weighted选择加权还是非加权网络
function random_Cluster_graph(N, paras; clusterDistance=250, distance=40, weighted=false, single=false, data=false)
    f_size(x) = 0.991 * exp(-0.66 * x) + 0.009 * exp(-0.017 * x)   #团簇大小概率函数
    #f_size(x) = 1. * exp(-0.966 * x) + 0 * exp(-0.017 * x)   #团簇大小概率函数
    sizeLimit = 100
    if single == true
        sizeLimit = 200
    end
    sizeV = sort(random_cluster_size(N, f_size, sizeLimit=sizeLimit), rev=true) #产生Nc个团簇的大小并排序
    if single == true
        sizeV = [N]
    end
    #plot_cluster_size(size,f_size) #画出随机产生团簇大小的结果
    if data == false
        @show sizeV, length(sizeV), sum(sizeV) #团簇个数,总节点数
        x, y, Ncs, n, lmargin = random_cluster_distribution(sizeV=sizeV, clusterDistance=clusterDistance, #0/sqrt(length(sizeV)),
            distance=distance, sigma=7.4) #产生分布
    else
        data = true
        x, y, Ncs, n, lmargin, sizeV, N = data(30)
        @show sizeV, length(sizeV), sum(sizeV) #团簇个数,总节点数
    end

    #net=load("results/Single0_net.jld","net") #读取cluster网络的数据
    #x,y,A,N,Ncs,n,lmargin,sizeV=net.x,net.y,net.A,net.N,net.Ncs,net.n,net.lmargin,net.sizeV

    if weighted == false
        net0 = euclidean2net(x, y)
    elseif weighted == true #转化为加权网络（weighted network）
        net0 = euclidean2net_w(x, y)
    end


    A::Matrix{Float64} = adjacency_matrix(net0)
    net=netWork(x,y,A,N,Ncs,n,lmargin,sizeV) 
    JLD2.save("results/" * paras * "_net.jld", "net", net)
    net
end


end
