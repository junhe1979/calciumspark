#----------------------------------------------------------------------------------------------------------
# 画图,可独立运行
#----------------------------------------------------------------------------------------------------------
using Plots,LaTeXStrings,JLD2,Colors,Compose,DelimitedFiles,Formatting
#网络信息的structure

module Cluster
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
end

#----------------------------------------------------------------------------------------------------------
function plot_res()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    for ii in 1:1
    paras=load("results/name.jld","name") #读取需要画图的数据文件名字
    #paras="CaN1000C"*string(ii)
    net=load("results/"*paras*"_net.jld","net") #读取cluster网络的数据
    N,size,x,y,A=net.N,net.sizeV,net.x,net.y,net.A
    x=x./1000.;y=y./1000. #将坐标缩小1000倍
    NC=length(size) #number of cluster
    dpi=1200 #画图精度
    markersize=1.;linewidth=0.00001; guifont=(8, :dark); tickfontsize=6 #点大小，线宽度，框架字体，刻度字体

    #画cluster网络图 
    #设置画图格式，空图
    pnet=plot(guidefont = guifont,framestyle=:box, ylabel=L"y~(\mu m)",dpi=dpi,
        xlims=(0.,maximum(x)),ylims=(0.,maximum(y)),tickfontsize=tickfontsize,legend = :none)
    #画连线 
    Amax=maximum(A) 
    lcolor=(Amax .- A)/Amax
    Aij=[A[i,j]  for i in 1:N for j in i:N]
    iv=[i  for i in 1:N for j in i:N]
    jv=[j  for i in 1:N for j in i:N]
    sort=sortperm(Aij)
    
    for s in sort 
        i=iv[s];j=jv[s]
        if A[i,j]>1000   
             plot!(pnet,[x[i],x[j]],[y[i],y[j]],legend = :none,linecolor=RGB(lcolor[i,j], lcolor[i,j], lcolor[i,j]),linewidth=A[i,j]*linewidth) #链接
        end
    end

    #println(lmargin)
    xx,yy=[],[]
    for i in eachindex(net.lmargin)
        if net.lmargin[i]==0
            push!(xx,x[i])
            push!(yy,y[i])
        end
    end

    #画点，不同团簇颜色不同

    NN=1
    for i in 1:NC
        scatter!(pnet,x[NN:NN+size[i]-1],y[NN:NN+size[i]-1], markercolor=i, 
            markersize=markersize .*(net.lmargin[NN:NN+size[i]-1].+1e-3),markerstrokewidth=0.,markershape=:circle) #节点
         NN=NN+size[i]
    end


    #paras=paras*string(ii)
    @show paras
    paras_plot = "results/" * paras * "1_plot.jld" #数据文件名字
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[3]
    Nnode = length(Ss[1])[1]
    ly = maximum(y) * 1.1
    lx = maximum(x) * 0.05
    p=[]    
    for i0 in 0:7
        i=convert(Int64, ceil(length(t)*i0/8*0.9))+10
        #println(i)
        S = Ss[i]
        xo, yo = [], []
        for inode = 1:Nnode
            if S[inode] == 1
                push!(xo, x[inode])
                push!(yo, y[inode])
                #println(xo,yo)
            end
        end
        p0=scatter(pnet, xo, yo, annotations=(lx, ly, Plots.text("t=" * string(format(t[i], precision=4)) * "s", 12, :left)),
            color=:red, markersize=markersize+1, markerstrokecolor=:yellow, markerstrokewidth=0.3, markershape=:star7)
        push!(p,p0)
    end


    #读取数据
    t,Ss=load(paras_plot,"t")[1],load(paras_plot,"Ss")[1] #时间，Ss中第一个数值
    #@show length(t), length(Ss[1])[1]
    Nnode=length(Ss[1])[1] 
    cg=cgrad([:white, :green,:black], [0, 0.5,0.75, 1.])
 
    #平均值， 
    x0,y0=t, [sum(S)/length(S)[1] for S in Ss] #时间，平均值
    xlim,ylim=maximum(x0),maximum(y0)*1.1 #取值范围
    p10=plot(x0,y0, markersize=markersize, markerstrokewidth=0, dpi=dpi,legend = :none,
        xlims=(0,xlim),ylims=(0,ylim), guidefont = guifont,tickfontsize=tickfontsize,
        framestyle=:box, ylabel=L"[Ca^{2+}]~(\mu M)") #画图
    #节点演化图，
    x0,y0,z0=t,[i/1000. for i in 1:Nnode],[S[i] for i in 1:Nnode, S in Ss] #时间，节点序号，浓度
    cmin,cmax=minimum(z0),maximum(z0)+1e-20 #浓度范围
    p1=heatmap(x0, y0, z0, c = cg,xlims=(0.,maximum(x0)),ylims=(0.,maximum(y0)),clims=(cmin,cmax), 
        framestyle=:box, guidefont = guifont,tickfontsize=tickfontsize,
        dpi=dpi,ylabel=L"n_{nodes}~(~\times 1000~)") #画图

    #读取数据
    t,Ss=load(paras_plot,"t")[2],load(paras_plot,"Ss")[2] 

    #平均值，  
    x0,y0=t, [sum(S)/length(S)[1] for S in Ss] #时间，平均值
    xlim,ylim=maximum(x0),maximum(y0)*1.1 #取值范围
    p20=plot(x0, y0,  markersize=markersize, markerstrokewidth=0, dpi=dpi, legend = :none,
        xlims=(0,xlim),ylims=(0,ylim), guidefont = guifont,tickfontsize=tickfontsize,
        framestyle=:box, ylabel=L"[Ca^{2+}]_{\rm SR}~(\mu M)",xlabel=L"t~(s)") #画图
    
    #节点演化图
    x0,y0,z0=t, [i/1000. for i in 1:Nnode],[S[i] for i in 1:Nnode, S in Ss]
    cmin,cmax=minimum(z0),maximum(z0)+1e-5

    p2=heatmap(x0, y0, z0, c = cg,xlims=(0.,maximum(x0)),ylims=(0.,maximum(y0)),clims=(cmin,cmax),
    framestyle=:box, guidefont = guifont,tickfontsize=tickfontsize,
    dpi=dpi,ylabel=L"n_{nodes}~(~\times 1000~)",xlabel=L"t~(s)")
 
    #读取数据
    t,Ss=load(paras_plot,"t")[1],load(paras_plot,"Ss")[3]
    #平均值    
    x0,y0=t, [sum(S)/length(S)[1] for S in Ss]
    xlim,ylim=maximum(x0),maximum(y0)*1.1
    p30=plot(x0, y0,  markersize=markersize, markerstrokewidth=0, dpi=dpi, legend = :none,
        xlims=(0,xlim),ylims=(0,ylim), guidefont = guifont,tickfontsize=tickfontsize,
        framestyle=:box, ylabel=L"P_O",xlabel=L"t~(s)")
          
    #节点演化图
    x0,y0,z0=t, [i/1000. for i in 1:Nnode],[S[i] for i in 1:Nnode, S in Ss]
    cmin,cmax=minimum(z0),maximum(z0)+1e-5

    p3=heatmap(x0, y0, z0, c = cg,xlims=(0.,maximum(x0)),ylims=(0.,maximum(y0)),clims=(cmin,cmax),
    framestyle=:box, guidefont = guifont,tickfontsize=tickfontsize,
    dpi=dpi,ylabel=L"n_{nodes}~(~\times 1000~)",xlabel=L"t~(s)")
    #绘制图像
    l = @layout[[[a b c d ];[a b c d ]]; [[e; f; g] [h; i;j]]]
    
    plot(p[1],p[2],p[3],p[4],p[5],p[6],p[7],p[8],p1,p2,p3,p10,p20,p30, layout = l,title = ["(a)" "(b)" "(c)" "(d)" "(e)" "(f)"],titleloc = :left, titlefont = 10,size=(1000,900),left_margin=6mm,right_margin=0mm, bottom_margin=2mm, top_margin=0mm, link=:all, dpi=dpi)
    
    Plots.savefig("results/"*paras*"_cluster.png")
    end

end

#----------------------------------------------------------------------------------------------------------
function plot_long()
    ENV["GKSwstype"] = "100" #不弹出图像窗口

    paras=load("results/name.jld","name") #读取需要画图的数据文件名字
    paras_plot="results/"*paras*"_plot.jld" #数据文件名字
    paras="Single0"
    @show paras

    dpi=1200 #画图精度
    markersize=0.5;linewidth=0.001; guifont=(8, :dark); tickfontsize=6 #点大小，线宽度，框架字体，刻度字体

     #第二三列第一个图
    #读取数据
    t,Ss=load(paras_plot,"t")[1],load(paras_plot,"Ss")[1] #时间，Ss中第一个数值

    Nnode=length(Ss[1])[1] 
    cg=cgrad([:white, :green,:black], [0, 0.5,0.75, 1.])
 
    #平均值，   
    x0,y0=t, [sum(S)/length(S)[1] for S in Ss] #时间，平均值
    xlim,ylim=maximum(x0),maximum(y0)*1.1 #取值范围
    xlims=(0,xlim)
    xticks0=0.02
    xticks0=xlim/10.
    p10=plot(x0,y0, markersize=markersize, markerstrokewidth=0, dpi=dpi,legend = :none,
        xlims=xlims,ylims=(0.1,ylim), guidefont = guifont,tickfontsize=tickfontsize,
        framestyle=:box, ylabel=L"[Ca^{2+}]~(\mu M)",xticks=[xticks0*i for i in 1:maximum(x0)/xticks0]) #画图

    x0,y0,z0=t,[i/1000. for i in 1:Nnode],[S[i] for i in 1:Nnode, S in Ss] #时间，节点序号，浓度

    cmin,cmax=minimum(z0),maximum(z0)+1e-20 #浓度范围
    p1=heatmap(x0, y0, z0, c = cg,xlims=xlims,ylims=(0.,maximum(y0)),clims=(cmin,cmax), 
        framestyle=:box, guidefont = guifont,tickfontsize=tickfontsize,
        dpi=dpi,ylabel=L"n_{nodes}~(~\times 1000~)") #画图
    
    #读取数据
    t,Ss=load(paras_plot,"t")[2],load(paras_plot,"Ss")[2] 

    x0,y0=t, [sum(S)/length(S)[1] for S in Ss] #时间，平均值
    xlim,ylim=maximum(x0),maximum(y0)*1.1 #取值范围
    p20=plot(x0, y0,  markersize=markersize, markerstrokewidth=0, dpi=dpi, legend = :none,
        xlims=xlims,ylims=(900,1010), guidefont = guifont,tickfontsize=tickfontsize,
        framestyle=:box, ylabel=L"[Ca^{2+}]_{\rm SR}~(\mu M)",xlabel=L"t~(s)",xticks=[xticks0*i for i in 1:maximum(x0)/xticks0]) #画图
    
     #节点演化图，第三列第二个图
    x0,y0,z0=t, [i/1000. for i in 1:Nnode],[S[i] for i in 1:Nnode, S in Ss]
    cmin,cmax=minimum(z0),maximum(z0)+1e-5

    p2=heatmap(x0, y0, z0, c = cg,xlims=xlims,ylims=(0.,maximum(y0)),clims=(cmin,cmax),
    framestyle=:box, guidefont = guifont,tickfontsize=tickfontsize,
    dpi=dpi,ylabel=L"n_{nodes}~(~\times 1000~)",xlabel=L"t~(s)")

    #读取数据
    t,Ss=load(paras_plot,"t")[1],load(paras_plot,"Ss")[3]

    x0,y0=t, [sum(S) for S in Ss]
    xlim,ylim=maximum(x0),maximum(y0)*1.1
    p30=plot(x0, y0,  markersize=markersize, markerstrokewidth=0, dpi=dpi, legend = :none,
        xlims=xlims,ylims=(0,ylim), guidefont = guifont,tickfontsize=tickfontsize,
        framestyle=:box, ylabel=L"P_O",xlabel=L"t~(s)",xticks=[xticks0*i for i in 1:maximum(x0)/xticks0])
          
     x0,y0,z0=t, [i/1000. for i in 1:Nnode],[S[i] for i in 1:Nnode, S in Ss]
    cmin,cmax=minimum(z0),maximum(z0)+1e-5

    p3=heatmap(x0, y0, z0, c = cg,xlims=xlims,ylims=(0.,maximum(y0)),clims=(cmin,cmax),
    framestyle=:box, guidefont = guifont,tickfontsize=tickfontsize,
    dpi=dpi,ylabel=L"n_{nodes}~(~\times 1000~)",xlabel=L"t~(s)")

    #绘制图像
    l = @layout [a;b;c]
    
    plot(p10,p20,p30, layout = l,title = ["(a)" "(b)" "(c)"],titleloc = :left, titlefont = 10,size=(1800,600),left_margin=3mm,right_margin=0mm, bottom_margin=2mm, top_margin=0mm, link=:all, dpi=dpi)
    
    Plots.savefig("results/"*paras*"_cluster0.png")

    plot(p1,p2,p3, layout = l,title = ["(a)" "(b)" "(c)"],titleloc = :left, titlefont = 10,size=(1800,600),left_margin=3mm,right_margin=0mm, bottom_margin=2mm, top_margin=0mm, link=:all, dpi=dpi)
    
    Plots.savefig("results/"*paras*"_cluster1.png")

end

#----------------------------------------------------------------------------------------------------------
function data()
    pointss=[]
    paras=load("results/name.jld","name") #读取需要画图的数据文件名字
    net=load("results/"*paras*"_net.jld","net") #读取cluster网络的数据
    net=NamedTuple{net[1]}(net[2])
    x,y,n=net.x,net.y,net.n
    pointss=[[x[i],y[i],n[i]] for i in eachindex(x)]
 
    i=1
    points_s = Vector{Vector{Float64}}[]
    temp = Vector{Float64}[]
    push!(temp,pointss[1])
    lpoints=length(pointss)
    while i < lpoints
        if pointss[i][3]==pointss[i+1][3]
            push!(temp,pointss[i+1])
            i += 1
        else
            push!(points_s,temp)
            temp =  Vector{Float64}[]
            i += 1
            push!(temp,pointss[i])
      end
    end
    push!(points_s,temp)
    points_s
end 
# RyR间距离
function sd(points_s,N)

    points_sd=Vector{Vector{Float64}}[] #整理好且筛选长度后的数组
    for points0 in points_s
        if length(points0) < N
            continue
        elseif length(points0) >=N
            push!(points_sd,points0)
    end
    end
    points_sd

end 
function d2d4()

    points_s=data()

    points_sd=sd(points_s,2)

    number = Float64[]
    for points in points_sd, n in range(1,length(points))
        min_dist = 999999
        for m in range(1,length(points))
            if points[n]!=points[m] && points[m][3]==points[n][3]
                temp_dist = sqrt((points[m][1]-points[n][1])^2+(points[m][2]-points[n][2])^2)
                if temp_dist < min_dist
                    min_dist = temp_dist
                end
            end
        end
        push!(number,min_dist)
    end
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    histogram(number, color="red", xlabel="distance(nm)",ylabel="% of RyRs", bins = range(0,100,length = 20),leg = false )
 
    points_sd=sd(points_s,5)

    #距离归类
    tempt=Float64[]#团簇中第i个点和其他所有点的距离[dis,...]
    temp_1 = Vector{Float64}[] #归纳一个团簇中每个点和其他所有点的距离[[dis1,...],[dis2,...]]
    temp_s = Vector{Vector{Float64}}[] #[[[dis1,..,],[dis2,..]]第一个团簇,[[dis1,...],[dis2,...]]第二个团簇]
    k=1
    while k <= length(points_sd)
        for n in range(1,length(points_sd[k]))
            for m in range(1,length(points_sd[k]))
                if points_sd[k][n]!=points_sd[k][m]
                    temp_dist = sqrt((points_sd[k][m][1]-points_sd[k][n][1])^2+(points_sd[k][m][2]-points_sd[k][n][2])^2)
                    push!(tempt,temp_dist)
                end
            end
            push!(temp_1,tempt)
            tempt=Float64[]
        end
        push!(temp_s,temp_1)
        temp_1 = Vector{Float64}[]
        k += 1
    end

    #选出最小的前四个求距离平均值
    number = Float64[]
    for t in temp_s
        for r in t
            d = sort(r)[1:4] #排序：每个团簇中第i个点到其他所有点的距离，然后只取前4个，求平均
            d_average=sum(d)/length(d)
            push!(number,d_average)
        end
    end
    histogram!(number, color="green", xlabel="distance(nm)",ylabel="% of RyRs", bins = range(0,100,length = 20),leg = false,alpha=0.5 )
    paras = load("results/name.jld", "name") #读取需要画图的数据文件名字
    Plots.pdf("results/"*paras*"_d2d4.pdf")
end

#----------------------------------------------------------------------------------------------------------
#画动态度
function plot_spark()
    ENV["GKSwstype"] = "100" #不弹出图像窗口

    paras = load("results/name.jld", "name") #读取需要画图的数据文件名字
    name1 = "results/" * paras * "_plot.jld" #数据文件名字

    net=load("results/"*paras*"_net.jld","net") #读取cluster网络的数据
    net=NamedTuple{net[1]}(net[2])
    N,size,x,y,A=net.N,net.sizeV,net.x,net.y,net.A
    x = x ./ 1000.0
    y = y ./ 1000.0 #将坐标缩小1000倍
    NC = length(size) #number of cluster
    dpi = 500 #画图精度
    markersize = 2
    linewidth = 0.001
    guifont = (8, :dark)
    tickfontsize = 6 #点大小，线宽度，框架字体，刻度字体

    #设置画图格式，空图
    pnet = plot(guidefont=guifont, framestyle=:box, ylabel=L"y~(\mu m)", xlabel=L"x~(\mu m)", dpi=dpi,
        xlims=(0.0, maximum(x)), ylims=(0.0, maximum(y)), tickfontsize=tickfontsize, legend=:none, size=(600, 600))
    #画连线
    xx, yy, net0 = [], [], []
    for i in 1:N, j in i:N
        if A[i, j] > 100
            push!(xx, [x[i], x[j]])
            push!(yy, [y[i], y[j]])
            push!(net0, A[i, j])
        end
    end
    plot!(pnet, xx, yy, legend=:none, linecolor=:grey, linewidth=transpose(net0) * linewidth) #链接
    #画点，不同团簇颜色不同
    NN = 1
    for i in 1:NC
        scatter!(pnet, x[NN:NN+size[i]-1], y[NN:NN+size[i]-1], markercolor=i, markersize=markersize, markershape=:rect, markerstrokewidth=0) #节点
        NN = NN + size[i]
    end
    t, Ss = load(name1, "t")[1], load(name1, "Ss")[3]
    Nnode = length(Ss[1])[1]
    ly = maximum(y) * 0.95
    lx = maximum(x) * 0.05
    anim = @animate for i in 1:100
        S = Ss[i]
        xo, yo = [], []
        for inode = 1:Nnode
            if S[inode] == 1
                push!(xo, x[inode])
                push!(yo, y[inode])
                #println(xo,yo)
            end
        end

        scatter(pnet, xo, yo, annotations=(lx, ly, Plots.text("t=" * string(format(t[i], precision=4)) * "s", 12, :left)), color=:red, markersize=7, markerstrokecolor=:yellow, markerstrokewidth=1, markershape=:star7)
    end
 

    @time gif(anim, "results/"*paras*"_spark.gif", fps=4)

end

#----------------------------------------------------------------------------------------------------------
function plot_net()
    ENV["GKSwstype"] = "100" #不弹出图像窗口



    paras = load("results/name.jld", "name") #读取需要画图的数据文件名字
    net=load("results/"*paras*"_net.jld","net") #读取cluster网络的数据
    net=NamedTuple{net[1]}(net[2])
    N,size,x,y,A=net.N,net.sizeV,net.x,net.y,net.A
    x=x./1000.;y=y./1000. #将坐标缩小1000倍
    NC=length(size) #number of cluster
    dpi=1200 #画图精度
    markersize=0.5;linewidth=0.0001; guifont=(12, :dark); tickfontsize=10 #点大小，线宽度，框架字体，刻度字体

    #画cluster网络图 
    #设置画图格式，空图
    pnet=plot(guidefont = guifont,framestyle=:box,dpi=dpi,xlims=(0.,maximum(x)),ylims=(0.,maximum(y)),tickfontsize=tickfontsize-2,legend = :none, xticks=0:0.5: maximum(x), yticks=0:0.5: maximum(x))
    #画连线
    Amax=maximum(A) 
    lcolor=(Amax .- A)/Amax
    for i in 1:N, j in i:N
        if A[i,j]>5000   
             plot!(pnet,[x[i],x[j]],[y[i],y[j]],legend = :none,linecolor=RGB(lcolor[i,j], lcolor[i,j], lcolor[i,j]),linewidth=A[i,j]*linewidth) #链接
        end
    end

    
    Plots.savefig("results/"*paras*"_net.pdf")



end

#----------------------------------------------------------------------------------------------------------
@time plot_res()
#@time plot_long()
#@time d2d4()
#@time plot_spark()
#@time plot_net()
