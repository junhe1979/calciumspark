using DelimitedFiles, Plots, LaTeXStrings, Colors, Compose, JLD2, Formatting, ProgressBars, Statistics

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
function readdata(filename)
    # 打开文件，逐行读取
    lines = readlines(filename)

    # 初始化一个空数组，用于存储每列的数组
    columns = []

    # 逐行处理文件内容
    for line in lines
        # 将每行数据按空格或其他分隔符分割
        values = split(line)

        # 如果当前列数比 columns 中的数组数多，动态添加新列
        for i in eachindex(values)
            if length(columns) < i
                push!(columns, Float64[])  # 创建新列数组
            end
            # 将当前行的第 i 列数据加入对应的列数组中
            push!(columns[i], parse(Float64, values[i]))
        end
    end

    return columns

end

function NpeakDuration(paras, t, ryr, CaSS, No, tp)


    Casp = [] #把所有钙信号钙闪烁钙火花都筛选出来分别形成数组在统一放在这个数组里
    CRE = []
    length0 = length(t)
    for i in 1:length0
        inext = i + 1
        if i == length0 && ryr[i] == 0
            continue
        elseif i == length0 && ryr[i] != 0
            push!(CRE, [t[i], ryr[i], CaSS[i]])
            push!(Casp, CRE)
            CRE = []
        elseif ryr[i] == 0 && ryr[inext] == 0
            continue
        elseif ryr[inext] != 0
            push!(CRE, [t[i], ryr[i], CaSS[i]])
        elseif ryr[i] != 0 && ryr[inext] == 0
            push!(CRE, [t[i], ryr[i], CaSS[i]])
            push!(Casp, CRE)
            CRE = []
        end
    end


    #把峰值相同的归类放进一个数组中
    Npeak = [] #把所有钙信号钙闪烁钙火花的峰值整理到这个数组
    CaSSpeak = [] #把所有钙信号钙闪烁钙火花的峰值整理到这个数组
    Duration = [] # 所有钙信号钙闪烁钙火花的时间 
    TTS = [] # 到2/3的时间 
    for Ca in Casp
        RyRO = [Ca0[2] for Ca0 in Ca] #open数量数组
        CaSS0 = [Ca0[3] for Ca0 in Ca] #时间数组
        time = [Ca0[1] for Ca0 in Ca] #时间数组
        pmax = findmax(RyRO) #open 峰值
        CaSSmax = findmax(CaSS0) #open 峰值
        #@show time
        TTS0 = 0
        for i in eachindex(RyRO)
            if RyRO[i] > pmax[1] * 2.0 / 3.0
                TTS0 = time[i] #2/3的时间
                break
            end
        end
        time0 = time[1]
        if No == 0
            lNo = (time0 < tp * 1.5)
        elseif No == 1
            lNo = (time0 < tp * 0.5)
        elseif No == 2
            lNo = (tp * 0.5 < time0 < tp * 1.5)
        end
        if lNo
            if SubString(paras, 1, 3) == "LCC"
                time0 = time[1] - mod(time[1] + 1e-10, 0.2) + 1e-10 #对LCC，第一个
            end
            TTS00 = (TTS0 - time0) * 10^3
            t = (last(time) - time0) * 10^3 #ms
            if pmax[1] < 5 && mod(time[1] + 1e-10, 0.2) > 8e-3 && SubString(paras, 1, 3) == "LCC"

            else
                #@show pmax[1], t
                push!(TTS, TTS00)
                push!(Npeak, pmax[1]) #保存峰值
                push!(CaSSpeak, CaSSmax[1]) #保存峰值
                push!(Duration, t)
            end
        end
    end


    return Npeak, Duration, TTS, CaSSpeak

end

function CRE(paras; No=0, tp=0.01, CaSS=false)
    f = open("results/" * paras * "_res.txt", "r")
    lines = readlines(f)
    close(f)

    str1, str2, str3 = split(lines[1], ",")
    Ndistribution = parse(Int64, str1)
    Nrepeat = parse(Int64, str2)
    datapoints = parse(Int64, str3)

    NpeakT = []
    DurationT = []
    TTST = []
    CaSSpeakT = []
    for idistribution in 1:Ndistribution
        NCluster = parse(Int64, lines[(idistribution-1)*Nrepeat*datapoints+idistribution+1])

        for irepeat in 1:Nrepeat
            n = (idistribution - 1) * Nrepeat * datapoints + (irepeat - 1) * datapoints + idistribution
            str = split.(lines[n+2:n+datapoints+1], ",")
            t = [parse(Float64, str0[1]) for str0 in str]
            for i in 1:NCluster
                #@show str[1]
                #@show str[1][i*2]
                ryr = [parse(Int64, str0[i*2]) for str0 in str]
                if CaSS == true
                    CaSS0 = [parse(Float64, str0[i*2+1]) for str0 in str]
                else
                    CaSS0 = [1.0 for str0 in str]
                end
                Npeak, Duration, TTS, CaSSpeak = NpeakDuration(paras, t, ryr, CaSS0, No, tp)
                append!(NpeakT, Npeak)
                append!(DurationT, Duration)
                append!(TTST, TTS)
                append!(CaSSpeakT, CaSSpeak)
            end
        end

    end
    return NpeakT, DurationT, TTST, CaSSpeakT
end

function plotSponE()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    paras = "SponE"
    paras_plot = "results/" * paras * "30_plot.jld" #数据文件名字

    dpi = 1200 #画图精度
    markersize = 0.5
    guifont = (11, :dark)
    tickfontsize = 9 #点大小，线宽度，框架字体，刻度字体

    #第二三列第一个图
    #读取数据
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[1] #时间，Ss中第一个数值
    x0, y0 = t, [sum(S) / length(S)[1] for S in Ss] #时间，平均值
    xlims0 = (0, 10)
    xticks0 = 0.0:1:10

    xlims = (5.698, 5.81)
    xticks = 5.7:0.03:5.8
    p10 = plot(x0, y0, markersize=markersize, markerstrokewidth=0, dpi=dpi, grid=:none,
        legend=:left, foreground_color_legend=nothing,
        legendfont=Plots.Font("sans-serif", 9, :right, :vcenter, 0.0, RGB(0.0, 0.0, 0.0)),
        xlims=xlims0, xticks=xticks0, ylims=(0.1, 1000), guidefont=guifont,
        tickfontsize=tickfontsize, framestyle=:box, xlabel=L"t~(s)",
        ylabel=L"[{\rm Ca}^{2+}]~(\mu \rm M)", label=L"[\textrm{ Ca}^{2+}]^\textrm{SS}")

    plot!(p10, x0, y0, xlims=xlims, xticks=xticks, ylims=(0, 1000), yticks=0:200:1000,
        tickfontsize=7, grid=:none,
        inset=(1, bbox(0.69, 0.25, 0.28, 0.5)),
        subplot=2,
        bg_inside=nothing,
        legend=:none,
        framestyle=:box)



    t, Ss = load(paras_plot, "t")[2], load(paras_plot, "Ss")[2]
    x0, y0 = t, [sum(S) / length(S)[1] for S in Ss] #时间，平均值
    plot!(p10, x0, y0, label=L"[\textrm{\rm Ca}^{2+}]^\textrm{ JSR}") #画图

    plot!(p10, x0, y0, xlims=xlims, xticks=xticks, ylims=(0, 1000), yticks=0:200:1000,
        tickfontsize=7, grid=:none, linecolor=:red,
        inset=(1, bbox(0.69, 0.25, 0.28, 0.5)),
        subplot=3,
        bg_inside=nothing,
        legend=:none,
        framestyle=:box)


    #读取数据
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[3]
    x0, y0 = t, [sum(S) for S in Ss]
    p20 = plot(x0, y0, markersize=markersize, markerstrokewidth=0, dpi=dpi, grid=:none,
        color=:black, legend=:none, xlims=xlims0, xticks=xticks0, ylims=(0, 50),
        guidefont=guifont, tickfontsize=tickfontsize, framestyle=:box,
        ylabel=L"N_{\rm O}", xlabel=L"t~(s)")

    plot!(p20, x0, y0, xlims=xlims, xticks=xticks, ylims=(0, 50), yticks=0:10:50,
        tickfontsize=7, grid=:none, color=:black,
        inset=(1, bbox(0.69, 0.25, 0.28, 0.5)),
        subplot=2,
        bg_inside=nothing,
        legend=:none,
        framestyle=:box
    )


    net = load("results/" * paras * "_net.jld", "net") #读取cluster网络的数据
    N, size, x, y, A = net.N, net.sizeV, net.x, net.y, net.A
    x = x ./ 1000.0 .- 0.2
    y = y ./ 1000.0 .- 0.9 #将坐标缩小1000倍
    xlims = (0.0, 0.6)
    ylims = (0.0, 0.5)
    inset = (1, bbox(0.06, 0.06, 0.2, 0.7))
    NC = length(size) #number of cluster
    dpi = 1200 #画图精度
    markersize = 1.5
    linewidth = 0.00001
    guifont = (8, :dark)
    tickfontsize = 8 #点大小，线宽度，框架字体，刻度字体

    #画cluster网络图 
    Amax = maximum(A)
    lcolor = (Amax .- A) / Amax
    Aij = [A[i, j] for i in 1:N for j in i:N]
    iv = [i for i in 1:N for j in i:N]
    jv = [j for i in 1:N for j in i:N]
    sort = sortperm(Aij)

    isp = 3
    for s in sort
        i = iv[s]
        j = jv[s]
        if A[i, j] > 1000
            plot!([x[i], x[j]], [y[i], y[j]], xlims=xlims, ylims=ylims, legend=:none, linecolor=RGB(lcolor[i, j], lcolor[i, j], lcolor[i, j]), linewidth=A[i, j] * linewidth, inset=inset, grid=:none, subplot=isp, tickfontsize=8, bg_inside=nothing, showaxis=false)
            isp += 1

        end
    end

    #println(lmargin)
    xx, yy = [], []
    for i in eachindex(net.lmargin)
        if net.lmargin[i] < 0
            push!(xx, x[i])
            push!(yy, y[i])
        end
    end
    scatter!(xx, yy, xlims=xlims, ylims=ylims, color=:blue, legend=:none, markersize=markersize, grid=:none, tickfontsize=6, markerstrokecolor=:black, markerstrokewidth=1, markershape=:circle, inset=inset, subplot=isp, bg_inside=nothing, showaxis=false)


    i = convert(Int64, ceil(length(t) / 10.0 * 5.71))
    println(i)
    S = Ss[i]
    xo, yo = [], []
    Nnode = length(S)[1]
    for inode = 1:Nnode
        if S[inode] == 1
            push!(xo, x[inode])
            push!(yo, y[inode])
            #println(xo,yo)
        end
    end


    scatter!(x, y, xlims=xlims, ylims=ylims, legend=:none, tickfontsize=6, grid=:none,
        inset=inset, subplot=isp + 1, bg_inside=nothing, showaxis=false, markercolor=:blue,
        markersize=markersize, markerstrokewidth=0)


    scatter!(xo, yo, xlims=xlims, ylims=ylims, color=:blue, legend=:none, grid=:none,
        tickfontsize=6, inset=inset, markercolor=:red, markersize=4, markerstrokecolor=:yellow,
        markerstrokewidth=0.5, markershape=:star7, subplot=isp + 2, bg_inside=nothing,
        showaxis=:x, xticks=0:0.2:0.49, annotations=(0.468, -0.048,
            Plots.text(L"(\mu\textrm{m})", 7, :left)), alpha=0.7) #节点



    l = @layout[a{0.99w}; a{0.99w}]
    title = fill("", 1, isp + 5)
    title[1] = "(a) Spontaneous"
    title[4] = "(b)"


    plot(p10, p20, layout=l, title=title, titleloc=:left, titlefont=11, size=(500, 400), left_margin=0mm, right_margin=0mm, bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)


    Plots.savefig("results/SponE.pdf")
end

function plotLCCE()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    paras = "LCCE"
    paras_plot = "results/" * paras * "1_plot.jld" #数据文件名字

    dpi = 1200 #画图精度
    markersize = 0.5
    guifont = (11, :dark)
    tickfontsize = 9 #点大小，线宽度，框架字体，刻度字体

    #读取数据
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[1] #时间，Ss中第一个数值
    x0, y0 = t, [sum(S) / length(S)[1] for S in Ss] #时间，平均值
    xlims0 = (0, 10)
    xticks0 = 0.0:1:10

    xlims = (5.998, 6.11)
    xticks = 6.0:0.03:6.14
    p10 = plot(x0, y0, markersize=markersize, markerstrokewidth=0, dpi=dpi, grid=:none,
        legend=:left, foreground_color_legend=nothing,
        legendfont=Plots.Font("sans-serif", 9, :right, :vcenter, 0.0, RGB(0.0, 0.0, 0.0)),
        xlims=xlims0, xticks=xticks0, ylims=(0.1, 1000), guidefont=guifont,
        tickfontsize=tickfontsize, framestyle=:box, xlabel=L"t~(s)",
        ylabel=L"[{\rm Ca}^{2+}]~(\mu \rm M)", label=L"[\textrm{ Ca}^{2+}]^\textrm{SS}")

    plot!(p10, x0, y0, xlims=xlims, xticks=xticks, ylims=(0, 1000), yticks=0:200:1000,
        tickfontsize=7, grid=:none,
        inset=(1, bbox(0.69, 0.25, 0.28, 0.5)),
        subplot=2,
        #bg_inside = nothing,
        legend=:none,
        framestyle=:box)

    t, Ss = load(paras_plot, "t")[2], load(paras_plot, "Ss")[2]
    x0, y0 = t, [sum(S) / length(S)[1] for S in Ss] #时间，平均值
    plot!(p10, x0, y0, label=L"[\textrm{\rm Ca}^{2+}]^\textrm{ JSR}") #画图


    plot!(p10, x0, y0, xlims=xlims, xticks=xticks, ylims=(0, 1000), yticks=0:200:1000,
        tickfontsize=7, grid=:none, linecolor=:red,
        inset=(1, bbox(0.69, 0.25, 0.28, 0.5)),
        subplot=3,
        bg_inside=nothing,
        legend=:none,
        framestyle=:box)

    #读取数据
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[3]
    x0, y0 = t, [sum(S) for S in Ss]
    p20 = plot(x0, y0, markersize=markersize, markerstrokewidth=0, dpi=dpi, grid=:none,
        color=:black, legend=:none, xlims=xlims0, xticks=xticks0, ylims=(0, 50),
        guidefont=guifont, tickfontsize=tickfontsize, framestyle=:box,
        ylabel=L"N_{\rm O}", xlabel=L"t~(s)")

    plot!(p20, x0, y0, xlims=xlims, xticks=xticks, ylims=(0, 50), yticks=0:10:50,
        tickfontsize=7, grid=:none, color=:black,
        inset=(1, bbox(0.69, 0.25, 0.28, 0.5)),
        subplot=2,
        #bg_inside = nothing,
        legend=:none,
        framestyle=:box
    )


    net = load("results/" * paras * "_net.jld", "net") #读取cluster网络的数据
    N, size, x, y, A = net.N, net.sizeV, net.x, net.y, net.A
    x = x ./ 1000.0 .- 0.2
    y = y ./ 1000.0 .- 0.9 #将坐标缩小1000倍
    xlims = (0.0, 0.6)
    ylims = (0.0, 0.5)
    inset = (1, bbox(0.06, 0.06, 0.2, 0.7))
    NC = length(size) #number of cluster
    dpi = 1200 #画图精度
    markersize = 1.5
    linewidth = 0.00001
    guifont = (8, :dark)
    tickfontsize = 8 #点大小，线宽度，框架字体，刻度字体

    #画cluster网络图 
    Amax = maximum(A)
    lcolor = (Amax .- A) / Amax
    Aij = [A[i, j] for i in 1:N for j in i:N]
    iv = [i for i in 1:N for j in i:N]
    jv = [j for i in 1:N for j in i:N]
    sort = sortperm(Aij)

    isp = 3
    for s in sort
        i = iv[s]
        j = jv[s]
        if A[i, j] > 1000
            plot!([x[i], x[j]], [y[i], y[j]], xlims=xlims, ylims=ylims, legend=:none, linecolor=RGB(lcolor[i, j], lcolor[i, j], lcolor[i, j]), linewidth=A[i, j] * linewidth, inset=inset, grid=:none, subplot=isp, tickfontsize=8, bg_inside=nothing, showaxis=false)
            isp += 1

        end
    end

    #println(lmargin)
    xx, yy = [], []
    for i in eachindex(net.lmargin)
        if net.lmargin[i] < 0
            push!(xx, x[i])
            push!(yy, y[i])
        end
    end
    scatter!(xx, yy, xlims=xlims, ylims=ylims, color=:blue, legend=:none, markersize=markersize, grid=:none, tickfontsize=6, markerstrokecolor=:black, markerstrokewidth=1, markershape=:circle, inset=inset, subplot=isp, bg_inside=nothing, showaxis=false)
    i = convert(Int64, ceil(length(t) / 10.0 * 6.01))
    println(i)
    S = Ss[i]
    xo, yo = [], []
    Nnode = length(S)[1]
    for inode = 1:Nnode
        if S[inode] == 1
            push!(xo, x[inode])
            push!(yo, y[inode])
            #println(xo,yo)
        end
    end


    scatter!(x, y, xlims=xlims, ylims=ylims, legend=:none, tickfontsize=6, grid=:none,
        inset=inset, subplot=isp + 1, bg_inside=nothing, showaxis=false, markercolor=:blue,
        markersize=markersize, markerstrokewidth=0)


    scatter!(xo, yo, xlims=xlims, ylims=ylims, color=:blue, legend=:none, grid=:none,
        tickfontsize=6, inset=inset, markercolor=:red, markersize=4, markerstrokecolor=:yellow,
        markerstrokewidth=0.5, markershape=:star7, subplot=isp + 2, bg_inside=nothing,
        showaxis=:x, xticks=0:0.2:0.49, annotations=(0.468, -0.048,
            Plots.text(L"(\mu\textrm{m})", 7, :left)), alpha=0.7) #节点



    l = @layout[a{0.99w}; a{0.99w}]
    title = fill("", 1, isp + 5)
    title[1] = "(c) LCC I"
    title[4] = "(d)"


    plot(p10, p20, layout=l, title=title, titleloc=:left, titlefont=11, size=(500, 400), left_margin=0mm, right_margin=0mm, bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)


    Plots.savefig("results/LCCE.pdf")
end

function plotSingle()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度
    guifont = (14, :dark)
    tickfontsize = 11
    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []

    for i in ["80", "50", "30", "M"]
        for Tr in ["Spon", "LCCI", "LCCII"]
            if i == "M"
                xlabel = L"\tau\ (\rm ms)"
                Nmax = 10
            else
                xlabel = ""
                Nmax = parse(Int64, i) / 4.0
            end

            if Tr in ["Spon"]
                ylabel = L"N^{\rm peak}_{\rm O}"
            else
                ylabel = ""
            end

            if Tr == "LCCII"
                colorbar = :legend
            else
                colorbar = :none
            end

            @show Tr * i

            Npeak, tCRE, TTS = CRE(Tr * i)


            println(Tr * i * "Npeak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(Tr * i * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(Tr * i * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            tlim = 60
            Nlim = 80
            histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
                xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.3],
                    [67], [Plots.text(L"N_\textrm{RyR}=%$i", 12)]), xlabel=xlabel,
                guidefont=guifont, framestyle=:box, grid=:none, ylabel=ylabel,
                dpi=dpi, colorbar_scale=:log10, colorbar=colorbar)
            #scatter!([mean(tCRE[Npeak.>Nmax])], [mean(Npeak[Npeak.>Nmax])], xerror=[std(tCRE[Npeak.>Nmax])], 
            #    yerror=[std(Npeak[Npeak.>Nmax])], markershape=:circle,markersize=4, markercolor=:black, 
            #    markerstrokewidth=0.5,markerstrokecolor=:black)
            t0 = plot!([0, tlim], [Nmax, Nmax], guidefont=guifont, framestyle=:box, grid=:none,
                tickfontsize=tickfontsize, color=:grey, legend=:none)

            push!(t, t0)

        end
    end


    l = @layout[[a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]]

    plot(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], t[11], t[12], layout=l,
        title=["(a) Spontaneous" "(b) LCC I" "(c) LCC II" "(d) Spontaneous" "(e) LCC I" "(f) LCC II" "(g) Spontaneous" "(h) LCC I" "(i) LCC II" "(j) Spontaneous" "(k) LCC I" "(l) LCC II"],
        titleloc=:left, titlefont=14, size=(800, 950), left_margin=3mm, right_margin=0mm,
        bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotSingle.pdf")
end

function plotCa0()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度

    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []
    Nmax = 10
    for N in ["1000", "5000", "5000s"]
        for C in ["", "nbf"]

            if C == "nbf"
                xlabel = L"\tau\ (\rm ms)"
            else
                xlabel = ""
            end
            if N == "1000"
                ylabel = L"N^{\rm peak}_{\rm O}"
            else
                ylabel = ""
            end

            if N == "5000s"
                if C == "nbf"
                    Nlabel = L"1800"
                else
                    Nlabel = L"5000"
                end
                Nlabels = L"~(\textrm{fix})"
                colorbar = :legend
            else
                Nlabel = latexstring(N)
                Nlabels = ""
                colorbar = :none
            end
            if N == "5000" && C == "nbf"
                Nlabel = L"1800"
            end


            name = "CSQ400C" * C * N
            #Cint = Int64(parse(Int64, C))
            nameN = L"\textrm{NSR}: " * Nlabel * Nlabels


            Npeak, tCRE, TTS = CRE(name)


            println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            tlim = 110
            Nlim = 60
            t0 = histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
                xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.05],
                    [60 * 0.8375], [Plots.text(nameN, 12, :left)]), xlabel=xlabel,
                guidefont=(18, :dark), framestyle=:box, grid=:none, ylabel=ylabel, tickfontsize=11,
                dpi=dpi, colorbar_scale=:log10, colorbar=colorbar)


            push!(t, t0)



        end
    end

    C2 = " Buffer function"
    C1 = " Binding to RyR"
    l = @layout[[a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]]
    plot(t[1], t[3], t[5], t[2], t[4], t[6], layout=l,
        title=["(a)" * C1 "(b)" * C1 "(c)" * C1 "(d)" * C2 "(e)" * C2 "(f)" * C2],
        titleloc=:left, titlefont=12, size=(800, 450), left_margin=3mm, right_margin=0mm,
        bottom_margin=4mm, top_margin=0mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotCa0.pdf")
end

function plotrate()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度

    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []
    Nmax = 10
    for N in ["C", "S", "K"]
        for C in ["rat", "sheep"]

            if C == "nbf"
                xlabel = L"\tau\ (\rm ms)"
            else
                xlabel = ""
            end
            if N == "1000"
                ylabel = L"N^{\rm peak}_{\rm O}"
            else
                ylabel = ""
            end

            if N == "K"
                if C == "nbf"
                    Nlabel = L"1800"
                else
                    Nlabel = L"5000"
                end
                Nlabels = L"~(\textrm{fix})"
                colorbar = :legend
            else
                Nlabel = latexstring(N)
                Nlabels = ""
                colorbar = :none
            end
            if N == "5000" && C == "nbf"
                Nlabel = L"1800"
            end


            name = C * N
            #Cint = Int64(parse(Int64, C))
            nameN = L"\textrm{NSR}: " * Nlabel * Nlabels


            Npeak, tCRE, TTS = CRE(name)


            println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            tlim = 60
            Nlim = 80
            t0 = histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
                xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.05],
                    [60 * 0.8375], [Plots.text(nameN, 12, :left)]), xlabel=xlabel,
                guidefont=(18, :dark), framestyle=:box, grid=:none, ylabel=ylabel, tickfontsize=11,
                dpi=dpi, colorbar_scale=:log10, colorbar=colorbar)


            push!(t, t0)



        end
    end

    C2 = " Buffer function"
    C1 = " Binding to RyR"
    l = @layout[[a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]]
    plot(t[1], t[3], t[5], t[2], t[4], t[6], layout=l,
        title=["(a)" * C1 "(b)" * C1 "(c)" * C1 "(d)" * C2 "(e)" * C2 "(f)" * C2],
        titleloc=:left, titlefont=12, size=(800, 450), left_margin=3mm, right_margin=0mm,
        bottom_margin=4mm, top_margin=0mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotrate.pdf")
end

function plotrateline()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度
    guifont = (11, :dark)
    tickfontsize = 14
    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    rate, Nt, Np, Npe, tau, taue = [], [], [], [], [], []
    Nmax = 10
    for N in ["rate_02", "rate_04", "rate_06", "rate_08", "rate_10", "rate_12", "rate_14", "rate_16", "rate_18", "rate_20", "rate_22", "rate_24", "rate_26", "rate_28", "rate_30", "rate_32", "rate_34", "rate_36", "rate_38"]
        #for N in ["rate_02", "rate_04", "rate_06", "rate_38", "rate_40"]
        name = N
        Npeak, tCRE, TTS = CRE(name)
        if length(Npeak[Npeak.>Nmax]) > 0
            println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            rate0 = parse(Float64, split(N, "_")[2]) / 10.0
            push!(rate, rate0)
            push!(Nt, length(Npeak[Npeak.>Nmax]))
            push!(Np, mean(Npeak[Npeak.>Nmax]))
            push!(Npe, std(Npeak[Npeak.>Nmax]))
            push!(tau, mean(tCRE[Npeak.>Nmax]))
            push!(taue, std(tCRE[Npeak.>Nmax]))
        end
    end
    a = scatter(rate, Nt ./ 1000.0,
        ylabel=L"N_\textrm{Event}\ (10^3)", xlabel=L"k_\textrm{O} (k_\textrm{O}^\textrm{rat})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(16, "Arial", :dark), yguidefont=(16, "Arial", :dark), xticks=0:1:4, yticks=0:0.2:2,
        xlim=(0, 4.0), ylim=(0, 1.1), legend=:none)

    b = scatter(rate, Np, yerror=Npe,
        ylabel=L"N_\textrm{O}^\textrm{peak}", xlabel=L"k_\textrm{O} (k_\textrm{O}^\textrm{rat})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(16, "Arial", :dark), yguidefont=(16, "Arial", :dark), xticks=0:1:4,
        xlim=(0, 4.0), ylim=(0, 50), legend=:none)

    c = scatter(rate, tau, yerror=taue,
        ylabel=L"\tau\ \textrm{(s)}", xlabel=L"k_\textrm{O} (k_\textrm{O}^\textrm{rat})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(16, "Arial", :dark), yguidefont=(16, "Arial", :dark), xticks=0:1:4,
        xlim=(0, 4.0), ylim=(0, 50), legend=:none)

    l = @layout[a b c]
    plot(a, b, c, layout=l,
        title=["(a)" "(b)" "(c)"],
        titleloc=:left, titlefont=14, size=(800, 250), left_margin=5mm, right_margin=0mm,
        bottom_margin=8mm, top_margin=1mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotrateline.pdf")

end

function plotNSRmyo()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度
    guifont = (11, :dark)
    tickfontsize = 13
    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    rate, Nt, Np, Npe, tau, taue = [], [], [], [], [], []
    Nmax = 10
    #for N in ["NSR_07", "NSR_14"]
    for N in ["NSR_06", "NSR_07", "NSR_08", "NSR_09", "NSR_10", "NSR_11", "NSR_12", "NSR_13", "NSR_14"]
        name = N
        Npeak, tCRE, TTS = CRE(name)
        if length(Npeak[Npeak.>Nmax]) > 0
            println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            rate0 = parse(Float64, split(N, "_")[2]) / 10.0
            push!(rate, rate0)
            push!(Nt, length(Npeak[Npeak.>Nmax]))
            push!(Np, mean(Npeak[Npeak.>Nmax]))
            push!(Npe, std(Npeak[Npeak.>Nmax]))
            push!(tau, mean(tCRE[Npeak.>Nmax]))
            push!(taue, std(tCRE[Npeak.>Nmax]))
        end
    end
    myo, mNt, mNp, mNpe, mtau, mtaue = [], [], [], [], [], []
    #for N in ["myo_01", "myo_12"]
    for N in ["myo_00", "myo_01", "myo_015", "myo_02", "myo_03", "myo_04", "myo_05", "myo_06", "myo_07", "myo_08", "myo_09", "myo_10", "myo_11", "myo_12", "myo_14"]
        name = N
        Npeak, tCRE, TTS = CRE(name)
        if length(Npeak[Npeak.>Nmax]) > 0
            println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            rate0 = parse(Float64, split(N, "_")[2]) / 10.0
            if N == "myo_015"
                rate0 = parse(Float64, split(N, "_")[2]) / 100.0
            end
            push!(myo, rate0)
            push!(mNt, length(Npeak[Npeak.>Nmax]))
            push!(mNp, mean(Npeak[Npeak.>Nmax]))
            push!(mNpe, std(Npeak[Npeak.>Nmax]))
            push!(mtau, mean(tCRE[Npeak.>Nmax]))
            push!(mtaue, std(tCRE[Npeak.>Nmax]))
        end
    end

    a = scatter(rate, Nt ./ 1000.0,
        ylabel=L"N_\textrm{Event}\ (10^3)", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{NSR} (10^3\mu\textrm{M})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(14, "Arial", :dark), yguidefont=(14, "Arial", :dark), xticks=0:0.2:2, yticks=0:0.2:2,
        xlim=(0.6, 1.42), ylim=(0, 1.0), legend=:none)

    b = scatter(rate, Np, yerror=Npe,
        ylabel=L"N_\textrm{O}^\textrm{peak}", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{NSR} (10^3\mu\textrm{M})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(14, "Arial", :dark), yguidefont=(14, "Arial", :dark), xticks=0:0.2:2,
        xlim=(0.6, 1.42), ylim=(0, 50), legend=:none)

    scatter(rate, tau, yerror=taue,
        ylabel=L"\tau\ \textrm{(s)}", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{NSR} (10^3\mu\textrm{M})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(14, "Arial", :dark), yguidefont=(14, "Arial", :dark), xticks=0:0.2:2,
        xlim=(0.6, 1.42), ylim=(0, 50), legend=:none)
    d = readdata("results/Bers.txt")
    c = plot!(d[1]./1000., d[2] )

    d = scatter(myo, mNt ./ 1000.0,
        ylabel=L"N_\textrm{Event}\ (10^3)", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{myo} (10^3\mu\textrm{M})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(14, "Arial", :dark), yguidefont=(14, "Arial", :dark), xticks=0:0.4:2, yticks=0:2:12,
        xlim=(-0.1, 1.3), ylim=(0, 11), legend=:none)

    e = scatter(myo, mNp, yerror=mNpe,
        ylabel=L"N_\textrm{O}^\textrm{peak}", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{myo} (10^3\mu\textrm{M})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(14, "Arial", :dark), yguidefont=(14, "Arial", :dark), xticks=0:0.4:2,
        xlim=(-0.1, 1.3), ylim=(0, 50), legend=:none)

    f = scatter(myo, mtau, yerror=mtaue,
        ylabel=L"\tau\ \textrm{(s)}", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{myo} (10^3\mu\textrm{M})",
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red, guifont=guifont,
        tickfontsize=tickfontsize, xguidefont=(14, "Arial", :dark), yguidefont=(14, "Arial", :dark), xticks=0:0.4:2,
        xlim=(-0.1, 1.3), ylim=(0, 70), legend=:none)

    l = @layout[a b c; d e f]
    plot(a, b, c, d, e, f, layout=l,
        title=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)"],
        titleloc=:left, titlefont=14, size=(800, 500), left_margin=3mm, right_margin=3mm,
        bottom_margin=4mm, top_margin=0mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotNSRmyo.pdf")

end

function plotCa()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度

    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []
    Nmax = 10
    for N in ["500", "1000", "1500"]
        for C in ["0", "100", "200"]

            if N == "1500"
                xlabel = L"\tau\ (\rm ms)"
            else
                xlabel = ""
            end
            if C in ["0"]
                ylabel = L"N^{\rm peak}_{\rm O}"
            else
                ylabel = ""
            end

            if C == "200"
                colorbar = :legend
            else
                colorbar = :none
            end



            name = "CaN" * N * "C" * C
            Nint = parse(Int64, N)
            Cint = Int64(parse(Int64, C))
            nameN = L"\textrm{NSR}: %$Nint, \textrm{myo}: %$Cint"


            Npeak, tCRE, TTS = CRE(name)


            println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
            println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
            println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
            tlim = 110
            Nlim = 60
            t0 = histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
                xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.05],
                    [60 * 0.8375], [Plots.text(nameN, 11, :left)]), xlabel=xlabel,
                guidefont=(16, :dark), framestyle=:box, grid=:none, ylabel=ylabel, tickfontsize=11,
                dpi=dpi, colorbar_scale=:log10, colorbar=colorbar)


            push!(t, t0)



        end
    end

    l = @layout[[a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]]
    plot(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], layout=l,
        title=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)"],
        titleloc=:left, titlefont=16, size=(800, 750), left_margin=3mm, right_margin=0mm,
        bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotCa.pdf")
end

function plotLCCti()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度

    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []
    Nmax = 10
    for N in ["LCCIt05i05", "LCCIt150i05", "LCCIt10005i05", "LCCIt67i01", "LCCIt67i15", "LCCIt67i151"]

        if N == "LCCIt67i151" || N == "LCCIt10005i05"
            colorbar = :legend
        else
            colorbar = :none
        end
        if N == "LCCIt05i05" || N == "LCCIt67i01"
            ylabel = L"N^{\rm peak}_{\rm O}"
        else
            ylabel = ""
        end

        if N == "LCCIt67i01" || N == "LCCIt67i15" || N == "LCCIt67i151"
            xlabel = L"\tau\ (\rm ms)"
        else
            xlabel = ""
        end



        name = N
        if N == "LCCIt05i05"
            name1 = L"t^\textrm{LCC}=0.5\ \textrm{ms}"
            name2 = L"i^\textrm{LCC}=0.5\ \textrm{pA}"
        elseif N == "LCCIt150i05"
            name1 = L"t^\textrm{LCC}=15\ \textrm{ms}"
            name2 = L"i^\textrm{LCC}=0.5\ \textrm{pA}"
        elseif N == "LCCIt10005i05"
            name1 = L"t^\textrm{LCC}:\textrm{O/C\ per}\ 0.5\ \textrm{ms}"
            name2 = L"i^\textrm{LCC}=0.5\ \textrm{pA}"
        elseif N == "LCCIt67i01"
            name1 = L"t^\textrm{LCC}=6.7\ \textrm{ms}"
            name2 = L"i^\textrm{LCC}=0.1\ \textrm{pA}"
        elseif N == "LCCIt67i15"
            name1 = L"t^\textrm{LCC}=6.7\ \textrm{ms}"
            name2 = L"i^\textrm{LCC}=15\ \textrm{pA}"
        elseif N == "LCCIt67i151"
            name1 = L"t^\textrm{LCC}=6.7\ \textrm{ms}"
            name2 = L"i^\textrm{LCC}=15(1-t/6.7)\ \textrm{pA}"
        end


        Npeak, tCRE, TTS = CRE(name)
        println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
        println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
        println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))
        tlim = 60
        Nlim = 60
        t0 = histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
            xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.07, tlim * 0.07],
                [60 * 0.8775, 50 * 0.8775], [Plots.text(name1, 10, :left), Plots.text(name2, 10, :left)]), xlabel=xlabel,
            guidefont=(16, :dark), framestyle=:box, grid=:none, ylabel=ylabel, tickfontsize=11,
            dpi=dpi, colorbar_scale=:log10, colorbar=colorbar)

        push!(t, t0)

    end

    l = @layout[[a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]]
    plot(t[1], t[2], t[3], t[4], t[5], t[6], layout=l,
        title=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)"],
        titleloc=:left, titlefont=12, size=(800, 450), left_margin=3mm, right_margin=0mm,
        bottom_margin=4mm, top_margin=0mm, link=:all, dpi=dpi)



    Plots.savefig("results/plotLCCti.pdf")
end

function plotCSQ0()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度


    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []
    Nmax = 50 / 4.0
    for C in ["0", "200", "400", "600", "800", "1000", "1200", "1400", "1600"]

        if C in ["1200", "1400", "1600"]
            xlabel = L"\tau\ (\rm ms)"
        else
            xlabel = ""
        end
        if C in ["0", "600", "1200"]
            ylabel = L"N^{\rm peak}_{\rm O}"
        else
            ylabel = ""
        end

        if C in ["400", "1000", "1600"]
            colorbar = :legend
        else
            colorbar = :none
        end



        name = "CSQ" * C
        Cint = Int64(parse(Int64, C))
        nameN = L"\textrm{CSQ}: %$Cint"


        Npeak, tCRE, TTS = CRE(name)
        Npeaknb, tCREnb, TTSnb = CRE("CSQnb" * C)


        println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
        println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
        println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))

        tlim = 70
        Nlim = 60
        t0 = histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
            xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.08],
                [Nlim * 0.8375], [Plots.text(nameN, 12, :left)]), xlabel=xlabel,
            guidefont=(16, :dark), framestyle=:box, grid=:none, ylabel=ylabel, tickfontsize=11,
            dpi=dpi, colorbar_scale=:log10, colorbar=colorbar)

        scatter!([mean(tCRE[Npeak.>Nmax])], [mean(Npeak[Npeak.>Nmax])], xerror=[std(tCRE[Npeak.>Nmax])],
            yerror=[std(Npeak[Npeak.>Nmax])], markershape=:circle, markersize=4, markercolor=:black,
            markerstrokewidth=0.5, markerstrokecolor=:black, label=" Full model", legend_position=:bottomright,
            legend_font_pointsize=8)
        scatter!([mean(tCREnb[Npeaknb.>Nmax])], [mean(Npeaknb[Npeaknb.>Nmax])], xerror=[std(tCREnb[Npeaknb.>Nmax])],
            yerror=[std(Npeaknb[Npeaknb.>Nmax])], markershape=:diamond, markersize=4, markercolor=:blue,
            markerstrokewidth=0.5, markerstrokecolor=:blue, label=" No CSQ buffer ", legend_position=:bottomright,
            legend_font_pointsize=8)

        push!(t, t0)



    end



    l = @layout[[a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]; [a{0.29w} b{0.29w} c{0.42w}]]

    plot(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], layout=l,
        title=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)"],
        titleloc=:left, titlefont=16, size=(800, 750), left_margin=3mm, right_margin=0mm,
        bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotCSQ.pdf")
end

function plotCSQ()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度


    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []
    Nmax = 50 / 4.0
    plotsc = scatter(xlims=(-100, 1700), xticks=0:400:1700, ylims=(0, 70), markershape=:circle,
        markercolor=:lightpink, markerstrokewidth=0.5, markerstrokecolor=:lightpink, size=(600, 500),
        legend_position=:topleft, legend_font_pointsize=10, guidefont=12, framestyle=:box, grid=:none,
        tickfontsize=11, xlabel=L"B_\textrm{CSQ} (\mu\textrm{M})", ylabel=L"\tau\ (\textrm{ms})")

    for C in ["0", "50", "100", "200", "400", "600", "800", "1000", "1200", "1400", "1600"]

        name = "CSQ" * C
        Cint = Int64(parse(Int64, C))
        nameN = L"\textrm{CSQ}: %$Cint"
        if C == "0"
            label1 = " Full model"
            label2 = " No buffer function "
            label3 = " No binding states "
        else
            label1 = ""
            label2 = ""
            label3 = ""
        end

        Npeak, tCRE, TTS = CRE(name)
        Npeaknb, tCREnb, TTSnb = CRE("CSQnb" * C)
        Npeaknbb, tCREnbb, TTSnbb = CRE("CSQ" * C * "nbb")


        println(name * "peak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
        println(name * "tCRE:  ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
        println(name * "TTS:  ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))

        scatter!(plotsc, [Cint], [mean(tCRE[Npeak.>Nmax])], yerror=[std(tCRE[Npeak.>Nmax])],
            markershape=:circle, markersize=6, markercolor=:black, markerstrokewidth=0.5,
            markerstrokecolor=:black, label=label1)

        scatter!(plotsc, [Cint], [mean(tCREnb[Npeaknb.>Nmax])], yerror=[std(tCREnb[Npeaknb.>Nmax])],
            markershape=:circle, markersize=6, markercolor=:white, markerstrokewidth=0.5,
            markerstrokecolor=:green, label=label2)

        scatter!(plotsc, [Cint], [mean(tCREnbb[Npeaknbb.>Nmax])], yerror=[std(tCREnbb[Npeaknbb.>Nmax])],
            markershape=:circle, markersize=6, markercolor=:white, markerstrokewidth=0.5,
            markerstrokecolor=:black, label=label3)
    end

    Plots.savefig("results/plotCSQ.pdf")
end

function plotCSQN()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    dpi = 1200 #画图精度
    guifont = (16, :dark)
    tickfontsize = 11
    Nt, N, t = [], [], []
    cg = cgrad([:darkgreen, :red, :yellow], [0.5, 0.7, 0.9])
    cg = cgrad([RGB(0.9, 0.9, 0.9), :darkgreen, :red, :yellow], [0.0, 0.3, 0.5, 0.9])

    t, Nt, Nta = [], [], []


    plotsc = scatter(xlims=(0, 320), ylims=(0, 1.2), markershape=:circle,
        markercolor=:lightpink, markerstrokewidth=0.5, markerstrokecolor=:lightpink,
        legend_position=:bottomright, legend_font_pointsize=10, guidefont=guifont, framestyle=:box, grid=:none,
        tickfontsize=tickfontsize, xlabel=L"\Delta t\ (\rm ms)", ylabel=L"\textrm{Fraction\ of\ first\ spark}")
    paras = "LCCE"
    paras_plot = "results/" * paras * "1_plot.jld" #数据文件名字
    tJSR, Ss = load(paras_plot, "t")[2], load(paras_plot, "Ss")[2]
    x0, y0 = tJSR, [sum(S) / length(S)[1] for S in Ss] #时间，平均值
    plot!(plotsc, (x0 .- 6.0) .* 1000.0, y0 ./ 1000.0, linewidth=2.0, color=:violet,
        label=L"[\textrm{\rm Ca}^{2+}]^\textrm{ JSR}/1000\mu\textrm{M}") #画图
    plot!(plotsc, [0, 300], [1, 1], linewidth=2.0, guidefont=guifont, framestyle=:box,
        grid=:none, tickfontsize=tickfontsize, color=:grey, label="")

    f = open("results/expBdata.txt", "r")
    lines = readlines(f)
    close(f)
    datat, datav, datae = Float64[], Float64[], Float64[]
    for line in lines
        datat0, datav1, datav2 = parse.(Float64, split(line, r"[,\s]+"))
        datav0 = (datav2 + datav1) / 200
        datae0 = (datav2 - datav1) / 200
        push!(datat, datat0)
        push!(datav, datav0)
        push!(datae, datae0)
    end
    #@show datat,datav,datae
    scatter!(plotsc, datat, datav, yerror=datae,
        markershape=:star, markercolor=:grey, markerstrokewidth=0.5, markersize=10,
        markerstrokecolor=:grey, label="Brochet et al., Data")

    f = open("results/expBl.txt", "r")
    lines = readlines(f)
    close(f)
    datat, datav = Float64[], Float64[], Float64[]
    for line in lines
        datat0, datav0 = parse.(Float64, split(line, r"[,\s]+"))
        datav0 = datav0 / 100
        push!(datat, datat0)
        push!(datav, datav0)
    end
    #@show datat,datav,datae
    plot!(plotsc, datat, datav, linewidth=2.0, label="Brochet et al.,fit", color=:grey)



    for C in ["100", "400", "1200"]
        ltt = 0
        for itt in 25:25:300
            tt = string(itt)

            name = "CSQN" * C * "t" * tt
            namenb = "CSQN" * C * "t" * tt * "nb"
            namenbb = "CSQN" * C * "t" * tt * "nbb"
            Cint = Int64(parse(Int64, C))
            tint = Int64(parse(Int64, tt))
            nameN = L"\textrm{CSQ}: %$Cint,~\Delta t: %$tint"

            Npeak, tCRE, TTS, CaSS = CRE(name, No=0, tp=tint / 1000.0, CaSS=false)
            Npeak1, tCRE1, TTS1, CaSS1 = CRE(name, No=1, tp=tint / 1000.0, CaSS=false)
            Npeak2, tCRE2, TTS2, CaSS2 = CRE(name, No=2, tp=tint / 1000.0, CaSS=false)
            if C == "400"
                Npeaknb1, tCREnb1, TTSnb1, CaSSnb1 = CRE(namenb, No=1, tp=tint / 1000.0, CaSS=false)
                Npeaknb2, tCREnb2, TTSnb2, CaSSnb2 = CRE(namenb, No=2, tp=tint / 1000.0, CaSS=false)
                Npeaknbb1, tCREnbb1, TTSnbb1, CaSSnbb1 = CRE(namenbb, No=1, tp=tint / 1000.0, CaSS=false)
                Npeaknbb2, tCREnbb2, TTSnbb2, CaSSnbb2 = CRE(namenbb, No=2, tp=tint / 1000.0, CaSS=false)
            end
            Nmax = 5
            if C == "1200" && itt == 25
                Nmax = 0
            end

            lengthN1 = length(Npeak1[Npeak1.>Nmax])
            lengthN2 = length(Npeak2[Npeak2.>Nmax])
            meanN1, stdN1 = mean(Npeak1[Npeak1.>Nmax]), std(Npeak1[Npeak1.>Nmax])
            meanN2, stdN2 = mean(Npeak2[Npeak2.>Nmax]), std(Npeak2[Npeak2.>Nmax])
            meant1, stdt1 = mean(tCRE1[Npeak1.>Nmax]), std(tCRE1[Npeak1.>Nmax])
            meant2, stdt2 = mean(tCRE2[Npeak2.>Nmax]), std(tCRE2[Npeak2.>Nmax])
            if C == "400"
                meanNnb1, stdNnb1 = mean(Npeaknb1[Npeaknb1.>Nmax]), std(Npeaknb1[Npeaknb1.>Nmax])
                meanNnb2, stdNnb2 = mean(Npeaknb2[Npeaknb2.>Nmax]), std(Npeaknb2[Npeaknb2.>Nmax])
                meanNnbb1, stdNnbb1 = mean(Npeaknbb1[Npeaknbb1.>Nmax]), std(Npeaknbb1[Npeaknbb1.>Nmax])
                meanNnbb2, stdNnbb2 = mean(Npeaknbb2[Npeaknbb2.>Nmax]), std(Npeaknbb2[Npeaknbb2.>Nmax])
            end
            @show nameN, lengthN1, lengthN2

            if C == "100"
                color = :blue
                markershape = :hexagon

            elseif C == "400"
                color = :black
                markershape = :circle
            elseif C == "1200"
                color = :red
                markershape = :utriangle
            end

            if ltt == 0
                label = "CSQ:" * C
                labelnb = "CSQ:" * C * ", no binding states "
                labelnbb = "CSQ:" * C * ", no buffer function "
                ltt = 1
            else
                label = ""
                labelnb = ""
                labelnbb = ""
            end

            if tt in ["50", "100", "200"]

                if C == "1200"
                    xlabel = L"\tau\ (\rm ms)"
                else
                    xlabel = ""
                end
                if tt == "50"
                    ylabel = L"N^{\rm peak}_{\rm O}"
                else
                    ylabel = ""
                end

                if tt == "200"
                    colorbar = :legend
                else
                    colorbar = :none
                end




                tlim = 70
                Nlim = 60
                t0 = histogram2d(tCRE, Npeak, bins=(range(0, tlim, step=1), range(0, Nlim, step=1)),
                    xlims=(0, tlim), ylims=(0, Nlim), clims=(1, 100), color=cg, annotations=([tlim * 0.07],
                        [Nlim * 0.87], [Plots.text(nameN, 12, :left)]), xlabel=xlabel,
                    guidefont=guifont, framestyle=:box, grid=:none, ylabel=ylabel, legend_font_pointsize=9,
                    dpi=dpi, colorbar_scale=:log10, colorbar=colorbar, tickfontsize=tickfontsize)

                scatter!([meant1], [meanN1], xerror=[stdt1], yerror=[stdN1],
                    markershape=:circle, markercolor=:black, markerstrokewidth=0.5, markersize=4,
                    markerstrokecolor=:black, legend_position=:bottomright, label=" 1st spark ")
                scatter!([meant2], [meanN2], xerror=[stdt2], yerror=[stdN2],
                    markershape=:diamond, markercolor=:blue, markerstrokewidth=0.5, markersize=4,
                    markerstrokecolor=:blue, legend_position=:bottomright, label=" 2nd spark ")

                push!(t, t0)

            end
            if C == "400"
                scatter!(plotsc, [tint], [meanNnb2 / meanNnb1], yerror=[stdNnb2 / meanNnb1],
                    markershape=:circle, markersize=6, markercolor=:white, markerstrokewidth=0.5,
                    markerstrokecolor=:black, label=labelnb)
                scatter!(plotsc, [tint], [meanNnbb2 / meanNnbb1], yerror=[stdNnbb2 / meanNnbb1],
                    markershape=:circle, markersize=6, markercolor=:white, markerstrokewidth=0.5,
                    markerstrokecolor=:green, label=labelnbb)
            end

            scatter!(plotsc, [tint], [meanN2 / meanN1], yerror=[stdN2 / meanN1],
                markershape=markershape, markersize=6, markercolor=color, markerstrokewidth=0.5,
                markerstrokecolor=color, label=label)


        end
    end



    l = @layout[[[a{0.26w} b{0.26w} c{0.48w}]; [a{0.26w} b{0.26w} c{0.48w}]; [a{0.26w} b{0.26w} c{0.48w}]] e{0.35w}]

    plot(t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], plotsc, layout=l,
        title=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)"],
        titleloc=:left, titlefont=16, size=(1500, 800), left_margin=6mm, right_margin=0mm,
        bottom_margin=2.9mm, top_margin=0.5mm, link=:all, dpi=dpi)

    Plots.savefig("results/plotCSQN.pdf")
end

#----------------------------------------------------------------------------------------------------------
function data(x, y, n)
    pointss = []
    pointss = [[x[i], y[i], n[i]] for i in eachindex(x)]

    i = 1
    points_s = Vector{Vector{Float64}}[]
    temp = Vector{Float64}[]
    push!(temp, pointss[1])
    lpoints = length(pointss)
    while i < lpoints
        if pointss[i][3] == pointss[i+1][3]
            push!(temp, pointss[i+1])
            i += 1
        else
            push!(points_s, temp)
            temp = Vector{Float64}[]
            i += 1
            push!(temp, pointss[i])
        end
    end
    push!(points_s, temp)
    points_s
end
function sd(points_s, N)

    points_sd = Vector{Vector{Float64}}[] #整理好且筛选长度后的数组
    for points0 in points_s
        if length(points0) < N
            continue
        elseif length(points0) >= N
            push!(points_sd, points0)
        end
    end
    points_sd

end
function d2d4(x, y, n)

    points_s = data(x, y, n)

    points_sd = sd(points_s, 2)

    d2 = Float64[]
    for points in points_sd, n in range(1, length(points))
        min_dist = 999999
        for m in range(1, length(points))
            if points[n] != points[m] && points[m][3] == points[n][3]
                temp_dist = sqrt((points[m][1] - points[n][1])^2 + (points[m][2] - points[n][2])^2)
                if temp_dist < min_dist
                    min_dist = temp_dist
                end
            end
        end
        push!(d2, min_dist)
    end


    points_sd = sd(points_s, 5)

    #距离归类
    tempt = Float64[]#团簇中第i个点和其他所有点的距离[dis,...]
    temp_1 = Vector{Float64}[] #归纳一个团簇中每个点和其他所有点的距离[[dis1,...],[dis2,...]]
    temp_s = Vector{Vector{Float64}}[] #[[[dis1,..,],[dis2,..]]第一个团簇,[[dis1,...],[dis2,...]]第二个团簇]
    k = 1
    while k <= length(points_sd)
        for n in range(1, length(points_sd[k]))
            for m in range(1, length(points_sd[k]))
                if points_sd[k][n] != points_sd[k][m]
                    temp_dist = sqrt((points_sd[k][m][1] - points_sd[k][n][1])^2 + (points_sd[k][m][2] - points_sd[k][n][2])^2)
                    push!(tempt, temp_dist)
                end
            end
            push!(temp_1, tempt)
            tempt = Float64[]
        end
        push!(temp_s, temp_1)
        temp_1 = Vector{Float64}[]
        k += 1
    end

    #选出最小的前四个求距离平均值
    d4 = Float64[]
    for t in temp_s
        for r in t
            d = sort(r)[1:4] #排序：每个团簇中第i个点到其他所有点的距离，然后只取前4个，求平均
            d_average = sum(d) / length(d)
            push!(d4, d_average)
        end
    end
    return d2, d4
end
function plot_net()
    ENV["GKSwstype"] = "100" #不弹出图像窗口



    paras = "SponStest" #读取需要画图的数据文件名字
    net = load("results/" * paras * "_net.jld", "net") #读取cluster网络的数据
    N, size, x, y, A = net.N, net.sizeV, net.x, net.y, net.A
    @show paras, N
    x = x ./ 1000.0
    y = y ./ 1000.0 #将坐标缩小1000倍
    NC = length(size) #number of cluster
    dpi = 1200 #画图精度
    markersize = 2.0
    linewidth = 0.00001
    guifont = (17, :dark)
    tickfontsize = 12 #点大小，线宽度，框架字体，刻度字体


    #设置画图格式，空图
    pnet = plot(guidefont=guifont, framestyle=:box, xlabel=L"x~(\mu \textrm{m})",
        ylabel=L"y~(\mu \textrm{m})", dpi=dpi, xlims=(0.0, maximum(x)), ylims=(0.0, maximum(y)),
        xticks=0:0.5:4, yticks=0:0.5:4,
        tickfontsize=tickfontsize, legend=:none, size=(300, 300), grid=:none)
    #画连线 
    Amax = maximum(A)
    lcolor = (Amax .- A) / Amax

    Aij = [A[i, j] for i in 1:N for j in i:N]
    iv = [i for i in 1:N for j in i:N]
    jv = [j for i in 1:N for j in i:N]
    sort = sortperm(Aij)

    for s in ProgressBar(sort)
        i = iv[s]
        j = jv[s]
        if A[i, j] > 1000
            plot!(pnet, [x[i], x[j]], [y[i], y[j]], legend=:none, linecolor=RGB(lcolor[i, j], lcolor[i, j], lcolor[i, j]), linewidth=A[i, j] * linewidth,) #链接
        end
    end


    #画点，不同团簇颜色不同
    NN = 1
    for i in 1:NC
        scatter!(pnet, x[NN:NN+size[i]-1], y[NN:NN+size[i]-1], markercolor=i, markersize=markersize,
            markerstrokewidth=0) #节点
        NN = NN + size[i]
    end




    x, y, n = net.x, net.y, net.n

    d2, d4 = d2d4(x, y, n)

    ENV["GKSwstype"] = "100" #不弹出图像窗口
    f_size(x) = (0.991 * exp(-0.66 * x) + 0.009 * exp(-0.017 * x))#/1.84*500
    println(size)

    histogram(size, bins=100, normalize=:pdf, grid=false, guidefont=13, xlabel="Cluster size",
        ylabel="Clusters", dpi=dpi, xlims=(0, 100), ylims=(0.00, 1), framestyle=:box,
        legend=:bottomright, label=:none, tickfontsize=tickfontsize)
    pc = plot!(f_size, linecolor=:red, label=:none, guidefont=13, tickfontsize=tickfontsize)


    histogram(d2, color="red", xlabel="distance (nm)", ylabel=" RyRs",
        xlims=(0.0, 100), ylims=(0.0, 0.4), bins=range(0, 100, length=20), normalize=:probability,
        leg=false, grid=:none, framestyle=:box, guidefont=13, tickfontsize=tickfontsize)

    pd = histogram!(d4, color="green", xlabel="distance (nm)", ylabel="RyRs", bins=range(0, 100, length=20),
        normalize=:probability, alpha=0.5, grid=:none, guidefont=13, framestyle=:box, tickfontsize=tickfontsize)

    xc, yc, nc = [], [], []

    ii = 1
    for i in net.sizeV
        push!(xc, sum(x[ii:ii+i-1]) / i)
        push!(yc, sum(y[ii:ii+i-1]) / i)
        push!(nc, 1)
        ii += i
    end

    d2, d4 = d2d4(xc, yc, nc)

    ENV["GKSwstype"] = "100" #不弹出图像窗口
    histogram(d2, color="red", xlabel="distance (nm)", ylabel="Clusters",
        xlims=(0.0, 500), ylims=(0.0, 0.4), bins=range(0, 500, length=20),
        normalize=:probability, leg=false, grid=:none, guidefont=13, framestyle=:box, tickfontsize=tickfontsize)

    pdc = histogram!(d4, color="green", xlabel="distance (nm)", ylabel="Cluesters",
        normalize=:probability, bins=range(0, 500, length=20), alpha=0.5, grid=:none, guidefont=13, framestyle=:box, tickfontsize=tickfontsize)

    l = @layout[a{0.7w} [e{0.33h}; f{0.33h}; f{0.33h}]]

    plot(pnet, pc, pd, pdc, layout=l, title=["(a)" "(b)" "(c)" "(d)"], titleloc=:left, titlefont=14, size=(800, 500), left_margin=0mm, right_margin=2mm, bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)


    Plots.savefig("results/net.pdf")



end

function plottest()
    ENV["GKSwstype"] = "100" #不弹出图像窗口
    #paras="Singletest"
    paras = load("results/name.jld", "name") #读取需要画图的数据文件名字
    dpi = 1200 #画图精度
    guifont = (18, :dark)
    tickfontsize = 12
    #演化图
    paras_plot = "results/" * paras * "1_plot.jld" #数据文件名字
    #读取数据
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[1] #时间，Ss中第一个数值
    @show length(t), length(Ss[1])[1]
    Nnode = length(Ss[1])[1]
    cg = cgrad([:white, :green, :black], [0, 0.5, 0.75, 1.0])

    #节点演化图，
    x0, y0, z0 = t, [i for i in 1:Nnode], [S[i] for i in 1:Nnode, S in Ss] #时间，节点序号，浓度
    cmin, cmax = minimum(z0), maximum(z0) + 1e-20 #浓度范围
    p1 = heatmap(x0, y0, z0, c=cg, xlims=(0.0, maximum(x0)), ylims=(0.0, maximum(y0)), clims=(cmin, cmax), framestyle=:box, guidefont=guifont, tickfontsize=tickfontsize, dpi=dpi, ylabel=L"\textrm{Node}", colorbar_title=L"[\textrm{\rm Ca}^{2+}]_\textrm{SS}", colorbar_titlefontrotation=0, colorbar_titlefontsize=14, colorbar_tickfontsize=8) #画图

    #读取数据
    t, Ss = load(paras_plot, "t")[2], load(paras_plot, "Ss")[2]


    #节点演化图
    x0, y0, z0 = t, [i for i in 1:Nnode], [S[i] for i in 1:Nnode, S in Ss]
    cmin, cmax = minimum(z0), maximum(z0) + 1e-5

    p2 = heatmap(x0, y0, z0, c=cg, xlims=(0.0, maximum(x0)), ylims=(0.0, maximum(y0)), clims=(cmin, cmax),
        framestyle=:box, guidefont=guifont, tickfontsize=tickfontsize,
        dpi=dpi, ylabel=L"\textrm{Node}", colorbar_title=L"[\textrm{\rm Ca}^{2+}]_\textrm{JSR}", colorbar_titlefontrotation=0, colorbar_titlefontsize=14, colorbar_tickfontsize=8)

    #读取数据
    t, Ss = load(paras_plot, "t")[1], load(paras_plot, "Ss")[3]

    #节点演化图
    x0, y0, z0 = t, [i for i in 1:Nnode], [S[i] for i in 1:Nnode, S in Ss]
    cmin, cmax = minimum(z0), maximum(z0) + 1e-5

    p3 = heatmap(x0, y0, z0, c=cg, xlims=(0.0, maximum(x0)), ylims=(0.0, maximum(y0)), clims=(cmin, cmax),
        framestyle=:box, guidefont=guifont, tickfontsize=tickfontsize,
        dpi=dpi, xlabel=L"t~(s)", ylabel=L"\textrm{Node}", colorbar_title=L"N_\textrm{O}", colorbar_titlefontrotation=0, colorbar_titlefontsize=14, colorbar_tickfontsize=8)



    Npeak, tCRE, TTS = CRE(paras)



    Nmax = 15
    println("Npeak:  ", mean(Npeak[Npeak.>Nmax]), "+-", std(Npeak[Npeak.>Nmax]), "    ", length(Npeak[Npeak.>Nmax]))
    println("tCRE: ", mean(tCRE[Npeak.>Nmax]), "+-", std(tCRE[Npeak.>Nmax]), "    ", length(tCRE[Npeak.>Nmax]))
    println("TTS: ", mean(TTS[Npeak.>Nmax]), "+-", std(TTS[Npeak.>Nmax]), "    ", length(TTS[Npeak.>Nmax]))


    Nt = histogram2d(tCRE, Npeak, bins=(range(0, 80, step=0.1), range(0, 60, step=1)), xlims=(0, 80), ylims=(0, 60), color=cgrad([:brown, :blue, :black], [0.3, 0.6, 0.99]), guidefont=guifont, framestyle=:box, grid=:none, ylabel=L"N^{\rm peak}_{\rm O}", xlabel=L"\tau\ (\rm ms)", dpi=dpi, colorbar_scale=:log10, colorbar_titlefontrotation=180, clims=(1, 100), colorbar_tickfontsize=6, tickfontsize=tickfontsize)


    N = histogram(Npeak, yaxis=:log, color="red", xlabel=L"N^{\rm peak}_{\rm O}", ylabel=L"\rm Number\  of\  CRE", bins=range(0, 80, step=1), leg=false, grid=:none, framestyle=:box, guidefont=guifont, xlims=(0, 80), ylims=(0.1, 1e5), tickfontsize=tickfontsize)


    t = histogram(tCRE, yaxis=:log, color="red", xlabel=L"\tau\ (\rm ms)", bins=range(0, 80, step=1), leg=false, grid=:none, framestyle=:box, guidefont=guifont, xlims=(0, 80), ylims=(0.1, 1e5), tickfontsize=tickfontsize)



    l = @layout[a{0.97w}; b; c; grid(1, 3, widths=(0.29, 0.29, 0.42)){0.32h}]

    plot(p1, p2, p3, N, t, Nt, layout=l, title=["(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)"], titleloc=:left, titlefont=12, size=(1000, 900), left_margin=0mm, right_margin=1mm, bottom_margin=0mm, top_margin=0mm, link=:all, dpi=dpi)


    Plots.savefig("results/test.pdf")
end


#----------------------------------------------------------------------------------------------------------
#@time plotSponE()
#@time plotLCCE()
#@time plotSingle()
#@time plotCa()
#@time plotCa0()
#@time plotrateline()
@time plotNSRmyo()
#@time plotLCCti()
#@time plotCSQ()
#@time plotCSQN()
#@time plot_net()
#@time plottest()
