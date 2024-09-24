using JuMP, Graphs, GraphPlot, Plots #最优化，图论，图论画图，画图
using HiGHS, Ipopt, QuadGK #线性最优化，非线性最优化，高斯积分
using Distributions, Colors, Cairo, Compose, JLD2, LaTeXStrings
using GLPK
using Juniper
import Base: redirect_stdout, redirect_stderr
function sigmoid(x, p)
    return p[1] * x^p[3] / (x^p[3] + p[2]^p[3])
end
function CC(x, p)
    return p[1]
end

#利用JuMP拟合
function fit(d, f)
    # 构建 JuMP 优化模型，使用 Juniper 作为求解器
    optimizer = Juniper.Optimizer
    model = Model(optimizer)

    # 设置连续求解器为 Ipopt，整数求解器为 GLPK
    set_optimizer_attribute(model, "nl_solver", Ipopt.Optimizer)
    set_optimizer_attribute(model, "mip_solver", GLPK.Optimizer)



    # 定义 ko, k 为实数，p3 为整数
    @variable(model, p1 >= 0)     # p1 对应 ko，非负实数
    @variable(model, p2 >= 0)     # p2 对应 k，非负实数
    @variable(model, p3 >= 0, Int) # p3 为非负整数

    # 将这些变量作为参数数组
    p = [p1, p2, p3]

    # 目标函数：最小化模型预测值与观测值之间的加权误差
    @objective(model, Min, sum((
        ((f(d[1][i], p) - d[2][i]) / ((d[3][i] - d[4][i]) / 2.0))^2
    ) for i in 1:length(d[1])))

    # 进行优化
    optimize!(model)
    # 获取拟合结果

    return value.(p)
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

function main()
    ENV["GKSwstype"] = "100" #不弹出图像窗口

    d = readdata("closedata2.txt")
    p = fit(d, CC)
    p_closerat = deepcopy(p)
    scatter(d[1], d[2], yerror=(d[2] - d[4], d[3] - d[2]),
        grid=:none, framestyle=:box, capsize=0, markerstrokewidth=0.5,
        markerstrokecolor=:red, color=:red, linecolor=:red,
        label="Rat: Data")

    d = readdata("closeline2.txt")
    d[1] = collect(1:1000)
    d[2] = CC.(d[1], Ref(p))
    plot!(d[1], d[2], color=:red, linecolor=:red, linestyle=:solid, label="Rat: Sigmiod")

    d[2] = min.(max.(7906 * d[1] .^ -0.5, 900.0), 200000.0)
    plot!(d[1], d[2], color=:red, linecolor=:red, linestyle=:dot, label="Rat: Cannell et al.")

    d = readdata("closedata1.txt")
    p = fit(d, CC)
    p_closesheep = deepcopy(p)
    # 绘制散点图，带不对称的 Y 轴误差条
    scatter!(d[1], d[2], yerror=(d[2] - d[4], d[3] - d[2]), label="Sheep: Data",
        xlabel=L"[\textrm{Ca}^{2+}]^\textrm{SS}~(\mu \textrm{M})", ylabel=L"k_\textrm{C}~(\textrm{s}^{-1})",marker=:diamond,
        markerstrokecolor=:blue, color=:blue, linecolor=:blue,
        xscale=:log10, yscale=:log10, xlim=(1, 1000), ylim=(100, 30000), legend=:topright)

    d = readdata("closeline1.txt")
    d[1] = collect(1:1000)
    d[2] = CC.(d[1], Ref(p))
    plot!(d[1], d[2], label="Sheep: Sigmiod", linecolor=:blue, linestyle=:solid)
    d[2] .= 1500
    plot!(d[1], d[2], label="Sheep: 1500/s", linecolor=:blue, linestyle=:dash)
    d[2] = 1581.85 .* d[1] .^ -0.27
    figclose = plot!(d[1], d[2], label="Sheep: Cannell et al.", linecolor=:blue, linestyle=:dot)


    ##########################

    d = readdata("opendata2.txt")
    p = fit(d, sigmoid)
    p_openrat = deepcopy(p)
    scatter(d[1], d[2], yerror=(d[2] - d[4], d[3] - d[2]),
        markercolor=:red, markerstrokewidth=0.5, markerstrokecolor=:red, color=:red, linecolor=:red,
        label="Rat: Data")

    d = readdata("openline2.txt")
    d[1] = collect(1:1000)
    d[2] = sigmoid.(d[1], Ref(p))
    plot!(d[1], d[2], color=:red, linecolor=:red, linestyle=:solid, label="Rat: Sigmiod")
    p[2] = 0.5 * p[2]
    d[2] = sigmoid.(d[1], Ref(p))
    plot!(d[1], d[2], label="Rat: 0.5 K", linecolor=:red, linestyle=:dash)
    d[2] = max.(min.(1.262e-3 .* d[1] .^ 2.8, 700.0), 0.0)
    plot!(d[1], d[2], color=:red, linecolor=:red, linestyle=:dot, label="Rat: Cannell et al.")

    d = readdata("opendata1.txt")
    p = fit(d, sigmoid)
    p_opensheep = deepcopy(p)
    # 绘制散点图，带不对称的 Y 轴误差条
    scatter!(d[1], d[2], yerror=(d[2] - d[4], d[3] - d[2]), label="Sheep: Data",
        xlabel=L"[\textrm{Ca}^{2+}]^\textrm{SS}~(\mu \textrm{M})", ylabel=L"k_\textrm{O}~(\textrm{s}^{-1})",marker=:diamond,
        grid=:none, framestyle=:box, color=:blue, linecolor=:blue, capsize=0,
        markercolor=:blue, markerstrokewidth=0.5, markerstrokecolor=:blue,
        xscale=:log10, yscale=:log10, xlim=(1, 1000), ylim=(0.01, 2000), legend=:bottomright)

    d = readdata("openline1.txt")
    d[1] = collect(1:1000)
    d[2] = sigmoid.(d[1], Ref(p))
    plot!(d[1], d[2], label="Sheep: Sigmiod", linecolor=:blue, linestyle=:solid)
    d[2] = max.(min.(0.199488 .* d[1] .^ 2.12, 800.0), 0.0)
    figopen = plot!(d[1], d[2], label="Sheep: Cannell et al.", linecolor=:blue, linestyle=:dot)


    @show p_openrat
    @show p_opensheep
    @show p_closerat
    @show p_closesheep

    dx = collect(1:1000)
    p_openrat[2] = p_openrat[2]
    dO = sigmoid.(dx, Ref(p_openrat))
    dC = CC.(d[1], Ref(p_closerat))
    rat = dO ./ (dO .+ dC)
    plot(dx, rat, grid=:none, framestyle=:box, xscale=:log10,
        linecolor=:red, linestyle=:solid, label="Rat: Sigmiod",
        ylabel=L"P_\textrm{O}", xlabel=L"[\textrm{Ca}^{2+}]^\textrm{SS}~(\mu \textrm{M})",
        xlim=(1, 1000), ylim=(0.0, 1.2))
    p_openrat[2] = p_openrat[2] / 2.0
    dO = sigmoid.(dx, Ref(p_openrat))
    dC = CC.(d[1], Ref(p_closerat))
    rat = dO ./ (dO .+ dC)
    plot!(dx, rat, linecolor=:red, linestyle=:dash, label="Rat: 0.5 K",)
    dO = max.(min.(1.262e-3 .* d[1] .^ 2.8, 700.0), 0.0)
    dC = min.(max.(7906 * d[1] .^ -0.5, 900.0), 200000.0)
    rat = dO ./ (dO .+ dC)
    plot!(dx, rat, grid=:none, framestyle=:box, xscale=:log10,
        linecolor=:red, linestyle=:dot, label="rate: Cannell et al.",
        xlim=(1, 1000), ylim=(0.0, 1.0))
    dO = sigmoid.(dx, Ref(p_opensheep))
    dC = CC.(d[1], Ref(p_closesheep))
    sheep = dO ./ (dO .+ dC)
    plot!(dx, sheep, grid=:none, framestyle=:box, xscale=:log10,
        linecolor=:blue, linestyle=:solid, label="Sheep: Sigmiod",
        xlim=(1, 1000), ylim=(0.0, 1.0))
    dO = sigmoid.(dx, Ref(p_opensheep))
    dC .= 1500.
    sheep = dO ./ (dO .+ dC)
    plot!(dx, sheep, grid=:none, framestyle=:box, xscale=:log10,
        linecolor=:blue, linestyle=:dash, label="Sheep: 1500/s",
        xlim=(1, 1000), ylim=(0.0, 1.0))

    dO = max.(min.(0.199488 .* dx .^ 2.12, 800.0), 0.0)
    dC = 1581.85 .* dx .^ -0.27
    sheep = dO ./ (dO .+ dC)
    figPo = plot!(dx, sheep, grid=:none, framestyle=:box, xscale=:log10,
        linecolor=:blue, linestyle=:dot, label="Sheep: Cannell et al.",
        xlim=(1, 1000), ylim=(0.0, 1.0))



    l = @layout [a b c]
    Plots.plot(figopen, figclose, figPo, layout=l, titleloc=:left, titlefont=14, size=(1200, 350),title=["(a)" "(b)" "(c)"],
        left_margin=6mm, right_margin=1mm, bottom_margin=8mm, top_margin=1mm, link=:all)




    # 保存图像为指定的文件
    savefig("rate.pdf")


end

main()
