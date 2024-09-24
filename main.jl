#--------------------------11--------------------------------------------------------------------------------
#主文件
#包含main函数和work函数。main函数包含参数与主要模拟流程。work函数包含模型及初始状态。
#----------------------------------------------------------------------------------------------------------
using Distributed, Formatting, DataFrames, CSV #调用并行包，格式化包，数据处理和CSV格式包
addprocs(10; exeflags="--project") # add new processes for parallel
@everywhere using Graphs, SimpleWeightedGraphs, JLD2, LinearAlgebra #图论，加权包，保存,线性代数
@everywhere include("sys.jl") #包含sys.jl文件: 核心计算
@everywhere using .sys  #调用sys.jl文件中的sys模块
@everywhere include("Cluster.jl") #产生cluster
@everywhere using .Cluster
const line = "***********************************************************************************************"
#----------------------------------------------------------------------------------------------------------
# set current work system including transition of states, calcium exchange and buffering, and initial state.
@everywhere function work(net; keys)
    project = keys.project
    A = net.A #邻接矩阵
    numberNodes = net.N #节点数

    #状态及钙离子浓度变化微分方程
    Func = Function[]          # vector for alcium exchange and buffering
    reactions = sys.Reaction[] # vecotr for trnasiiton of states
    # basic states of RyR
    if project == "spark" #钙火花，在JSR与LCC的subspace 
        #设置状态
        numberStates = 4 #单节点状态数
        C = [[1, 0, 0, 0], [0, 0, 1, 0]] #矢量化状态 close-bound close-unbound
        O = [[0, 1, 0, 0], [0, 0, 0, 1]] # open-bound open-unbound
        RyR = sys.RyRstates(numberStates, C, O)
        #transition of states
        v700 = fill(700.0, numberNodes, 1)
        v0 = fill(0, numberNodes, 1)
        v900 = fill(900.0, numberNodes, 1)
        v200000 = fill(200000.0, numberNodes, 1)
        para=split(keys.p.para,"_")
        if para[1]=="rate"
        ratevar=parse(Float64,split(keys.p.para,"_")[2])/10.
        else
            ratevar=1.
        end
        rateCuOu_ratC(S) = ratevar*max.(min.(1.262e-3 * S.Ca .^ 2.8, v700), v0) #ratC
        rateCuOu_ratS(S) = 816.36 ./ (1.0 .+ (86.96^4) ./ S.Ca .^ 4) #ratS
        rateCuOu_ratK(S) = 816.36 ./ (1.0 .+ ((86.96*0.5)^4) ./ S.Ca .^ 4) #ratS
        rateCuOu_sheepC(S) = max.(min.(0.199488 .* S.Ca .^ 2.12, 800.0), 0.0) #ratC
        rateCuOu_sheepS(S) = 1259.1 ./ (1.0 .+ (63.35^2) ./ S.Ca .^ 2) #ratS
        rateCuOu_sheepK(S) = 1259.1 ./ (1.0 .+ (63.35^2) ./ S.Ca .^ 2) #ratS
          if keys.p.para == "ratC"
            rateCuOu = rateCuOu_ratC
        elseif keys.p.para == "ratS"
            rateCuOu = rateCuOu_ratS
        elseif keys.p.para == "ratK"
            rateCuOu = rateCuOu_ratK
        elseif keys.p.para == "sheepC"
            rateCuOu = rateCuOu_sheepC
        elseif keys.p.para == "sheepS"
            rateCuOu = rateCuOu_sheepS
        elseif keys.p.para == "sheepK"
            rateCuOu = rateCuOu_sheepK
        else
            rateCuOu = rateCuOu_ratC
        end
        push!(reactions, sys.Reaction(rateCuOu, C[2], O[2]))
        rateOuCu_ratC(S) = min.(max.(7906 * S.Ca .^ -0.5, v900), v200000) #ratC
        rateOuCu_ratS(S) = fill(1066.8, numberNodes, 1) #ratS
        rateOuCu_ratK(S) = fill(1066.8, numberNodes, 1) #ratS
        rateOuCu_sheepC(S) = 1581.85 .* S.Ca .^ -0.27  #ratC
        rateOuCu_sheepS(S) = fill(810.017, numberNodes, 1) #ratS
        rateOuCu_sheepK(S) = fill(1000.0, numberNodes, 1) #ratS
        if keys.p.para == "ratC"
            rateOuCu = rateOuCu_ratC
        elseif keys.p.para == "ratS"
            rateOuCu = rateOuCu_ratS
        elseif keys.p.para == "ratK"
            rateOuCu = rateOuCu_ratK
        elseif keys.p.para == "sheepC"
            rateOuCu = rateOuCu_sheepC
        elseif keys.p.para == "sheepS"
            rateOuCu = rateOuCu_sheepS
        elseif keys.p.para == "sheepK"
            rateOuCu = rateOuCu_sheepK
        else
            rateOuCu = rateOuCu_ratC
        end
        push!(reactions, sys.Reaction(rateOuCu, O[2], C[2]))
        rateCbOb_ratC(S) = (1.0 / 7.6) .* rateCuOu_ratC(S) #ratC
        rateCbOb_ratS(S) = (1.0 / 7.6) .* rateCuOu_ratS(S) #ratS
        rateCbOb_ratK(S) = (1.0 / 7.6) .* rateCuOu_ratK(S) #ratS
        rateCbOb_sheepC(S) = (1.0 / 7.6) .* rateCuOu_sheepC(S) #ratC
        rateCbOb_sheepS(S) = (1.0 / 7.6) .* rateCuOu_sheepS(S) #ratS
        rateCbOb_sheepK(S) = (1.0 / 7.6) .* rateCuOu_sheepK(S) #ratS
        if keys.p.para == "ratC"
            rateCbOb = rateCbOb_ratC
        elseif keys.p.para == "ratS"
            rateCbOb = rateCbOb_ratS
        elseif keys.p.para == "ratK"
            rateCbOb = rateCbOb_ratK
        elseif keys.p.para == "sheepC"
            rateCbOb = rateCbOb_sheepC
        elseif keys.p.para == "sheepS"
            rateCbOb = rateCbOb_sheepS
        elseif keys.p.para == "sheepK"
            rateCbOb = rateCbOb_sheepK
        else
            rateCbOb = rateCbOb_ratC
        end
        push!(reactions, sys.Reaction(rateCbOb, C[1], O[1]))

        if keys.p.para == "ratC"
            rateObCb = rateOuCu_ratC
        elseif keys.p.para == "ratS"
            rateObCb = rateOuCu_ratS
        elseif keys.p.para == "ratK"
            rateObCb = rateOuCu_ratK
        elseif keys.p.para == "sheepC"
            rateObCb = rateOuCu_sheepC
        elseif keys.p.para == "sheepS"
            rateObCb = rateOuCu_sheepS
        elseif keys.p.para == "sheepK"
            rateObCb = rateOuCu_sheepK
        else
            rateObCb = rateOuCu_ratC
        end
        push!(reactions, sys.Reaction(rateObCb, O[1], C[1]))

        KK = 1000.0
        BCSQN = keys.p.CSQ
        rho(Cas) = @. 5000.0 * Cas^23 / (KK^23 + Cas^23)
        function Mhat(Cas)
            @fastmath rh = @. 5000.0 * Cas^23 / (KK^23 + Cas^23)
            @fastmath M = @. (sqrt(1.0 + 8.0 * rh * BCSQN) - 1.0) / (4.0 * rh * BCSQN)
            return M
        end
        rateCuCb(S) = Mhat(S.Cas) .* (BCSQN / 400.0 / 0.002) #k14
        push!(reactions, sys.Reaction(rateCuCb, C[2], C[1]))
        rateCbCu(S) = fill(8.0, numberNodes, 1) #k41
        push!(reactions, sys.Reaction(rateCbCu, C[1], C[2]))

        rateOuOb(S) = Mhat(S.Cas) .* (BCSQN / 400.0 / 0.002) #k23
        push!(reactions, sys.Reaction(rateOuOb, O[2], O[1]))
        rateObOu(S) = fill(60.8, numberNodes, 1)  #k32
        push!(reactions, sys.Reaction(rateObOu, O[1], O[2]))

        v1 = 0.6e-15 / 2.0 / 9.6485e-2 / 1e-19 #RyR release rate
        tauefflux, taurefill = 1e-5, 0.005   # time constants
        tauefflux, taurefill = 2e-5, 0.005   # time constants

        Ca0, Cas0 = keys.p.Ca0, keys.p.Cas0 # initial state
        Camyo0, CaNSR = keys.p.Ca0, keys.p.Cas0 # resting state
        Buffer1T, Buffer2T, Buffer3T = 24.0, 47.0, 1124

        lmargin = net.lmargin #边界节点的程度
        sumA = sum(A, dims=2)
        ASR = copy(A)
        for i in 1:numberNodes, j in 1:numberNodes
            if net.n[i] != net.n[j]
                ASR[i, j] = 0
            end
        end
        sumASR = sum(ASR, dims=2)

        function Fspark!(S, dt::Float64)
            @fastmath open = @. S.R[:, 2] + S.R[:, 4]
            open[S.LCC.==1] .= 1
            @fastmath JRyR = @. v1 * (S.Cas - S.Ca) * open

            @fastmath JRRcyt = (A * S.Ca .- sumA .* S.Ca) #.* 50. ./ net.Ncs 
            @fastmath JRRcyt = @. JRRcyt * open * 70.0 / (net.Ncs + 20) + JRRcyt * (1 - open)

            #JRRcyt[JRyR.>2e7] .= 0.0
            #@fastmath JRRcyt= @. JRRcyt *(1 - (JRyR/3.11e7))
            @fastmath JRRSR = (ASR * S.Cas .- sumASR .* S.Cas) #.* 50. ./ net.Ncs 

            @fastmath Jefflux = @. (S.Ca - Camyo0) / tauefflux * lmargin #efflux only on margin 
            @fastmath Jrefill = @. (CaNSR - S.Cas) / taurefill

            @fastmath JBuffer1 = @. (100.0 * S.Buffer[1] * S.Ca - 38.0 * (Buffer1T - S.Buffer[1]))
            @fastmath JBuffer2 = @. (115.0 * S.Buffer[2] * S.Ca - 100.0 * (Buffer2T - S.Buffer[2]))
            @fastmath JBuffer3 = @. (115.0 * S.Buffer[3] * S.Ca - 1000.0 * (Buffer3T - S.Buffer[3]))

            @fastmath S.Buffer[1] = @. S.Buffer[1] - JBuffer1 * dt
            @fastmath S.Buffer[2] = @. S.Buffer[2] - JBuffer2 * dt
            @fastmath S.Buffer[3] = @. S.Buffer[3] - JBuffer3 * dt

            @fastmath S.Ca = @. S.Ca + (JRyR + JRRcyt - Jefflux - JBuffer1 - JBuffer2 - JBuffer3) * dt
            KC = 600.0
            x = (4.0 * BCSQN) .* rho(S.Cas)
            @fastmath xp = @. (4.0 * BCSQN * 5000.0 * KK^23 * 23.0) * S.Cas^22 / (KK^23 + S.Cas^23)^2
            @fastmath sqrt12x = @. sqrt(1.0 + 2.0 * x)
            @fastmath pMh = @. (x / sqrt12x - sqrt12x + 1.0) / x^2 * xp
            nc = -20.0 .* Mhat(S.Cas) .+ 35.0
            @fastmath pnc = @. -20 * pMh
            @fastmath beta = @. 1.0 / (1 + (KC * BCSQN * nc + pnc * (S.Cas * KC + S.Cas^2)) / (KC + S.Cas)^2)
            @fastmath S.Cas = @. S.Cas + ((-JRyR + JRRSR) * 0.05 + Jrefill) * beta * dt
        end

        push!(Func, Fspark!)

    end
    #设置初始值
    Ri = zeros(Int64, numberNodes, RyR.numberStates)
    for i in 1:numberNodes
        Ri[i, 1:RyR.numberStates] = C[2] #设置初始状态
        ran = rand()
        if ran < 0.01
            Ri[i, 1:RyR.numberStates] = C[2] #设置初始状态
        end
    end
    Cai = fill(Ca0, numberNodes, 1)  #初始SS钙浓度
    Casi = fill(Cas0, numberNodes, 1) #初始内质网钙浓度
    Camyoi = fill(Camyo0, numberNodes, 1) #初始内质网钙浓度
    Bufferi = [fill(19.0, numberNodes, 1), fill(42.1525, numberNodes, 1), fill(1111.22, numberNodes, 1)]
    LCCi = zeros(Int64, numberNodes) #初始内质网钙浓度
    S = sys.State(Ri, Cai, Casi, Bufferi, Camyoi, LCCi)

    #运行
    t, Ss = sys.evolute(reactions, Func, S, RyR, net, keys)
    return t, Ss
end

#----------------------------------------------------------------------------------------------------------
#主函数
function main()
    println(line, "\n", "Program start...")
    #参数设置 
    println(line, "\n", "Set parameter...")
    pall = CSV.read("p.csv", DataFrame) #读取参数


    for i in findall(pall.onoff .== 1)
        #Ndistribution: 团簇重新产生次数,Nrepeat: 每个团簇重复次数,timeMax: 总时间
        #triggerPeriod 触发（打开RyR或者LCC释放）的间隔,  
        #Spontaneous 为自发钙火花参数，0为关闭，小数为每个时间点RyR打开概率，整数为打开的RyR个数(仅对单团簇,取1)
        #LCC: LCC诱发参数，0为关闭，整数为每触发时间点每个团簇里打开的LCC个数,小数为概率;NoCluster: 多团簇时LCC触发的团簇个数
        #numberNodes 节点数. single 单团簇还是多团簇,  data实验数据还是产生的网络
        p = pall[i, :]
        println(p) #输出参数
        #----------------------------------------------------
        #保存团簇重新产生次数,每个团簇重复次数，和每次保存时间点数(一般0.1ms间隔，+1为最后一个时间点，计算的时间间隔为0.001ms)
        f = open("results/" * String(p.para) * "_res.txt", "w") #打开保存数据的文件
        println(f, p.Ndistribution, ",", p.Nrepeat, ",", Int64(round(p.timeMax * 1e4)) + 1)
        #将参数放到keys元组里
        #plotpoints：画图的间隔数
        #method=Gillespie or Fixed; model=open, spark
        keys = (timeMax=p.timeMax, steps=Int64(round(p.timeMax * 1e6)),
            plotPoints=Int64(round(p.timeMax * 1e4)), triggerPeriod=p.triggerPeriod,
            Spontaneous=p.Spontaneous, LCC=p.LCC, NoCluster=p.NoCluster,
            method="Fixed", project="spark", p=p)
        for icluster in 1:p.Ndistribution
            println(line, "\n", "Generate network:", icluster)
            net = Cluster.random_Cluster_graph(p.numberNodes, p.para, weighted=true, single=p.single, data=p.data, clusterDistance=p.clusterDistance, distance=p.distance) #产生几何网络 
            println("done")
            net = load("results/" * String(p.para) * "_net.jld", "net") #读取网络

            println(f, length(net.sizeV))
            for i10 in 1:Int(p.Nrepeat / (nprocs() - 1))

                println(line, "\n", "Simulation:", i10, "    (Total 2-", nprocs(), " workers, but only worker 2 is shown)")
                net = load("results/" * String(p.para) * "_net.jld", "net") #读取网络

                res = [fetch(@async remotecall(work, iNK + 1, net, keys=keys)) for iNK in 1:nprocs()-1]
                @everywhere GC.gc() #memory release


                for iNK in 1:nprocs()-1
                    #保存结果为txt文件。
                    t = [res[iNK][1], res[iNK][1]]
                    Ss = [[res12.Ca for res12 in res[iNK][2]], [res12.Cas for res12 in res[iNK][2]],
                        [res12.R[:, 2] + res12.R[:, 4] for res12 in res[iNK][2]]]

                    for i in eachindex(t[1])
                        print(f, format(t[1][i], precision=6))
                        Nt = 1
                        for isizeV in eachindex(net.sizeV)
                            size = net.sizeV[isizeV]
                            print(f, ",", sum(Ss[3][i][Nt:Nt+size-1]), ",", sum(Ss[1][i][Nt:Nt+size-1]) / size)
                            Nt = Nt + size
                        end
                        print(f, "\n")
                    end
                    #保存结果为name+.jld文件以便画图。只保存第一个 
                    if icluster == 1 && i10 == 1 && iNK == 1
                        print("Save results for plot...")

                        JLD2.save("results/" * String(p.para) * string(iNK) * "_plot.jld", "t", t, "Ss", Ss)
                        JLD2.save("results/name.jld", "name", String(p.para))
                        #println(name)
                        println("done")
                    end
                end
                res, net = nothing, nothing #memory release
                @everywhere GC.gc() #memory release
            end
        end
        close(f)
    end

    #关闭多余线程
    for i in 1:nprocs()-1
        rmprocs(i + 1)
    end
    println(line, "\n", "program end", "\n", line)
end



@time main() #运行主程序
print(line, "\n", "Plot", "\n")
include("plot.jl") #画图
println("Plot end", "\n", line)

