#----------------------------------------------------------------------------------------------------------
#核心计算
#----------------------------------------------------------------------------------------------------------
module sys
using Distributed, JLD2, Printf, ProgressBars #格式输出包

#----------------------------------------------------------------------------------------------------------
# transitions btween different states of reacteptor
struct Reaction{F<:Function,V<:Vector{Int64}}
    rate::F #transition rate
    stateInitial::V
    stateFinal::V
end
#state of reacteptor and Calcium
mutable struct State{MI<:Matrix{Int64},MF<:Matrix{Float64},V<:Vector{Matrix{Float64}},I<:Vector{Int64}}
    R::MI
    Ca::MF
    Cas::MF
    Buffer::V
    Camyo::MF
    LCC::I
end
#Basic states of RyR
mutable struct RyRstates{I<:Int64,V<:Vector{Vector{Int64}}}
    numberStates::I
    C::V
    O::V
end


#copy state
function Base.copy(s::State)::State
    State(copy(s.R), copy(s.Ca), copy(s.Cas), copy(s.Buffer), copy(s.Camyo), copy(s.LCC))
end

#----------------------------------------------------------------------------------------------------------
#倾向函数 （核心函数）
function propensity(state, left, numberStates)::Float64
    return prod([binomial(state[i], left[i]) for i in 1:numberStates])
end
#选择触发的RyR
function initRyR(keys, net)
    sizeV = net.sizeV
    lmargin = net.lmargin
    No = Int64[]
    Ntemp = 0
    if keys.p.para == "Wave" || keys.p.para == "Wavetest"

        xc, yc, nc = [], [], []

        ii = 1
        for i in net.sizeV
            push!(xc, sum(net.x[ii:ii+i-1]) / i)
            push!(yc, sum(net.y[ii:ii+i-1]) / i)
            push!(nc, 1)
            ii += i
        end

        NNN = 1
        dd0 = 100000
        for i in eachindex(xc)
            dd = sqrt((xc[i] - 1)^2 + (yc[i] - 1)^2)
            if dd < dd0 && sizeV[i] > 40
                dd0 = dd
                NNN = i
            end
        end

        no1 = NNN
        NoCluster = NNN
        Ntemp = sum(sizeV[1:NNN-1])
    else
        no1 = 1
        NoCluster = keys.NoCluster
        if keys.NoCluster > length(net.sizeV)
            NoCluster = length(net.sizeV)
        end
    end

    for size in sizeV[no1:NoCluster]
        if mod(keys.LCC, 1) > 1e-10
            for i in Ntemp+1:Ntemp+size
                ran = rand()
                if ran < keys.LCC
                    push!(No, i)
                    #@show i, Ntemp+1:Ntemp+size
                end
            end
        elseif mod(keys.LCC, 1) < 1e-10

            temp = sortperm(lmargin[Ntemp+1:Ntemp+size]) .+ Ntemp

            io = 0
            NLCC = keys.LCC
            if keys.LCC > size
                NLCC = size
            end
            while io < NLCC && temp != Int[]
                ran = rand(temp)


                io += 1
                push!(No, ran)
                #println(size," ",ran)
            end
        end
        Ntemp += size
    end
    #exit()
    return No
end
#模拟演化
function evolute(reactions::Array{Reaction}, Func::Array{Function}, S::State, RyR::RyRstates, net, keys)
    timeMax, steps, method = keys.timeMax, keys.steps, keys.method
    numberNodes, numberStates = size(S.R)
    if keys.p.para == "LCCIt67i01"
        tmax, iLCC = 6.7, 0.1
    elseif keys.p.para == "LCCIt67i15"
        tmax, iLCC = 6.7, 1.5
    elseif keys.p.para == "LCCIt67i151"
        tmax = 6.7
    elseif keys.p.para == "LCCIt05i05"
        tmax, iLCC = 0.5, 0.5
    elseif keys.p.para == "LCCIt150i05"
        tmax, iLCC = 15.0, 0.5
    elseif keys.p.para == "LCCIt10005i05"
        tmax = 15.
    else
        tmax,iLCC=6.7,0.5
    end

    t = Float64[0.0] # 时间轨迹，t[0]是初始时间
    Ss = sys.State[copy(S)] # fractions of states，n[0] is the initial fractions
    t00 = 0.0 #初始时间
    IP = Float64[] #倾向函数数组，元素为每个节点和反应的转换倾向。下面计算初始的倾向
    for react in reactions, iN in 1:numberNodes
        push!(IP, propensity(view(S.R, iN, 1:numberStates), react.stateInitial, numberStates))
    end
    tplot, tTrigger, tLCC, lTrigger = 0.0, timeMax / steps, 0.0, 1 # for saving the results of fractions
    NLcc = initRyR(keys, net)
    isteps = 0   # 步数
    dt = timeMax / steps #步长
    #以下为演化
    #while isteps<steps && t00<timeMax #满足最大步数或时间则停止
    iter = 0:steps
    if myid() == 2
        iter = ProgressBar(0:steps)
    end
    for isteps in iter
        #计算RyR的状态变化，具体变换机制见main.jl
        #两种方法，一种固定步长Fixed，一种Gillespie算法
        if method == "Fixed"
            t00 += dt #增加步长
            #isteps+=1
            tplot += dt
            iN0 = Int64[] #用于记录状态变化的节点以便后面改变其转换倾向
            ireact = 1
            for react in reactions
                rate0::Array{Float64} = react.rate(S)
                for iN in 1:numberNodes
                    A::Float64 = rate0[iN] * IP[(ireact-1)*numberNodes+iN] #速率乘倾向为反应概率
                    if A * dt > rand(Float64)  #判断在dt时间内可否发生
                        #@show S.R[iN, 1:4]
                        if in(iN, iN0) == false
                            for inumberStates in 1:numberStates #改变状态
                                S.R[iN, inumberStates] += react.stateFinal[inumberStates] - react.stateInitial[inumberStates]
                                #if S.R[iN, inumberStates]!=0 && S.R[iN, inumberStates]!=1
                                #@show iN,inumberStates,react.stateFinal[inumberStates] - react.stateInitial[inumberStates]
                                #end
                            end
                            #@show S.R[iN, 1:4]
                            push!(iN0, iN)
                        end #记录状态改变的节点以便后面改变其转换倾向
                    end
                end
                ireact += 1
                #@show ireact
            end
            #改变倾向
            ireact = 1
            for react in reactions
                for iN in iN0 #记录的节点
                    IP[(ireact-1)*numberNodes+iN] = propensity(view(S.R, iN, 1:numberStates), react.stateInitial, numberStates)  #重新计算倾向
                end
                ireact += 1
            end
        end
        # Gillespie算法
        if method == "Gillespie"
            r1 = rand(Float64)
            r2 = rand(Float64)

            ireact = 1
            rIP = Float64[]
            for react in reactions
                rate0::Array{Float64} = react.rate(S.Ca)
                for iN in 1:numberNodes
                    push!(rIP, rate0[iN] * IP[(ireact-1)*numberNodes+iN])  #计算每个节点和状态的反应概率 这里将节点和状态统一排序
                end
                ireact += 1
            end

            cum = cumsum(rIP) #将得到的反应概率累加
            alpha = last(cum) #总反应概率
            dt = (1.0 / alpha) * log(1.0 / r1) #计算反应时间
            t00 += dt #新时间

            index = searchsortedfirst(cum, alpha * r2) #判断哪一个反应发生
            index_react, index_N = fldmod(index - 1, numberNodes) #有统一排序的序号反推得到反应和节点的序号
            index_react = index_react + 1 #序号需要加1
            index_N = index_N + 1  #序号需要加1
            for inumberStates in 1:numberStates
                S.R[index_N, inumberStates] += reactions[index_react].stateFinal[inumberStates] - reactions[index_react].stateInitial[inumberStates] #改变状态
            end
            ireact = 1
            for react in reactions
                IP[(ireact-1)*numberNodes+index_N] = propensity(S.R[index_N, 1:numberStates], react.stateInitial, numberStates) #改变倾向
                ireact += 1
            end
        end
        #计算连续变化变量，比如钙浓度等具体函数见main.jl
        for func in Func
            func(S, dt)
        end
        #当钙离子浓度变为负值时，强制为0。
        S.Ca[S.Ca.<0] .= 0.0
        S.Cas[S.Ca.<0] .= 0.0

        #触发自发钙火花 
        if keys.Spontaneous != 0
            if abs(tTrigger - keys.triggerPeriod) < dt / 10
                lTrigger = 1
                tTrigger = 0.0
            end

            if lTrigger == 1
                if keys.Spontaneous != 1
                    for i in 1:numberNodes
                        ran = rand()
                        if ran < keys.Spontaneous

                            if S.R[i, 1] == 1
                                S.R[i, 1:numberStates] = RyR.O[1] #设置为初始打开状态
                            end
                            if S.R[i, 3] == 1
                                S.R[i, 1:numberStates] = RyR.O[2] #设置为初始打开状态
                            end
                            ireact = 1
                            for react in reactions
                                IP[(ireact-1)*numberNodes+i] = propensity(view(S.R, i, 1:numberStates), react.stateInitial, numberStates)  #重新计算倾向
                                ireact += 1
                            end
                        end
                    end
                end
                if mod(keys.Spontaneous, 1) < 1e-10 #for single cluster
                    i = rand((1:numberNodes))

                    if S.R[i, 1] == 1
                        S.R[i, 1:numberStates] = RyR.O[1] #设置为初始打开状态
                    end
                    if S.R[i, 3] == 1
                        S.R[i, 1:numberStates] = RyR.O[2] #设置为初始打开状态
                    end
                    ireact = 1
                    for react in reactions
                        IP[(ireact-1)*numberNodes+i] = propensity(view(S.R, i, 1:numberStates), react.stateInitial, numberStates)  #重新计算倾向
                        ireact += 1
                    end
                    #println(t00," ", sum(S.R[:,2].+S.R[:,4]), "   ",tRyR) 
                end
                lTrigger = 0
            end
            tTrigger += dt
        end
        #触发LCC钙火花 
        if keys.LCC != 0

            if abs(tTrigger - keys.triggerPeriod) < dt / 10
                lTrigger = 1
                tTrigger = 0.0
            end

            if lTrigger == 1

                if tLCC < tmax*1e-3

                    for i in NLcc
                        if keys.p.para == "LCCIt67i151"
                            iLCC = (tmax - tLCC) / tmax*15.
                        end
                        if keys.p.para == "LCCIt10005i05"
                            if tLCC % 1.0 < 0.5
                                iLCC=0.5
                            else
                                iLCC=0.
                            end
                        end
                        S.Ca[i] += iLCC * 1e-12 / 2.0 / 9.6485e-2 / 1e-19 * dt  #设置为初始打开状态 #注意i的单位pA和RyR的i单位pA/mM不同
                        S.LCC[i] = 1
                    end
                    tLCC += dt
                else
                    S.LCC[:] .= 0
                    NLcc = initRyR(keys, net)
                    tLCC = 0.0
                    lTrigger = 0

                end
            end
            tTrigger += dt
        end
        #记录数值
        if abs(tplot - timeMax / keys.plotPoints) < timeMax / keys.steps * 0.1
            #@printf(" step: %8.0f   time:%9.5f  Ca:%8.4f  Cas:%8.4f PO:%6.4f \n",
            #    isteps,t00,sum(S.Ca)/numberNodes,  sum(S.Cas)/numberNodes, sum(S.R[:,2])/numberNodes) #显示结果
            push!(Ss, copy(S))  #记录状态
            push!(t, t00) #记录时间
            tplot = 0.0

        end
    end
    return t, Ss
end

end
