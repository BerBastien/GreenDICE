#Desired variables for RESULTS
d_v = ["climatedynamics" "TATM"; "emissions" "E"; "green_naturalcapital" "NC"; "green_naturalcapital" "nonUV"; "grosseconomy" "ES"; "grosseconomy" "YGROSS"; "neteconomy" "S"; "neteconomy" "invNCfrac"; "emissions" "MIU"; "neteconomy" "CPC"; "neteconomy" "YNET"; "grosseconomy" "K"]


#defining functions (start)
function foo()
    A = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    for i = 2:size(d_v)[1]
        B = getdataframe(GreenDICE,Symbol(d_v[i,1]),Symbol(d_v[i,2]))
        A = join(A, B, on= :time)
    end
    scc = DataFrame(time = Int64[], scc = Float64[])
    for ii in 1:60
        year = 2005 + ii * 5
        last_year = 2305
        eta = 1.45
        prtp = prtp = round((1/GreenDICE[:welfare, :rr][2])^(1/5)-1,digits = 4)
        theta = GreenDICE[:welfare,:theta]
        share = GreenDICE[:welfare,:share]
        theta2 = GreenDICE[:welfare,:theta2]
        share2 = GreenDICE[:welfare,:share2]
        Extraton = fill(0.,60)
        set_param!(GreenDICE,:emissions,:ExtraCO2,Extraton)        
        Extracons = fill(0.,60)
        set_param!(GreenDICE,:neteconomy,:ExtraC,Extracons)
        run(GreenDICE)
        C = GreenDICE[:neteconomy, :CPC]
        ES = GreenDICE[:neteconomy, :ESPC]
        nonUV = GreenDICE[:green_naturalcapital, :nonUV]
        l = GreenDICE[:welfare, :l]
        df = [zeros(ii-1)..., (1/(1+prtp)^(t-year) for (i,t) in enumerate(model_years) if year<=t<=last_year)...]
        U = ((((1.0).*((1.0 - share).* C .^ theta) + (share .* ES .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        
        W = sum(U.*df.*l)
        Extraton = fill(0.,60)
        Extraton[ii] = 1
        set_param!(GreenDICE,:emissions,:ExtraCO2,Extraton)
        run(GreenDICE)
        C2 = GreenDICE[:neteconomy, :CPC]
        ES2 = GreenDICE[:neteconomy, :ESPC]
        nonUV2 = GreenDICE[:green_naturalcapital, :nonUV]
        U2 = ((((1.0).*((1.0 - share).* C2 .^ theta) + (share .* ES2 .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV2) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        
        W2 = sum(U2.*df .* l)
        dWdE = - (W2 - W)
        #Now compute dWdC
        Extraton = fill(0.,60)
        set_param!(GreenDICE,:emissions,:ExtraCO2,Extraton)
        Extracons = fill(0.,60)
        Extracons[ii] = 1
        set_param!(GreenDICE,:neteconomy,:ExtraC,Extracons)
        run(GreenDICE)
        ES3 = GreenDICE[:neteconomy, :ESPC]
        nonUV3 = GreenDICE[:green_naturalcapital, :nonUV]
        C3 = GreenDICE[:neteconomy, :CPC]
        U3 = ((((1.0).*((1.0 - share).* C3 .^ theta) + (share .* ES3 .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV3) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        W3 = sum(U3.*df .* l)
        dWdC = (W3 - W)
        
        scc_i = 1000 * dWdE / dWdC
        year_i = 2005 + ii * 5
                push!(scc, [year_i, scc_i])
                end
                A = join(A,scc, on= :time)
                prtp = round((1/GreenDICE[:welfare, :rr][2])^(1/5)-1,digits = 4)
                df = DataFrame(time = A.time, gama3 = ones(60).*GreenDICE[:grosseconomy, :gama3], ratioKNC = ones(60).*GreenDICE[:grosseconomy,:ratioNC],g4 = ones(60).*GreenDICE[:green_naturalcapital,:g4],theta1 = ones(60).*GreenDICE[:welfare,:theta],share1 = ones(60).*GreenDICE[:welfare,:share],climsen = ones(60).*GreenDICE[:climatedynamics,:t2xco2],prtp = ones(60)*prtp,share2 = ones(60).*GreenDICE[:welfare,:share2],theta2 = ones(60).*GreenDICE[:welfare,:theta2],YGreen = GreenDICE[:neteconomy,:YGreen],C = GreenDICE[:neteconomy,:C],damageNC = ones(60).*GreenDICE[:damages,:a4])
                A = join(A,df, on= :time)
    return A
end


function eval_dice(x)
    m = x[1:60]
    s = x[61:end]
    set_param!(GreenDICE,:emissions,:MIU,m)
    set_param!(GreenDICE,:neteconomy,:S,s)
    run(GreenDICE)
    return -GreenDICE[:welfare, :UTILITY]
end


function eval_dice_inv(x)
    m = x[1:60]
    s = x[61:120]
    inv = x[121:180]
    set_param!(GreenDICE,:emissions,:MIU,m)
    set_param!(GreenDICE,:neteconomy,:S,s)
    set_param!(GreenDICE,:neteconomy,:invNCfrac,inv)
    run(GreenDICE)
    
    return -GreenDICE[:welfare, :UTILITY]
end

function choice(a::Array)
    n = length(a)
    idx = rand(1:n)
    return a[idx]
end
#defining functions (end)



# Finding damage parameters (start)
function damage_params(x)
    a4 = x[1]
    #Natural capital parameters
    r = GreenDICE[:grosseconomy,:ratioNC] 
    g3 = GreenDICE[:grosseconomy,:gama3]
    g2 = 0.3 - g3
    g4 = GreenDICE[:green_naturalcapital,:g4]
    #Stock
    C0 = GreenDICE[:neteconomy, :CPC][1]
    ES = GreenDICE[:neteconomy, :ESPC][1]
    K = GreenDICE[:grosseconomy, :K][1]
    L = GreenDICE[:grosseconomy, :l][1]
    YNET = GreenDICE[:neteconomy, :YNET][1]
    YNETpc = YNET* 1000 / L
    s = GreenDICE[:neteconomy, :S][1]
    tfp = GreenDICE[:grosseconomy, :al][1]
    #Utility parameters
    eta = 1.45
    #Consumption
    #C0 = (tfp * L ^ 0.7 + K ^ g2 + (r*K)^g3) * 1000 / L #consumption per capita 
    theta = GreenDICE[:welfare,:theta]
    theta2 = GreenDICE[:welfare,:theta2]
    share = GreenDICE[:welfare,:share]
    share2 = GreenDICE[:welfare,:share2]
    a_d = GreenDICE[:damages,:a_d]  #0.0026686075 #Damage parameter in Nordhaus and Sztorc (2013)
    a_k = GreenDICE[:damages,:a2]
    T=2.5
    C0_d = C0 * (1-a_d * T^2) #Damaged consumption according to DICE 
    U_aggregated = ((((1.0).*((1.0 - share).* ((1-a_d*T^2)*((1-s)*(tfp*L^0.7*K^0.3*1000/L))*(r/(1+0*T^2))^g3) .^ theta) + (share .* (YNETpc * r^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K)^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
    U_explicit = ((((1.0).*((1.0 - share).* ((1-a_k*T^2)*((1-s)*(tfp*L^0.7*K^0.3*1000/L))*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (YNETpc*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
    dif_U = abs(U_aggregated - U_explicit)
    
    return dif_U
end

function pricesNC()
    priceN = DataFrame(time = Int64[], inv = Float64[])
    for ii in 1:60
        year = 2005 + ii * 5
        last_year = 2305
        eta = 1.45
        prtp = prtp = round((1/GreenDICE[:welfare, :rr][2])^(1/5)-1,digits = 4)
        theta = GreenDICE[:welfare,:theta]
        share = GreenDICE[:welfare,:share]
        theta2 = GreenDICE[:welfare,:theta2]
        share2 = GreenDICE[:welfare,:share2]
        ExtraNAsset = fill(0.,60)
        set_param!(GreenDICE,:green_naturalcapital,:ExtraN,ExtraNAsset)       
        run(GreenDICE)
        C = GreenDICE[:neteconomy, :CPC]
        ES = GreenDICE[:neteconomy, :ESPC]
        nonUV = GreenDICE[:green_naturalcapital, :nonUV]
        l = GreenDICE[:welfare, :l]
        df = [zeros(ii-1)..., (1/(1+prtp)^(t-year) for (i,t) in enumerate(model_years) if year<=t<=last_year)...]
        U = ((((1.0).*((1.0 - share).* C .^ theta) + (share .* ES .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        W = sum(U.*df.*l)

        ExtraNAsset = fill(0.,60)
        ExtraNAsset[ii] = 5
        set_param!(GreenDICE,:green_naturalcapital,:ExtraN,ExtraNAsset)   
        run(GreenDICE)
        C2 = GreenDICE[:neteconomy, :CPC]
        ES2 = GreenDICE[:neteconomy, :ESPC]
        nonUV2 = GreenDICE[:green_naturalcapital, :nonUV]
        U2 = ((((1.0).*((1.0 - share).* C2 .^ theta) + (share .* ES2 .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV2) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        W2 = sum(U2.*df .* l)
        dWdN = (W2 - W)


        #Now compute dWdC
        ExtraNAsset = fill(0.,60)
        set_param!(GreenDICE,:green_naturalcapital,:ExtraN,ExtraNAsset)
        ExtraKAsset = fill(0.,60)
        ExtraKAsset[ii] = 5
        set_param!(GreenDICE,:grosseconomy,:ExtraK,ExtraKAsset)
        run(GreenDICE)
        ES3 = GreenDICE[:neteconomy, :ESPC]
        nonUV3 = GreenDICE[:green_naturalcapital, :nonUV]
        C3 = GreenDICE[:neteconomy, :CPC]
        U3 = ((((1.0).*((1.0 - share).* C3 .^ theta) + (share .* ES3 .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV3) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        W3 = sum(U3.*df .* l)
        dWdK = (W3 - W)
        
        pricen_i = dWdN / dWdK
        year_i = 2005 + ii * 5
        push!(priceN, [year_i, pricen_i])
    end
    return priceN 
end

function iterations_damagereduction()
    while mc < mc_max
        priceN = pricesNC()
        set_param!(GreenDICE,:neteconomy,:priceNC,priceN.inv)
        resinv = bboptimize(eval_dice_inv;SearchRange=(0.,1.0), NumDimensions=180, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
        best_candidate(resinv) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(resinv)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(resinv)[61:120])
        set_param!(GreenDICE,:neteconomy,:invNCfrac,best_candidate(resinv)[121:180])
        run(GreenDICE)
        results = foo()
        global Results_investment_iterations = join(Results_investment_iterations,results,on= :time, makeunique = true)
        global mc = mc + 1
        #get results (end)
    end
    return Results_investment_iterations
end
