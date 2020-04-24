using DataFrames
using BlackBoxOptim
using Distributions
using CSVFiles
using CSV


#Different estimations of  the role of natural capital in production
# ALL: Adjusted_tfp_i = [1.0203639621,0.9565217391,1.006993007,1,1.1052631579,1,1.0090909091,0.96,1,1,1,1.0555555556,1.0070422535,1,1,1,1.0526315789,0.8644067797,0.9705882353,1,1,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #26 elasticities
Adjusted_tfp_i = [1.0203639621,0.9565217391,1.006993007,1,1.1052631579,1.0090909091,0.96,1.0555555556,1.0070422535,1.0526315789,0.8644067797,0.9705882353,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities
# only positive numbers: Adjusted_tfp_i = [1.0203639621,1.006993007,1,1.1052631579,1.0090909091,1.0555555556,1.0070422535,1.0526315789,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities

k_nc_i = [44760/15841, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
#[world, low-income, lower-middle income, upper/middle income, high income non/OECD, high income OECD]
# Proportion of NC in initial composite = [15841/(44760+15841), 6421/(1967+6421), 6949/(6531+6949), 18960/(28527+18960), 80104/(59096+80104), 19525/(195929+19525)] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B

share_i = [0.1, 0.01, 0.05, 0.15, 0.2] #invented
#0.1 lower bound invented (almost perfect complements)
#0.5 Sterner and Persson
#1/(1-0.28) lower error range Drupp 2018 limits to substitution
#1/(1-0.57) mean Drupp 2018 limits to substitution
#1/(1-0.86) upper error range Drupp 2018 limits to substitution
#10^10 Standard DICE (perfect substitutes)


#sigma_subs_i = [0.5, 0.1, 1/(1-0.28), 1/(1-0.57), 1/(1-0.86)]
#0.1 lower bound invented (almost perfect complements)
#0.5 Sterner and Persson
#1/(1-0.28) lower error range Drupp 2018 limits to substitution
#1/(1-0.57) mean Drupp 2018 limits to substitution
#1/(1-0.86) upper error range Drupp 2018 limits to substitution
#10^10 Standard DICE (perfect substitutes)

sigma_subs_i = [0.5, 10^10, 1/(1-0.57), 1/(1-0.74),1/(1-0.86),1/(1-0.63),1/(1-0.68),1/(1-0.69),1/(1-0.78),1/(1-0.32),1/(1-0.62),1/(1-0.27),1/(1-0.58),1/(1+0.16),1/(1-0.41),1/(1-0.73),1/(1-0.79)]
#0.5 Sterner and Persson
#following table 1 of Drupp 2018 Limits to substitution
#[Sterner, mean empirical, other empirical estimates in Drupp 2018]


#Desired variables for RESULTS
d_v = ["climatedynamics" "TATM"; "emissions" "E"; "green_naturalcapital" "NC"; "green_naturalcapital" "nonUV"; "grosseconomy" "ES"; "grosseconomy" "YGROSS"; "neteconomy" "S"; "neteconomy" "invNCfrac"; "emissions" "MIU"; "neteconomy" "CPC"; "neteconomy" "YNET"; "grosseconomy" "K"]

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
        ES = GreenDICE[:neteconomy, :ES]
        nonUV = GreenDICE[:green_naturalcapital, :nonUV]
        l = GreenDICE[:welfare, :l]
        df = [zeros(ii-1)..., (1/(1+prtp)^(t-year) for (i,t) in enumerate(model_years) if year<=t<=last_year)...]
        #U = (((((1.0 - share) .* (C .* 1000 ./ l) .^ (1.0 - 1.0 ./ sigma_subs)) + ((share) * (ES.* 1000 ./ l) .^ (1.0 - 1.0 ./ sigma_subs))).^((1 - eta).*(sigma_subs ./ (sigma_subs - 1))) .- 1) ./ (1-eta) .- 1)
        U = ((((1.0).*((1.0 - share).* C .^ theta) + (share .* ES .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        
        W = sum(U.*df.*l)
        Extraton = fill(0.,60)
        Extraton[ii] = 1
        set_param!(GreenDICE,:emissions,:ExtraCO2,Extraton)
        run(GreenDICE)
        C2 = GreenDICE[:neteconomy, :CPC]
        ES2 = GreenDICE[:neteconomy, :ES]
        nonUV2 = GreenDICE[:green_naturalcapital, :nonUV]
        #U2 = (((((1.0 - share) .* (C2 .* 1000 ./ l) .^ (1.0 - 1.0 ./ sigma_subs)) + ((share) * (ES2 .* 1000 ./ l) .^ (1.0 - 1.0 ./ sigma_subs))).^((1 - eta).*(sigma_subs ./ (sigma_subs - 1))) .- 1) ./ (1-eta) .- 1)
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
        ES3 = GreenDICE[:neteconomy, :ES]
        nonUV3 = GreenDICE[:green_naturalcapital, :nonUV]
        C3 = GreenDICE[:neteconomy, :CPC]
        #U3 = (((((1.0 - share) .* (C3 .* 1000 ./ l) .^ (1.0 - 1.0 ./ sigma_subs)) + ((share) * (ES3 .* 1000 ./ l) .^ (1.0 - 1.0 ./ sigma_subs))).^((1 - eta).*(sigma_subs ./ (sigma_subs - 1))) .- 1) ./ (1-eta) .- 1)
        #U3 = ((((1.0 - share2).*((1.0 - share).* C3 .^ theta) + (share .* ES3 .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV3) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        U3 = ((((1.0).*((1.0 - share).* C3 .^ theta) + (share .* ES3 .^ theta)) .^ (theta2 / theta)+ share2 .* (nonUV3) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        W3 = sum(U3.*df .* l)
        dWdC = (W3 - W)
        
        scc_i = 1000 * dWdE / dWdC
        year_i = 2005 + ii * 5
                push!(scc, [year_i, scc_i])
                end
                A = join(A,scc, on= :time)
                prtp = round((1/GreenDICE[:welfare, :rr][2])^(1/5)-1,digits = 4)
                df = DataFrame(time = A.time, gama3 = ones(60).*GreenDICE[:grosseconomy, :gama3], ratioKNC = ones(60).*GreenDICE[:grosseconomy,:ratioNC],XX = ones(60).*0,theta1 = ones(60).*GreenDICE[:welfare,:theta],share1 = ones(60).*GreenDICE[:welfare,:share],climsen = ones(60).*GreenDICE[:climatedynamics,:t2xco2],prtp = ones(60)*prtp,share2 = ones(60).*GreenDICE[:welfare,:share2],theta2 = ones(60).*GreenDICE[:welfare,:theta2],YGreen = GreenDICE[:neteconomy,:YGreen],C = GreenDICE[:neteconomy,:C],damageNC = ones(60).*GreenDICE[:damages,:a4])
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

function eval_DICE(x)
    m = x[1:60]
    s = x[61:end]
    set_param!(DICE,:emissions,:MIU,m)
    set_param!(DICE,:neteconomy,:S,s)
    run(DICE)
    return -DICE[:welfare, :UTILITY]
end

function eval_dice_multi(x)
    m = x[1:60]
    s = x[61:120]
    set_param!(GreenDICE,:emissions,:MIU,m)
    set_param!(GreenDICE,:neteconomy,:S,s)
    run(GreenDICE)
    
    return (-GreenDICE[:neteconomy,:C][19],-GreenDICE[:green_naturalcapital,:ES][19])
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


function f()
    paretopoints = DataFrame(x = Float64[], y = Float64[])
    for x in pareto_frontier(res)
    X = fitness(x)[1]
    Y = fitness(x)[2]
    push!(paretopoints, [X, Y])
    end
    return paretopoints
end


function choice(a::Array)
    n = length(a)
    idx = rand(1:n)
    return a[idx]
end



#Vectors of estimations (start)
    Adjusted_tfp_i = [1.0203639621,0.9565217391,1.006993007,1,1.1052631579,1.0090909091,0.96,1.0555555556,1.0070422535,1.0526315789,0.8644067797,0.9705882353,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities
    k_nc_i = [44760/15841, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
    latfp = log.(Adjusted_tfp_i)
    lknc = log.(k_nc_i)
    lknc_inv = 1 ./ lknc
    gama3_iM = latfp * lknc_inv'
    gama3_i = [gama3_iM[1:end,1];gama3_iM[1:end,2];gama3_iM[1:end,3];gama3_iM[1:end,4];gama3_iM[1:end,5];gama3_iM[1:end,6]]
    gama3_i = gama3_i[gama3_i.>0]
    gama3_i = gama3_i[gama3_i.<0.2999999]
    share1_i = [0.05, 0.1, 0.15] #invented 
    share2_i = [0.05, 0.1, 0.15] #invented 
    theta1_i = [0.74,0.86,0.68,0.69,0.78,0.32,0.62,0.27,0.58,-0.16,0.41,-0.1,0.73,0.79,0.76,0.80]
    theta2_i = [0.63,0.78,0.62,0.27]
    damage_i = [0.00480515,0.01604]#0 damages in NC; 0.00181 25% of market damages; 0.01604 according to S&P
    prtp_i = [0.001,0.015,0.03]
#Vectors of estimations (end)

#Distributions (start)
    cs_lnd = LogNormal(exp(3.2), 0.12) #based on Roe and Baker 2007, using the mean of Nordhaus
    prtp_ud = Uniform(0.001, 0.03) 
    gama3_nd = Normal(log(1.0203639621)/log(44760/15841), std(gama3_i)) #normal distributions
    k_nc_nd = Normal(44760/15841,std(k_nc_i))
#Distributions (end) 




## Explanation of the parameters (START)
##
#set_param!(GreenDICE,:welfare,:utility_path,1)   # Is ES in Utility Function?
#set_param!(GreenDICE,:damages,:path,1)           # Damage path: 0 standard DICE; 1 damage to NC; 2 damate to ES
#set_param!(GreenDICE,:grosseconomy,:greenGDP,1)  # Is GDP in the production function?
#set_param!(GreenDICE,:welfare,:share,0.1)                  # Share of ES in the Utility function
#set_param!(GreenDICE,:welfare,:share2,0)                  # Share of nonUV in the Utility function
# 0.1 according to Sterner and Persson
# 0 to return to standard DICE
#set_param!(GreenDICE,:grosseconomy,:share,0.1)                  # Share of ES in the Utility function
#set_param!(GreenDICE,:welfare,:sigma_subs,1/(1-0.57))             # Elasticity of substitution
# 0.5 according to Sterner and Persson
# 10^100 to return to standard DICE 
# 1/(1-0.57) According to Drupp 2018 Limits to substitution
#set_param!(GreenDICE,:welfare,:theta,0.57)             # Elasticity of substitution theta1 of ES to C
#set_param!(GreenDICE,:welfare,:theta2,0.57)             # Elasticity of substitution theta 2 of nonUV to UV
#set_param!(GreenDICE,:damages,:a4,0.00480515)                # NC damage coefficient
# Damages to NC according to Drupp and Hansel (2018): 0.01604
# Damages as if taken out the 25%: 0.00480515
#set_param!(GreenDICE,:damages,:a5,0)                # ES damage coefficient
# Damages to ES according to Drupp (2018): 0.01604
# Damages as if taken out the 25%: 0.00480515
#set_param!(GreenDICE,:damages,:a2,0.00181)                 # YGROSS damage coefficient
# Standard DICE: 0.00236; Minus non-market damages: 0.00181
#K_NC = 44760/15841
#set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   # Ratio of NC to K0 
#set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)   # Ratio of NC to K0 
# Produced Capital vs. Natural Capital (K/NC)
# Default value (World estimate): 27/9 Ref: Changing Wealth of Nations World Bank
# Other possibilities - Rich countries: 28/3; Poor countries: 11/10 Ref: Changing Wealth of Nations World Bank, Table 2.2
# Other possibilities: [world, low-income, lower-middle income, upper/middle income, high income non/OECD, high income OECD]
# K_NC = [44760/15841, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
#Adjusted_tfp = 1.0203639621
#Adjusted TFP once NC is explicitely in the production function
# Default value (Weighted average of estimates): 1.0203639621
# Other posssibilities (individual countries): [0.9565217391,1.006993007,1,1.1052631579,1,1.0090909091,...
#    0.96,1,1,1,1.0555555556,1.0070422535,1,1,1,1.0526315789,0.8644067797,...
#    0.9705882353,1,1,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716,1.0203639621,1.0203639621
#elasticity_nc = log(Adjusted_tfp)/log(K_NC)
#set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
   # set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)
             # Elasticity of YGROSS to NC
##
## Esxplanation of the parameters (END)


##
## GreenDICE with exogenous investments
    set_param!(GreenDICE,:welfare,:utility_path,1)
    set_param!(GreenDICE,:damages,:path,1)        
    set_param!(GreenDICE,:grosseconomy,:greenGDP,1) 
    set_param!(GreenDICE,:welfare,:share,0.1)
    set_param!(GreenDICE,:welfare,:share2,0.1)
    set_param!(GreenDICE,:grosseconomy,:share,0.1) 
    set_param!(GreenDICE,:welfare,:sigma_subs,1/(1-0.57)) 
    set_param!(GreenDICE,:welfare,:theta,mean(theta1_i)) 
    set_param!(GreenDICE,:welfare,:theta2,mean(theta2_i))
    set_param!(GreenDICE,:damages,:a4,0.00480515)
    set_param!(GreenDICE,:neteconomy,:a4,0.00480515)
    set_param!(GreenDICE,:damages,:a5,0)
    set_param!(GreenDICE,:damages,:a2,0.00181) 
    K_NC = 44760/15841
    set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    Adjusted_tfp = 1.0203639621
    elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)
    ExtraNAsset = fill(0.,60)
    set_param!(GreenDICE,:green_naturalcapital,:ExtraN,ExtraNAsset)  
    ExtraKAsset = fill(0.,60)
    set_param!(GreenDICE,:grosseconomy,:ExtraK,ExtraKAsset)
    set_param!(GreenDICE,:neteconomy,:priceNC,fill(10e100,60))
    set_param!(GreenDICE,:neteconomy,:invNCfrac,fill(0.,60))
    
    run(GreenDICE)

    global Results_spread_combined_mcs_opt = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    while mc < 20
        resinv = bboptimize(eval_dice_inv;SearchRange=(0.,1.0), NumDimensions=180, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=100_000)
                    best_candidate(resinv) # optimal vector of miu emissions trajectories
                    set_param!(GreenDICE,:emissions,:MIU,best_candidate(resinv)[1:60])
                    set_param!(GreenDICE,:neteconomy,:S,best_candidate(resinv)[61:120])
                    set_param!(GreenDICE,:neteconomy,:invNCfrac,best_candidate(resinv)[121:180])
                    run(GreenDICE)
                    results = foo()
                    
                global Results_spread_combined_mcs_opt = join(Results_spread_combined_mcs_opt,results,on= :time, makeunique = true)
                global mc = mc + 1
    end

                    CSV.write("C:/Users/bastien/Documents/DICE/MimiGreenDICE/Results/conservation.csv",Results_spread_combined_mcs_opt)

