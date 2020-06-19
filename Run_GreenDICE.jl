using DataFrames
using BlackBoxOptim
using Distributions
using CSVFiles
using CSV

#setting up (start)
    #Do you wnat to perform optimization of welfare?
    welf_opt = true
    spread_combined_mcs_opt=true
    combination_NC_tfp = true
    investments_NC_mcs = true
    spread_combined_mcs_opt_inv=true
    sensitivity_per_parameter_opt = false
    uncertainty_per_parameter = false

    uncertainty_per_parameter_opt = false

    #Do you want to plot pareto frontier?
    par_front = false

    uncertainty_montecarlo = false
    uncertainty_analysis = false

    #Do you want to run GreenDICE with different values for natural capital estimations and total factor productivity asjustments?
    uncertainty_analysis = false


    #Do you want to run GreenDICE with different values for the share of non-use values in utility functions and elasticities of substitution?
    uncertainty_analysis_utility = false

    #Different estimations of  the role of natural capital in production
    # ALL Adjusted_tfp_i = [1.014353,0.9565217391,1.006993007,1,1.1052631579,1,1.0090909091,0.96,1,1,1,1.0555555556,1.0070422535,1,1,1,1.0526315789,0.8644067797,0.9705882353,1,1,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #26 elasticities
    Adjusted_tfp_i = [1.014353,0.9565217391,1.006993007,1,1.1052631579,1.0090909091,0.96,1.0555555556,1.0070422535,1.0526315789,0.8644067797,0.9705882353,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities
    # only positive numbers Adjusted_tfp_i = [1.014353,1.006993007,1,1.1052631579,1.0090909091,1.0555555556,1.0070422535,1.0526315789,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities

    k_nc_i = [3.87, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
    #[world, low-income, lower-middle income, upper/middle income, high income non/OECD, high income OECD]
    #
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
#setting up (end)

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

    function eval_DICE(x)
        m = x[1:60]
        s = x[61:end]
        set_param!(DICE,:emissions,:MIU,m)
        set_param!(DICE,:neteconomy,:S,s)
        run(DICE)
        return -DICE[:welfare, :UTILITY]
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
        C = GreenDICE[:neteconomy, :CPC][1]
        ES = GreenDICE[:neteconomy, :ES][1]
        K = GreenDICE[:grosseconomy, :K][1]
        L = GreenDICE[:grosseconomy, :l][1]
        tfp = GreenDICE[:grosseconomy, :al][1]
        #Utility parameters
        eta = 1.45
        #Consumption
        C0 = tfp * L ^ 0.7 + K ^ g2 + (r*K)^g3
        theta = GreenDICE[:welfare,:theta]
        theta2 = GreenDICE[:welfare,:theta2]
        share = GreenDICE[:welfare,:share]
        share2 = GreenDICE[:welfare,:share2]
        a_d = GreenDICE[:damages,:a_d]  #0.0026686075 #Damage parameter in Nordhaus and Sztorc (2013)
        a_k = GreenDICE[:damages,:a2]
        T=2.5
        C0_d = C0 * (1/(1+a_d * T^2)) #Damaged consumption according to DICE 
        U_aggregated = ((((1.0).*((1.0 - share).* C0_d .^ theta) + (share .* (C0 * r^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K)^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        #a = (((1 + a_k * T^2)(r/(1+a4 * T^2))^g3) / T^2)
        #U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+a*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+a_k*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        #U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+(((1 + a_k * T^2)*(r/(1+a4 * T^2))^g3) / T^2)*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        dif_U = abs(U_aggregated - U_explicit)
        
        return dif_U
    end

    function damage_params_usevalues(x)
        a4 = x[1]
        #Natural capital parameters
        r = GreenDICE[:grosseconomy,:ratioNC] 
        g3 = GreenDICE[:grosseconomy,:gama3]
        g2 = 0.3 - g3
        g4 = GreenDICE[:green_naturalcapital,:g4]
        #Stock
        C = GreenDICE[:neteconomy, :CPC][1]
        ES = GreenDICE[:neteconomy, :ES][1]
        K = GreenDICE[:grosseconomy, :K][1]
        L = GreenDICE[:grosseconomy, :l][1]
        tfp = GreenDICE[:grosseconomy, :al][1]
        #Utility parameters
        eta = 1.45
        #Consumption
        C0 = tfp * L ^ 0.7 + K ^ g2 + (r*K)^g3
        theta = GreenDICE[:welfare,:theta]
        theta2 = GreenDICE[:welfare,:theta2]
        share = GreenDICE[:welfare,:share]
        share2 = GreenDICE[:welfare,:share2]      
        a_k = GreenDICE[:damages,:a2] #0.000708454 From Howard and Sterner (2017) or 
        a_d =  GreenDICE[:damages,:ad] 
        0.0026686075 #Damage parameter in Nordhaus and Sztorc (2013)
        T=3
        C0_d = C0 * (1/(1+a_d * T^2)) #Damaged consumption according to DICE 
        U_aggregated = ((((1.0).*((1.0 - share).* C0_d .^ theta) + (share .* (C0 * r^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((0)^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        #a = (((1 + a_k * T^2)(r/(1+a4 * T^2))^g3) / T^2)
        #U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+a*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+a_k*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((0/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        #U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+(((1 + a_k * T^2)*(r/(1+a4 * T^2))^g3) / T^2)*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        dif_U = abs(U_aggregated - U_explicit)
        
        return dif_U
    end

    function damage_params_market(x)
        a4 = x[1]
        #Natural capital parameters
        r = GreenDICE[:grosseconomy,:ratioNC] 
        g3 = GreenDICE[:grosseconomy,:gama3]
        g2 = 0.3 - g3
        g4 = GreenDICE[:green_naturalcapital,:g4]
        #Stock
        C = GreenDICE[:neteconomy, :CPC][1]
        ES = GreenDICE[:neteconomy, :ES][1]
        K = GreenDICE[:grosseconomy, :K][1]
        L = GreenDICE[:grosseconomy, :l][1]
        tfp = GreenDICE[:grosseconomy, :al][1]
        #Utility parameters
        eta = 1.45
        #Consumption
        C0 = tfp * L ^ 0.7 + K ^ g2 + (r*K)^g3
        theta = GreenDICE[:welfare,:theta]
        theta2 = GreenDICE[:welfare,:theta2]
        share = GreenDICE[:welfare,:share]
        share2 = GreenDICE[:welfare,:share2]      
        a_k = GreenDICE[:damages,:a2] #0.000708454 From Howard and Sterner (2017) or 
        a_d = 0.0026686075 #Damage parameter in Nordhaus and Sztorc (2013)
        T=3
        C0_d = C0 * (1/(1+a_d * T^2)) #Damaged consumption according to DICE 
        U_aggregated = ((((1.0).*((1.0 - share).* C0_d .^ theta) + (share .* (0 * r^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((0)^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        #a = (((1 + a_k * T^2)(r/(1+a4 * T^2))^g3) / T^2)
        #U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+a*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+a_k*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((0/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        #U_explicit = ((((1.0).*((1.0 - share).* ((1/(1+(((1 + a_k * T^2)*(r/(1+a4 * T^2))^g3) / T^2)*T^2))*C0*(r/(1+a4*T^2))^g3) .^ theta) + (share .* (C0*(r/(1+a4*T^2))^g2) .^ theta)) .^ (theta2 / theta)+ share2 .* ((r*K/(1+a4*T^2))^g4) .^ theta2) .^ ((1 - eta) / theta2) .- 1) ./ (1 - eta) .-1 #nested utility function with constant elasticities of substiution
        dif_U = abs(U_aggregated - U_explicit)
        
        return dif_U
    end
# Finding damage parameters (end)

#Vectors of estimations (start)
    Adjusted_tfp_i = [1.014353,0.9565217391,1.006993007,1,1.1052631579,1.0090909091,0.96,1.0555555556,1.0070422535,1.0526315789,0.8644067797,0.9705882353,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities
    k_nc_i = [3.87, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
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
    damagek_i = [1,2,3]
    damaged_i = [0.0026686075,0.0026686075, (0.0026686075/1.25)*2]
    prtp_i = [0.001,0.015,0.03]
    damagek3 = 0.624 #From Nordhaus 1994 Survey
#Vectors of estimations (end)

#Distributions (start)
    cs_lnd = LogNormal(log(3.2), 0.12) #based on Roe and Baker 2007, using the mean of Nordhaus
    prtp_ud = Uniform(0.001, 0.03) 
    tfp_param = Normal(1.014353,0.028783) #based on wighted adjusted TFP by country GDP (see gamma3_computation.R)
    gama3_nd = Normal(log(1.014353)/log(3.87), std(gama3_i)) #normal distributions
    k_nc_nd = Normal(3.87,2.11)
    damagek1_nd = Normal(0.151, 0.147) #From regression using Howard and Sterner (2017) data
    damagek2_ud = Uniform(0.3,0.5) #From Howard and Sylvan (2015)
    atfp_nd = Normal(1.0144,0.028783) #From Brandt et al (2017) weighted mean and standard deviation
    gamma4_ud = Uniform(0.6,1)
    
#Distributions (end) 


## Explanation of the parameters (START)
    #set_param!(GreenDICE,:welfare,:utility_path,1)   # Is ES in Utility Function?
    #set_param!(GreenDICE,:damages,:path,1)           # Damage path: 0 standard DICE; 1 damage to NC; 2 damate to ES
    #set_param!(GreenDICE,:grosseconomy,:greenGDP,1)  # Is GDP in the production function?
    #set_param!(GreenDICE,:welfare,:share,0.1)                  # Share of ES in the Utility function
    #set_param!(GreenDICE,:welfare,:share2,0.1)                  # Share of nonUV in the Utility function
    # 0.1 according to Sterner and Persson
    # 0 to return to standard DICE
    #set_param!(GreenDICE,:grosseconomy,:share,0.1)                  # Share of ES in the Utility function
    #set_param!(GreenDICE,:welfare,:sigma_subs,1/(1-0.57))             # Elasticity of substitution
    # 0.5 according to Sterner and Persson
    # 10^100 to return to standard DICE 
    # 1/(1-0.57) According to Drupp 2018 Limits to substitution
    #set_param!(GreenDICE,:welfare,:theta,0.57)             # Elasticity of substitution theta1 of ES to C
    #set_param!(GreenDICE,:welfare,:theta2,0.57)             # Elasticity of substitution theta 2 of nonUV to UV
    #
    #set_param!(GreenDICE,:damages,:a4,0.00480515)
    #set_param!(GreenDICE,:neteconomy,:a4,0.00480515)                # NC damage coefficient
    # Damages to NC according to Drupp and Hansel (2018): 0.01604
    # Damages as if taken out the 25%: 0.00480515
    #set_param!(GreenDICE,:damages,:a5,0)                # ES damage coefficient
    # Damages to ES according to Drupp (2018): 0.01604
    # Damages as if taken out the 25%: 0.00480515
    #set_param!(GreenDICE,:damages,:a2,0.00181)                 # YGROSS damage coefficient
    # Standard DICE: 0.00236; Minus non-market damages: 0.00181
    #K_NC = 3.87
    #set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   # Ratio of NC to K0 
    #set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)   # Ratio of NC to K0 
    # Produced Capital vs. Natural Capital (K/NC)
    # Default value (World estimate): 27/9 Ref: Changing Wealth of Nations World Bank
    # Other possibilities - Rich countries: 28/3; Poor countries: 11/10 Ref: Changing Wealth of Nations World Bank, Table 2.2
    # Other possibilities: [world, low-income, lower-middle income, upper/middle income, high income non/OECD, high income OECD]
    # K_NC = [3.87, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
    #Adjusted_tfp = 1.014353
    #Adjusted TFP once NC is explicitely in the production function
    # Default value (Weighted average of estimates): 1.014353
    # Other posssibilities (individual countries): [0.9565217391,1.006993007,1,1.1052631579,1,1.0090909091,...
    #    0.96,1,1,1,1.0555555556,1.0070422535,1,1,1,1.0526315789,0.8644067797,...
    #    0.9705882353,1,1,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716,1.014353,1.014353];
    #set_param!(GreenDICE,:grosseconomy,:atfp,Adjusted_tfp)             # Elasticity of YGROSS to NC
    #elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    #set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)             # Elasticity of YGROSS to NC
    ##
## Esxplanation of the parameters (END)


##
## GreenDICE main specification (START)
        #set parameters of specification (start)
        set_param!(GreenDICE,:welfare,:utility_path,1)
        set_param!(GreenDICE,:damages,:path,1)        
        set_param!(GreenDICE,:grosseconomy,:greenGDP,1) 
        set_param!(GreenDICE,:welfare,:share,0.1)
        set_param!(GreenDICE,:welfare,:share2,0.1)
        set_param!(GreenDICE,:grosseconomy,:share,0.1) 
        set_param!(GreenDICE,:welfare,:sigma_subs,1/(1-0.57)) 
        set_param!(GreenDICE,:welfare,:theta,mean(theta1_i)) 
        set_param!(GreenDICE,:welfare,:theta2,mean(theta2_i))
        set_param!(GreenDICE,:damages,:a5,0)
        K_NC = 3.878362
        set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
        set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
        Adjusted_tfp = 1.014353
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
       
        
        # option 1
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        # # option 2
        # a_k = 0.0026686075/1.25 #Corresponding to market damages only, 25% of Norduas and Sztorc (2013)
        # a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        # # option 3
        # a_k = 0.0026686075/1.25 #Corresponding to market damages only, 25% of Norduas and Sztorc (2013)
        # a_d = (0.0026686075/1.25)*2 #Corresponding to total aggregated damages, which account for another 100% according to Stern Review
        
        set_param!(GreenDICE,:damages,:a2,a_k)
        set_param!(GreenDICE,:damages,:a_d,a_d) 
        set_param!(GreenDICE,:damages,:a4,0)
        set_param!(GreenDICE,:neteconomy,:a4,0)

        # r_nk = 15841 / (15841+44760)
        # r_kn = 44760 / (15841+44760)
        # set_param!(GreenDICE,:grosseconomy,:r_kn,r_kn)
        # set_param!(GreenDICE,:green_naturalcapital,:r_nk,r_nk)
        # r_nk = 0.1
        # g3 = -log(1/((1-r_nk)^0.3))/(log(1-r_nk)-log(r_nk))
        # ((1-r_nk)^(0.3-g3)) * r_nk^(g3)
        # k_in_K = 1
        # set_param!(GreenDICE,:grosseconomy,:k_in_K,k_in_K)
        # set_param!(GreenDICE,:green_naturalcapital,:k_in_K,k_in_K)
        # if k_in_K ==1
        # set_param!(GreenDICE,:grosseconomy,:gama3,g3)
        # set_param!(GreenDICE,:neteconomy,:gama3,g3)
        # else
            # set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
            # set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
        # end
        g4 = 0.8
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)

        #Finding the damage parameters (start)
            Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
            a_4 = best_candidate(Resd) #
            set_param!(GreenDICE,:damages,:a4,a_4[1])
            set_param!(GreenDICE,:neteconomy,:a4,a_4[1])
        #Finding the damage parameters (end)
    #set parameters of specification (end)


        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
        run(GreenDICE)
        results = foo()
        CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV.csv",results)



    # if spread_combined_mcs_opt==true
    #     mc_max = 2000
        
    #     global Results_spread_combined_mcs_opt = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     global mc = 0
    #     while mc < mc_max


    #         #Draw 8 parameters (start)
    #             D_d = choice(damaged_i)
    #             D_k = choice(damagek_i)
    #             K_NC = rand(k_nc_nd)
    #             Share = choice(share1_i)
    #             Share2 = choice(share2_i)
    #             gamma4 = rand(gamma4_ud)
    #             Theta = choice(theta1_i)
    #             Theta2 = choice(theta2_i)
    #             atfp = rand(atfp_nd)
    #         #Draw 8 parameters (end)

    #         #choose damage parameters (start)
    #             if D_k == 1
    #                 perck = rand(damagek1_nd)/0.222 #Percentage of damages corresponding to market impacts
    #             end
    #             if D_k == 2
    #                 perck = rand(damagek2_ud)
    #             end
    #             if D_k == 3
    #                 perck = damagek3
    #             end
    #             a_k = perck * D_d
    #         #choose damage parameters (end)

    #             if K_NC < 0
    #                 continue
    #             end
    #             if atfp < 1
    #                 continue
    #             end
    #             elasticity_nc = log(atfp) / log(K_NC)
    #             if elasticity_nc < 0
    #                 continue #skip to the next step in the for loop
    #             end
    #             if elasticity_nc > 0.2999999
    #                 continue #skip to the next step in the for loop
    #             end
    #         #conditions (end)

    #         #set parameters (start)
    #             set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   
    #             set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)    
    #             set_param!(GreenDICE,:welfare,:share,Share)  
    #             set_param!(GreenDICE,:welfare,:share2,Share2)   
    #             set_param!(GreenDICE,:green_naturalcapital,:g4,gamma4)   
    #             set_param!(GreenDICE,:welfare,:theta,Theta)    
    #             set_param!(GreenDICE,:welfare,:theta2,Theta2)          
    #             set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)   
    #             set_param!(GreenDICE,:damages,:a_d,D_d)
    #             set_param!(GreenDICE,:damages,:a2,a_k)
    #             run(GreenDICE)
    #             Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
    #             a_4 = best_candidate(Resd) 
    #             set_param!(GreenDICE,:damages,:a4,a_4[1])
    #             set_param!(GreenDICE,:neteconomy,:a4,a_4[1]) 
    #         #set parameters (end)
            
    #         run(GreenDICE)
    #         #optimize
    #             res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
    #         #optimize

    #         #get results (start)
    #             best_candidate(res) # optimal vector of miu emissions trajectories
    #             set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #             set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #             run(GreenDICE)
    #             results = foo()
    #             global Results_spread_combined_mcs_opt = join(Results_spread_combined_mcs_opt,results,on= :time, makeunique = true)
    #             global mc = mc + 1
    #         #get results (end)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_spread_mcs_opt.csv",Results_spread_combined_mcs_opt)
    #     #going back to preferred specification (start)
    #         set_param!(GreenDICE,:welfare,:utility_path,1)
    #         set_param!(GreenDICE,:damages,:path,1)        
    #         set_param!(GreenDICE,:grosseconomy,:greenGDP,1) 
    #         set_param!(GreenDICE,:welfare,:share,0.1)
    #         set_param!(GreenDICE,:welfare,:share2,0.1)
    #         set_param!(GreenDICE,:grosseconomy,:share,0.1) 
    #         set_param!(GreenDICE,:welfare,:sigma_subs,1/(1-0.57)) 
    #         set_param!(GreenDICE,:welfare,:theta,mean(theta1_i)) 
    #         set_param!(GreenDICE,:welfare,:theta2,mean(theta2_i))
    #         set_param!(GreenDICE,:damages,:a5,0)
    #         K_NC = 3.878362
    #         set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    #         set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    #         Adjusted_tfp = 1.014353
    #         elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    #         set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    #         set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
    #         a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
    #         k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
    #         mortality_perc = 0.2 #Percentage corresponding to mortality damages
    #         a_k = a_d * (k_perc + mortality_perc)
    #         set_param!(GreenDICE,:damages,:a2,a_k)
    #         set_param!(GreenDICE,:damages,:a_d,a_d) 
    #         set_param!(GreenDICE,:damages,:a4,0)
    #         set_param!(GreenDICE,:neteconomy,:a4,0)
    #         g4 = 0.8
    #         set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
    #         run(GreenDICE)
    #         Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
    #         a_4 = best_candidate(Resd) #
    #         set_param!(GreenDICE,:damages,:a4,a_4[1])
    #         set_param!(GreenDICE,:neteconomy,:a4,a_4[1])
    #     #going back to preferred specification (end)
    #     run(GreenDICE)
    # end

    # combination_NC_tfp = true
    # if combination_NC_tfp == true
    #     global Results_combination_NC_tfp = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     mc_max = 20
    #     global mc = 0
    #     while mc < mc_max
    #             K_NC = rand(k_nc_nd)
    #             atfp = rand(atfp_nd)

    #             if K_NC < 0
    #                 continue
    #             end
    #             if atfp < 1
    #                 continue
    #             end
    #             elasticity_nc = log(atfp) / log(K_NC)
    #             if elasticity_nc < 0
    #                 continue #skip to the next step in the for loop
    #             end
    #             if elasticity_nc > 0.2999999
    #                 continue #skip to the next step in the for loop
    #             end
                            
    #             set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   # Ratio of NC to K0 
    #             set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)   # Ratio of NC to K0 
    #             set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    #             set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
    #             run(GreenDICE)
    #             res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
    #             best_candidate(res)
    #             set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #             set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #             run(GreenDICE)
    #             results = foo()
    #             global Results_combination_NC_tfp = join(Results_combination_NC_tfp,results,on= :time, makeunique = true)
    #         end
    #     results = Results_combination_NC_tfp
    #     #back to standard for future runs
    #     K_NC = 3.878362
    #     set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    #     set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    #     Adjusted_tfp = 1.014353
    #     elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    #     set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    #     set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_combination_NC_tfp.csv",results)
    #     run(GreenDICE)
    # end

    # sensitivity_per_parameter_opt = true
    # if sensitivity_per_parameter_opt == true
    #     mc_max = 20
        
    #     global Results_uncertainty_damage = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     global mc = 0
    #     mc_max = 20
    #     while mc < mc_max
    #         D_d = choice(damaged_i)
    #         D_k = choice(damagek_i)
    #         if D_k == 1
    #             perck = rand(damagek1_nd)/0.222 #Percentage of damages corresponding to market impacts
    #         end
    #         if D_k == 2
    #             perck = rand(damagek2_ud)
    #         end
    #         if D_k == 3
    #             perck = damagek3
    #         end
    #         a_k = perck * D_d
    #         set_param!(GreenDICE,:damages,:a_d,D_d)
    #         set_param!(GreenDICE,:damages,:a2,a_k)
    #         run(GreenDICE)
    #         Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
    #         a_4 = best_candidate(Resd) 
    #         set_param!(GreenDICE,:damages,:a4,a_4[1])
    #         set_param!(GreenDICE,:neteconomy,:a4,a_4[1]) 
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_damage = join(Results_uncertainty_damage,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_damage.csv",Results_uncertainty_damage)
    #     set_param!(GreenDICE,:damages,:a4,0.050642536)
    #     run(GreenDICE)
        
    #     global Results_uncertainty_prtp = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     for prtpi in (1:size(prtp_i)[1])
    #         prtp = prtp_i[prtpi]
    #         RR = [1/(1+prtp)^(t*5) for t in 0:59]
    #         set_param!(GreenDICE,:welfare,:rr,RR)   # Ratio of NC to K0 
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_prtp = join(Results_uncertainty_prtp,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_prtp.csv",Results_uncertainty_prtp)
    #     RR = [1/(1+0.015)^(t*5) for t in 0:59]
    #     set_param!(GreenDICE,:welfare,:rr,RR)   # Ratio of NC to K0 
    #     run(GreenDICE)

    #     global Results_uncertainty_cs = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     global mc = 0
    #     mc_max = 20
    #     while mc < mc_max
    #         cs = rand(cs_lnd)
    #         if cs < 0.01
    #             cs = 3.2
    #         end
    #         set_param!(GreenDICE,:climatedynamics,:t2xco2,cs) 
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_cs = join(Results_uncertainty_cs,results,on= :time, makeunique = true)
    #         global mc = mc + 1
    #     end
    #     set_param!(GreenDICE,:climatedynamics,:t2xco2,3.2) 
    #     run(GreenDICE)
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_cs.csv",Results_uncertainty_cs)

        global Results_uncertainty_gama3 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
        global mc = 0
        mc_max = 20
        while mc < mc_max              
            tfpi = rand(tfp_param)
            if tfpi < 1
                continue
            end
            elasticity_nc = log(tfpi) / log(K_NC)
            if elasticity_nc < 0
                continue
            end
            set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
            set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc) 
            res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
            best_candidate(res) # optimal vector of miu emissions trajectories
            set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
            set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])       
            run(GreenDICE)
            results = foo()
            global Results_uncertainty_gama3 = join(Results_uncertainty_gama3,results,on= :time, makeunique = true)
            global mc = mc + 1
        end
        CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_gama3.csv",Results_uncertainty_gama3)
        K_NC = 3.878362
        set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
        set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
        Adjusted_tfp = 1.014353
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
        run(GreenDICE)

        
    #     global Results_uncertainty_gama4 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     global mc = 0
    #     global mc = 0
    #     mc_max = 20
    #     while mc < mc_max
    #         g4 = gamma4_ud
    #         set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])       
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_gama4 = join(Results_uncertainty_gama4,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_gama4.csv",Results_uncertainty_gama4)
    #     set_param!(GreenDICE,:green_naturalcapital,:g4,0.8) #gamma4
    #     run(GreenDICE)
            
    #     global Results_uncertainty_share2 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     for sharei in (1:size(share2_i)[1])
    #         Share = share2_i[sharei]
    #         set_param!(GreenDICE,:welfare,:share2,Share)   
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_share2 = join(Results_uncertainty_share2,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_share2.csv",Results_uncertainty_share2)
    #     set_param!(GreenDICE,:welfare,:share2,0.1)
    #     run(GreenDICE)
        


    #     global Results_uncertainty_nc = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     global mc = 0
    #     mc_max = 20
    #     while mc < mc_max
    #         K_NC = rand(k_nc_nd)
    #         set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   # Ratio of NC to K0 
    #         set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)   # Ratio of NC to K0    
    #         elasticity_nc = log(atfp) / log(K_NC)
    #         set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    #         set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc) 
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])   
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_nc = join(Results_uncertainty_nc,results,on= :time, makeunique = true)
    #     end
    #     K_NC = 3.878362
    #     set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    #     set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    #     Adjusted_tfp = 1.014353
    #     elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    #     set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    #     set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)    
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_ratio.csv",Results_uncertainty_nc)
    #     run(GreenDICE)
        


    #     global Results_uncertainty_share = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     for sharei in (1:size(share1_i)[1])
    #         Share = share1_i[sharei]
    #         set_param!(GreenDICE,:welfare,:share,Share)   
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_share = join(Results_uncertainty_share,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_share1.csv",Results_uncertainty_share)
    #     set_param!(GreenDICE,:welfare,:share,0.1)
    #     run(GreenDICE)
        

    #     global Results_uncertainty_theta2 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     for theta2i in (1:size(theta2_i)[1])
    #         Theta2 = theta2_i[theta2i]
    #         set_param!(GreenDICE,:welfare,:theta2,Theta2)   
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_theta2 = join(Results_uncertainty_theta2,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_theta2.csv",Results_uncertainty_theta2)
    #     set_param!(GreenDICE,:welfare,:theta2,mean(theta2_i))
    #     run(GreenDICE)

    #     global Results_uncertainty_theta = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    #     for theta1i in (1:size(theta1_i)[1])
    #         Theta1 = theta1_i[theta1i]
    #         set_param!(GreenDICE,:welfare,:theta,Theta1)   
    #         res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    #         best_candidate(res) # optimal vector of miu emissions trajectories
    #         set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
    #         set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
    #         run(GreenDICE)
    #         results = foo()
    #         global Results_uncertainty_theta = join(Results_uncertainty_theta,results,on= :time, makeunique = true)
    #     end
    #     CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVnonUV_sens_opt_theta.csv",Results_uncertainty_theta)
    #     set_param!(GreenDICE,:welfare,:theta,mean(theta1_i))
    #     run(GreenDICE)


        
    # end
#GreenDICE main specification (end)



##
## All Use values specification (Start)

    set_param!(GreenDICE,:welfare,:utility_path,1)
    set_param!(GreenDICE,:damages,:path,1)        
    set_param!(GreenDICE,:grosseconomy,:greenGDP,1) 
    set_param!(GreenDICE,:welfare,:share,0.1)
    set_param!(GreenDICE,:welfare,:share2,0.)
    set_param!(GreenDICE,:grosseconomy,:share,0.1) 
    set_param!(GreenDICE,:welfare,:sigma_subs,1/(1-mean(theta1_i))) 
    set_param!(GreenDICE,:welfare,:theta,mean(theta1_i)) 
    set_param!(GreenDICE,:welfare,:theta2,1)

    # #Finding the damage parameters (start)
    # run(GreenDICE)
    # Res = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
    # a_4 = best_candidate(Res) #
    # set_param!(GreenDICE,:damages,:a4,a_4[1])
    # set_param!(GreenDICE,:neteconomy,:a4,a_4[1])
    # set_param!(GreenDICE,:damages,:a2,a_k) 
    # #Finding the damage parameters (end)
        run(GreenDICE)
        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])

    run(GreenDICE)
    results = foo()
    CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UV.csv",results)


##
## 
##
## All Use values specification (end)


##
## Market only specification (START)
    set_param!(GreenDICE,:welfare,:utility_path,1)
    set_param!(GreenDICE,:damages,:path,1)        
    set_param!(GreenDICE,:grosseconomy,:greenGDP,1) 
    set_param!(GreenDICE,:welfare,:share,0.)
    set_param!(GreenDICE,:welfare,:share2,0.)
    set_param!(GreenDICE,:grosseconomy,:share,0.) 
    set_param!(GreenDICE,:welfare,:sigma_subs,10^10) 
    set_param!(GreenDICE,:welfare,:theta,1) 
    set_param!(GreenDICE,:welfare,:theta2,1)
    set_param!(GreenDICE,:damages,:a4,a4[1])
    set_param!(GreenDICE,:damages,:a5,0)
    
    # #Finding the damage parameters (start)
    #  run(GreenDICE)
    # Res = bboptimize(damage_params_market;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
    # a_4 = best_candidate(Res) #
    # set_param!(GreenDICE,:damages,:a4,a_4[1])
    # set_param!(GreenDICE,:neteconomy,:a4,a_4[1])
    # set_param!(GreenDICE,:damages,:a2,a_k) 
    # #Finding the damage parameters (end)

        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
    run(GreenDICE)
    results = foo()
    CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/GreenDICE_UVmkt.csv",results)


##
## Market only specification (END)


##DICE

set_param!(GreenDICE,:welfare,:utility_path,0)
set_param!(GreenDICE,:damages,:path,1)        
set_param!(GreenDICE,:grosseconomy,:greenGDP,0) 
set_param!(GreenDICE,:welfare,:share,0.)
set_param!(GreenDICE,:welfare,:share2,0.)
set_param!(GreenDICE,:grosseconomy,:share,0.) 
set_param!(GreenDICE,:welfare,:sigma_subs,10^10) 
set_param!(GreenDICE,:welfare,:theta,1) 
set_param!(GreenDICE,:welfare,:theta2,1)
set_param!(GreenDICE,:damages,:a4,0)
set_param!(GreenDICE,:damages,:a5,0)
a_all = 0.0026686075 #All damages,Nordhaus and Sztorc (2013)
set_param!(GreenDICE,:damages,:a2,a_all) 

    Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
    best_candidate(Res) # optimal vector of miu emissions trajectories
    set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
    set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
run(GreenDICE)
results = foo()
CSV.write("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/DICE.csv",results)

