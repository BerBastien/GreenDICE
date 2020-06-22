
    
    global Results_spread_combined_mcs_opt = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    while mc < mc_max


        #Draw 8 parameters (start)
            D_d = choice(damaged_i)
            D_k = choice(damagek_i)
            K_NC = rand(k_nc_nd)
            Share = choice(share1_i)
            Share2 = choice(share2_i)
            gamma4 = rand(gamma4_ud)
            Theta = choice(theta1_i)
            Theta2 = choice(theta2_i)
            atfp = rand(atfp_nd)
        #Draw 8 parameters (end)

        #choose damage parameters (start)
            if D_k == 1
                perck = rand(damagek1_nd)/0.222 #Percentage of damages corresponding to market impacts
            end
            if D_k == 2
                perck = rand(damagek2_ud)
            end
            if D_k == 3
                perck = damagek3
            end
            a_k = perck * D_d
        #choose damage parameters (end)

            if K_NC < 0
                continue
            end
            if atfp < 1
                continue
            end
            elasticity_nc = log(atfp) / log(K_NC)
            if elasticity_nc < 0
                continue #skip to the next step in the for loop
            end
            if elasticity_nc > 0.2999999
                continue #skip to the next step in the for loop
            end
        #conditions (end)

        #set parameters (start)
            set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   
            set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)    
            set_param!(GreenDICE,:welfare,:share,Share)  
            set_param!(GreenDICE,:welfare,:share2,Share2)   
            set_param!(GreenDICE,:green_naturalcapital,:g4,gamma4)   
            set_param!(GreenDICE,:welfare,:theta,Theta)    
            set_param!(GreenDICE,:welfare,:theta2,Theta2)          
            set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)   
            set_param!(GreenDICE,:damages,:a_d,D_d)
            set_param!(GreenDICE,:damages,:a2,a_k)
            run(GreenDICE)
            Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
            a_4 = best_candidate(Resd) 
            set_param!(GreenDICE,:damages,:a4,a_4[1])
            set_param!(GreenDICE,:neteconomy,:a4,a_4[1]) 
        #set parameters (end)
        
        run(GreenDICE)
        #optimize
            res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        #optimize

        #get results (start)
            best_candidate(res) # optimal vector of miu emissions trajectories
            set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
            set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
            run(GreenDICE)
            results = foo()
            global Results_spread_combined_mcs_opt = join(Results_spread_combined_mcs_opt,results,on= :time, makeunique = true)
            global mc = mc + 1
        #get results (end)
    end
    CSV.write(string(dir,"/Results/spread/GreenDICE_UVnonUV_spread_mcs_opt.csv"),Results_spread_combined_mcs_opt)
    #going back to preferred specification (start)
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
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        set_param!(GreenDICE,:damages,:a2,a_k)
        set_param!(GreenDICE,:damages,:a_d,a_d) 
        set_param!(GreenDICE,:damages,:a4,0)
        set_param!(GreenDICE,:neteconomy,:a4,0)
        g4 = 0.8
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)
        Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
        a_4 = best_candidate(Resd) #
        set_param!(GreenDICE,:damages,:a4,a_4[1])
        set_param!(GreenDICE,:neteconomy,:a4,a_4[1])
    #going back to preferred specification (end)
    run(GreenDICE)