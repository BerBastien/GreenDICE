## GreenDICE main specification (START)
        #set parameters of specification (start)
        
        set_param!(GreenDICE,:welfare,:share,share1_param)
        set_param!(GreenDICE,:welfare,:share2,share2_param)
        set_param!(GreenDICE,:grosseconomy,:share,share1_param)
        set_param!(GreenDICE,:welfare,:theta,theta1_param) 
        set_param!(GreenDICE,:welfare,:theta2,theta2_param)
        set_param!(GreenDICE,:damages,:a5,0)
        K_NC = k_nc_param
        set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
        set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
        Adjusted_tfp = atfp_param
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

        g4 = gamma4_param
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)

        #Finding the damage parameters (start)
            Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
            a_4 = best_candidate(Resd) #
            set_param!(GreenDICE,:damages,:a4,a_4[1])
            set_param!(GreenDICE,:neteconomy,:a4,a_4[1])
        #Finding the damage parameters (end)
    #set parameters of specification (end)


    
    run(GreenDICE)

    global Results_damagereduction = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))

    resinv = bboptimize(eval_dice_inv;SearchRange=(0.,1.0), NumDimensions=180, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
                best_candidate(resinv) # optimal vector of miu emissions trajectories
                set_param!(GreenDICE,:emissions,:MIU,best_candidate(resinv)[1:60])
                set_param!(GreenDICE,:neteconomy,:S,best_candidate(resinv)[61:120])
                set_param!(GreenDICE,:neteconomy,:invNCfrac,best_candidate(resinv)[121:180])
                run(GreenDICE)
                results = foo()
                
            global Results_damagereduction = join(Results_damagereduction,results,on= :time, makeunique = true)
       


CSV.write(string(dir,"Results/investment/ReducedDamages/GreenDICE_ReducedDamages.csv"),Results_damagereduction)
