## GreenDICE main specification (START)
        #set parameters of specification (start)
        update_param!(GreenDICE,:elasmu,elasmu_param)
        set_param!(GreenDICE,:share,share1_param)
        update_param!(GreenDICE,:welfare,:share2,share2_param)
        set_param!(GreenDICE,:share,share1_param)
        update_param!(GreenDICE,:welfare,:theta,theta1_param) 
        update_param!(GreenDICE,:welfare,:theta2,theta2_param)
        update_param!(GreenDICE,:damages,:a5,0)
        K_NC = k_nc_param
        update_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
        update_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
        Adjusted_tfp = atfp_param
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        update_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        update_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
       
        
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.32  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        
        set_param!(GreenDICE,:a2,a_k)
        update_param!(GreenDICE,:damages,:a_d,a_d) 
        update_param!(GreenDICE,:damages,:a4,0)
        update_param!(GreenDICE,:neteconomy,:a4,0)

        g4 = gamma4_param
        update_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)

        #Finding the damage parameters (start)
            Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
            a_4 = best_candidate(Resd) #
            update_param!(GreenDICE,:damages,:a4,a_4[1])
            update_param!(GreenDICE,:neteconomy,:a4,a_4[1])
        #Finding the damage parameters (end)


            
    #set parameters of specification (end)
