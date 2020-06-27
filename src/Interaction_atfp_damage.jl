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
        Adjusted_tfp = 1.05
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
       
        
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        
        set_param!(GreenDICE,:damages,:a2,a_k)
        set_param!(GreenDICE,:damages,:a_d,a_d) 
        set_param!(GreenDICE,:damages,:a4,0.08)
        set_param!(GreenDICE,:neteconomy,:a4,0.08)

        g4 = gamma4_param
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)
        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
        run(GreenDICE)
        results_hatfp_hd = foo()



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
        Adjusted_tfp = 1.05
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
       
        
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        
        set_param!(GreenDICE,:damages,:a2,a_k)
        set_param!(GreenDICE,:damages,:a_d,a_d) 
        set_param!(GreenDICE,:damages,:a4,0.0008)
        set_param!(GreenDICE,:neteconomy,:a4,0.0008)

        g4 = gamma4_param
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)
        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
        run(GreenDICE)
        results_hatfp_ld = foo()


       
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
        Adjusted_tfp = 1.00000000000000000005
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
       
        
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        
        set_param!(GreenDICE,:damages,:a2,a_k)
        set_param!(GreenDICE,:damages,:a_d,a_d) 
        set_param!(GreenDICE,:damages,:a4,0.0008)
        set_param!(GreenDICE,:neteconomy,:a4,0.0008)

        g4 = gamma4_param
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)
        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
        run(GreenDICE)
        results_latfp_ld = foo()


        set_param!(GreenDICE,:welfare,:share,share1_param)
        set_param!(GreenDICE,:welfare,:share2,share2_param)
        set_param!(GreenDICE,:grosseconomy,:share,share1_param)
        set_param!(GreenDICE,:welfare,:theta,theta1_param) 
        set_param!(GreenDICE,:welfare,:theta2,theta2_param)
        set_param!(GreenDICE,:damages,:a5,0)
        K_NC = k_nc_param
        set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
        set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
        Adjusted_tfp = 1.00000000000000000005
        elasticity_nc = log(Adjusted_tfp)/log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
       
        
        a_d = 0.0026686075 #Corresponding to total aggregated damages Norhaus and Sztorc (2013)
        k_perc = 0.31  #Percentage Corresponding to market damages only, Howard and Sterner (2017)
        mortality_perc = 0.2 #Percentage corresponding to mortality damages
        a_k = a_d * (k_perc + mortality_perc)
        
        set_param!(GreenDICE,:damages,:a2,a_k)
        set_param!(GreenDICE,:damages,:a_d,a_d) 
        set_param!(GreenDICE,:damages,:a4,0.08)
        set_param!(GreenDICE,:neteconomy,:a4,0.08)

        g4 = gamma4_param
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        run(GreenDICE)
        Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
        best_candidate(Res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
        run(GreenDICE)
        results_latfp_hd = foo()

        results_latfp_hd.scc[3]
        results_hatfp_hd.scc[3]
        results_hatfp_ld.scc[3]
        results_latfp_ld.scc[3]