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
