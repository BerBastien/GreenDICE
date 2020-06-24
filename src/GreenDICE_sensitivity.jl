    global Results_uncertainty_damage = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    
    while mc < sens_max
        D_d = choice(damaged_i)
        D_k = choice(damagek_i)
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
        set_param!(GreenDICE,:damages,:a_d,D_d)
        set_param!(GreenDICE,:damages,:a2,a_k)
        run(GreenDICE)
        Resd = bboptimize(damage_params;SearchRange=(0.,1.), NumDimensions=1, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=3000)
        a_4 = best_candidate(Resd) 
        set_param!(GreenDICE,:damages,:a4,a_4[1])
        set_param!(GreenDICE,:neteconomy,:a4,a_4[1]) 
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_damage = join(Results_uncertainty_damage,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_damage.csv"),Results_uncertainty_damage)
    set_param!(GreenDICE,:damages,:a4,0.0362)
    run(GreenDICE)
    
    global Results_uncertainty_prtp = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    for prtpi in (1:size(prtp_i)[1])
        prtp = prtp_i[prtpi]
        RR = [1/(1+prtp)^(t*5) for t in 0:59]
        set_param!(GreenDICE,:welfare,:rr,RR)   # Ratio of NC to K0 
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_prtp = join(Results_uncertainty_prtp,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_prtp.csv"),Results_uncertainty_prtp)
    RR = [1/(1+0.015)^(t*5) for t in 0:59]
    set_param!(GreenDICE,:welfare,:rr,RR)   # Ratio of NC to K0 
    run(GreenDICE)

    global Results_uncertainty_cs = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    
    while mc < sens_max
        cs = rand(cs_lnd)
        if cs < 0.01
            cs = 3.2
        end
        set_param!(GreenDICE,:climatedynamics,:t2xco2,cs) 
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_cs = join(Results_uncertainty_cs,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    set_param!(GreenDICE,:climatedynamics,:t2xco2,3.2) 
    run(GreenDICE)
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_cs.csv"),Results_uncertainty_cs)

    global Results_uncertainty_gama3 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    while mc < sens_max              
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
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])       
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_gama3 = join(Results_uncertainty_gama3,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_gama3.csv"),Results_uncertainty_gama3)
    K_NC = k_nc_param
    set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    Adjusted_tfp = atfp_param
    elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
    run(GreenDICE)

    
    global Results_uncertainty_gama4 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    
    while mc < sens_max
        g4 = rand( gamma4_ud)
        set_param!(GreenDICE,:green_naturalcapital,:g4,g4) #gamma4
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])       
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_gama4 = join(Results_uncertainty_gama4,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_gama4.csv"),Results_uncertainty_gama4)
    set_param!(GreenDICE,:green_naturalcapital,:g4,0.8) #gamma4
    run(GreenDICE)
        
    global Results_uncertainty_share2 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    for sharei in (1:size(share2_i)[1])
        Share = share2_i[sharei]
        set_param!(GreenDICE,:welfare,:share2,Share)   
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_share2 = join(Results_uncertainty_share2,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_share2.csv"),Results_uncertainty_share2)
    set_param!(GreenDICE,:welfare,:share2,0.1)
    run(GreenDICE)
    


    global Results_uncertainty_nc = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    
    while mc < sens_max
        K_NC = rand(k_nc_nd)
        set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   # Ratio of NC to K0 
        set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)   # Ratio of NC to K0    
        elasticity_nc = log(atfp_param) / log(K_NC)
        set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
        set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc) 
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])   
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_nc = join(Results_uncertainty_nc,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    K_NC = k_nc_param
    set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    Adjusted_tfp = atfp_param
    elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)    
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_ratio.csv"),Results_uncertainty_nc)
    run(GreenDICE)
    


    global Results_uncertainty_share = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    for sharei in (1:size(share1_i)[1])
        Share = share1_i[sharei]
        set_param!(GreenDICE,:welfare,:share,Share)   
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_share = join(Results_uncertainty_share,results,on= :time, makeunique = true)
        global mc = mc + 1
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_share1.csv"),Results_uncertainty_share)
    set_param!(GreenDICE,:welfare,:share,0.1)
    run(GreenDICE)
    

    global Results_uncertainty_theta2 = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    for theta2i in (1:size(theta2_i)[1])
        Theta2 = theta2_i[theta2i]
        set_param!(GreenDICE,:welfare,:theta2,Theta2)   
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_theta2 = join(Results_uncertainty_theta2,results,on= :time, makeunique = true)
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_theta2.csv"),Results_uncertainty_theta2)
    set_param!(GreenDICE,:welfare,:theta2,mean(theta2_i))
    run(GreenDICE)

    global Results_uncertainty_theta = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    for theta1i in (1:size(theta1_i)[1])
        Theta1 = theta1_i[theta1i]
        set_param!(GreenDICE,:welfare,:theta,Theta1)   
        res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
        best_candidate(res) # optimal vector of miu emissions trajectories
        set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
        set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
        run(GreenDICE)
        results = foo()
        global Results_uncertainty_theta = join(Results_uncertainty_theta,results,on= :time, makeunique = true)
    end
    CSV.write(string(dir,"Results/sensitivity/GreenDICE_UVnonUV_sens_opt_theta.csv"),Results_uncertainty_theta)
    set_param!(GreenDICE,:welfare,:theta,mean(theta1_i))
    run(GreenDICE)

