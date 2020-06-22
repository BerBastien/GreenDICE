
    global Results_combination_NC_tfp = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
    global mc = 0
    while mc < sens_max
            K_NC = rand(k_nc_nd)
            atfp = rand(atfp_nd)

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
                        
            set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC)   # Ratio of NC to K0 
            set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)   # Ratio of NC to K0 
            set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
            set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
            run(GreenDICE)
            res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=49999)
            best_candidate(res)
            set_param!(GreenDICE,:emissions,:MIU,best_candidate(res)[1:60])
            set_param!(GreenDICE,:neteconomy,:S,best_candidate(res)[61:120])
            run(GreenDICE)
            results = foo()
            global Results_combination_NC_tfp = join(Results_combination_NC_tfp,results,on= :time, makeunique = true)
        end
    results = Results_combination_NC_tfp
    #back to standard for future runs
    K_NC = 3.878362
    set_param!(GreenDICE,:grosseconomy,:ratioNC,K_NC) 
    set_param!(GreenDICE,:green_naturalcapital,:ratioNC,K_NC)  
    Adjusted_tfp = 1.014353
    elasticity_nc = log(Adjusted_tfp)/log(K_NC)
    set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
    set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)        
    CSV.write(string(dir,"Results/NCtfp/GreenDICE_UVnonUV_combination_NC_tfp.csv"),results)
    run(GreenDICE)