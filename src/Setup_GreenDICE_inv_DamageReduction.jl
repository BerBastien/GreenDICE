    include(string(dir,"src/Setup_GreenDICE_mainSpecification.jl"))
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
