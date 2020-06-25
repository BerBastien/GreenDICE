include(string(dir,"src/Setup_GreenDICE_mainSpecification.jl"))
run(GreenDICE)
function iterations_damagereduction()
    while mc < mc_max
                resinv = bboptimize(eval_dice_inv;SearchRange=(0.,1.0), NumDimensions=180, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
                    best_candidate(resinv) # optimal vector of miu emissions trajectories
                    set_param!(GreenDICE,:emissions,:MIU,best_candidate(resinv)[1:60])
                    set_param!(GreenDICE,:neteconomy,:S,best_candidate(resinv)[61:120])
                    set_param!(GreenDICE,:neteconomy,:invNCfrac,best_candidate(resinv)[121:180])
                    run(GreenDICE)
                    results = foo()
                    global Results_investment_iterations = join(Results_investment_iterations,results,on= :time, makeunique = true)
            global mc = mc + 1
        #get results (end)
    end
    return Results_investment_iterations
end
mc_max = 100
global Results_investment_iterations = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
global mc = 0
iterations_damagereduction()
CSV.write(string(dir,"Results/investment/AssetInvestment/GreenDICE_investment.csv"),Results_investment_iterations)
