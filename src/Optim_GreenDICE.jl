Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
best_candidate(Res) # optimal vector of miu emissions trajectories
set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
run(GreenDICE)
results = foo()
CSV.write(string(dir,"Results/GreenDICE_UVnonUV.csv"),results)

