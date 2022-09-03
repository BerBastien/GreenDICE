Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=optim_steps)
best_candidate(Res) # optimal vector of miu emissions trajectories
update_param!(GreenDICE,:MIU,best_candidate(Res)[1:60])
update_param!(GreenDICE,:S,best_candidate(Res)[61:120])
run(GreenDICE)
results = foo()
CSV.write(joinpath(dir,"Results/GreenDICE_UVnonUV.csv"),results)