#set parameters of specification (start)
include(string(dir,"src/Setup_GreenDICE_mainSpecification.jl"))
set_param!(GreenDICE,:welfare,:share,0)
set_param!(GreenDICE,:welfare,:share2,0)
set_param!(GreenDICE,:grosseconomy,:share,0)
set_param!(GreenDICE,:welfare,:theta,1) 
set_param!(GreenDICE,:welfare,:theta2,1)
set_param!(GreenDICE,:damages,:a5,0)
K_NC = k_nc_param
set_param!(GreenDICE,:grosseconomy,:ratioNC,k_nc_param) 
set_param!(GreenDICE,:green_naturalcapital,:ratioNC,k_nc_param)  
Adjusted_tfp = 1
elasticity_nc = 0
set_param!(GreenDICE,:grosseconomy,:gama3,elasticity_nc)
set_param!(GreenDICE,:neteconomy,:gama3,elasticity_nc)
set_param!(GreenDICE,:damages,:a2,0.0026686075)        


run(GreenDICE)

Res = bboptimize(eval_dice;SearchRange=(0.,1.), NumDimensions=120, Method=:adaptive_de_rand_1_bin_radiuslimited,MaxSteps=99999)
best_candidate(Res) # optimal vector of miu emissions trajectories
set_param!(GreenDICE,:emissions,:MIU,best_candidate(Res)[1:60])
set_param!(GreenDICE,:neteconomy,:S,best_candidate(Res)[61:120])
run(GreenDICE)
results = foo()
CSV.write(string(dir,"Results/GreenDICE_standard.csv"),results)
