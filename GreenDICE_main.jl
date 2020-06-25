#Building and running GreenDICE
#set directory
dir = "C:/Users/bastien/Documents/GitHub/GreenDICE/"

#max number of steps for optimization
optim_steps = 99999 #default = 99999
#number of monte carlo samples
mc_max = 10
#number of samples for the one at a time sensitivity analysis
sens_max = 10

#load required packages:
include(string(dir,"src/Load_req_Packages.jl"))

#build GreenDICE model:
include(string(dir,"src/GreenDICE_model.jl"))

#Build required functions for evaluating GreenDICE model:
include(string(dir,"src/Functions_GreenDICE.jl"))

#Load parameters of the model and its uncertainty ranges/PDFs
include(string(dir,"src/Parameters_Uncertainty.jl"))

#Setup GreenDICE main specification
include(string(dir,"src/Setup_GreenDICE_mainSpecification.jl"))

run(GreenDICE)

#Optimize GreenDICE for welfare maximizing
include(string(dir,"src/Optim_GreenDICE.jl"))

#Run one by one sensitivity analysis:
include(string(dir,"src/GreenDICE_sensitivity.jl"))

#Run combinations of natural capital and total factor productivity
include(string(dir,"src/GreenDICE_NC_TFP.jl"))

#Run Monte Carlo analysis
include(string(dir,"src/MonteCarlo_GreenDICE.jl"))



#Setup and optimize GreenDICE - All use values specification
include(string(dir,"src/Setup_GreenDICE_AllUseValues.jl"))

#Setup and optimize GreenDICE - market-only specification
include(string(dir,"src/Setup_GreenDICE_MarketOnly.jl"))

#Run Main GreenDICE specification with investments in damages reduction
include(string(dir,"src/Setup_GreenDICE_inv_DamageReduction.jl"))

#Build special GreenDICE version that includes investments in natural capital assets
include(string(dir,"src/GreenDICE_model_AssetInvestment.jl"))

#Run GreenDICE main specification with investments in natural capital assets
include(string(dir,"src/Setup_GreenDICE_inv_AssetInvestment.jl"))
