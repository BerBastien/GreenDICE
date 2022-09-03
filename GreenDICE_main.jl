#Building and running GreenDICE
#set directory
dir = @__DIR__

#max number of steps for optimization
optim_steps = 99 #default = 99999
#number of monte carlo samples
mc_max = 10
#number of samples for the one at a time sensitivity analysis
sens_max = 10

#load required packages:
include(joinpath(dir,"src/Load_req_Packages.jl"))

#build GreenDICE model:
include(joinpath(dir,"src/GreenDICE_model.jl"))

#Build required functions for evaluating GreenDICE model:
include(joinpath(dir,"src/Functions_GreenDICE.jl"))

#Load parameters of the model and its uncertainty ranges/PDFs
include(joinpath(dir,"src/Parameters_Uncertainty.jl"))

#Setup GreenDICE main specification
include(joinpath(dir,"src/Setup_GreenDICE_mainSpecification.jl"))

run(GreenDICE)

#Optimize GreenDICE for welfare maximizing
include(joinpath(dir,"src/Optim_GreenDICE.jl"))

#Run one by one sensitivity analysis:
include(joinpath(dir,"src/GreenDICE_sensitivity.jl"))

#Run combinations of natural capital and total factor productivity
include(joinpath(dir,"src/GreenDICE_NC_TFP.jl"))

#Run Monte Carlo analysis
include(joinpath(dir,"src/MonteCarlo_GreenDICE.jl"))



#Setup and optimize GreenDICE - All use values specification
include(joinpath(dir,"src/Setup_GreenDICE_AllUseValues.jl"))

#Setup and optimize GreenDICE - market-only specification
include(joinpath(dir,"src/Setup_GreenDICE_MarketOnly.jl"))

#Setup and optimize GreenDICE - market-only specification
include(joinpath(dir,"src/Setup_StandardDICE.jl"))

#Run Main GreenDICE specification with investments in damages reduction
include(joinpath(dir,"src/Setup_GreenDICE_inv_DamageReduction.jl"))

#Build special GreenDICE version that includes investments in natural capital assets
include(joinpath(dir,"src/GreenDICE_model_AssetInvestment.jl"))

#Run GreenDICE main specification with investments in natural capital assets
include(joinpath(dir,"src/Setup_GreenDICE_inv_AssetInvestment.jl"))
