include(string(dir,"src/Setup_GreenDICE_mainSpecification.jl"))
run(GreenDICE)
mc_max = 60
global Results_investment_iterations = getdataframe(GreenDICE,Symbol(d_v[1,1]),Symbol(d_v[1,2]))
global mc = 0
iterations_damagereduction()
CSV.write(string(dir,"Results/investment/AssetInvestment/GreenDICE_AssetInvestment.csv"),Results_investment_iterations)
