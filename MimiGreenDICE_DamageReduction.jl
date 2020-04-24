using Mimi
using MimiDICE2013
GreenDICE = MimiDICE2013.get_model()
const model_years = 2010:5:2305 #necesary to compute SCC

@defcomp green_grosseconomy begin
    K       = Variable(index=[time])    #Capital stock (trillions 2005 US dollars)
    YGROSS  = Variable(index=[time])    #Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
    ES  = Variable(index=[time])    #Ecosystem services
    
    NC      = Parameter(index=[time])   #Natural capital; although is a variable, in this component is a parameter
    al      = Parameter(index=[time])   #Level of total factor productivity
    I       = Parameter(index=[time])   #Investment (trillions 2005 USD per year)
    l       = Parameter(index=[time])   #Level of population and labor
    dk      = Parameter()               #Depreciation rate on capital (per year)
    gama    = Parameter()               #Standard DICE gama of K
    gama3   = Parameter()               #GreenDICE: natural Capital elasticity in production function
    k0      = Parameter()               #Initial capital value (trill 2005 USD)
    ratioNC = Parameter()               #GreenDICE: ratio of NC to K0
    greenGDP= Parameter()               #GreenDICE: Is NC in the production function?
    k0         = Parameter()               #Initial capital value (trill 2005 USD)
    ratioNC    = Parameter()               #GreenDICE: ratio of NC to K0
    share           = Parameter()               #GreenDICE: share of ES produced by Green Production Function
    ExtraK = Parameter(index=[time])          #GreenDICE: Year to add an extra asset of K, to compute investments 
    


    function run_timestep(p, v, d, t)
        #Define function for K
        if is_first(t)
            v.K[t] = p.k0 + p.ExtraK[1]
        else
            v.K[t] = (1 - p.dk)^5 * v.K[t-1] + 5 * p.I[t-1]  + p.ExtraK[t]
        end

        #Define function for YGROSS
        if p.greenGDP == 1
            if is_first(t)
                v.YGROSS[t] = (p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama - p.gama3)) * ((p.k0 / p.ratioNC)^p.gama3)
                v.ES[t] = (p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama3)) * ((p.k0 / p.ratioNC)^(p.gama - p.gama3))
            else
                v.YGROSS[t] = ((p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama - p.gama3)) * (p.NC[t-1]^p.gama3))
                v.ES[t] = ((p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama3)) * (p.NC[t-1]^(p.gama - p.gama3)))
                #In this case, the Green Production function yields green GDP, composed by market and non-market goods
            end
                else
                    if is_first(t)
                        v.YGROSS[t] = (p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^p.gama)
                        v.ES[t] = ((p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama3)) * ((p.k0 / p.ratioNC)^(p.gama - p.gama3)))
                     else
                        v.YGROSS[t] = (p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^p.gama)
                        v.ES[t] = ((p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama3)) * (p.NC[t-1]^(p.gama - p.gama3)))
                        #In this case, the standard production function, the output is purely economic market goods
                     end
            end

    end
end

replace_comp!(GreenDICE,green_grosseconomy,:grosseconomy)
set_param!(GreenDICE,:grosseconomy,:k0,135.)
set_param!(GreenDICE,:grosseconomy,:ExtraK,fill(0.,60))


#Redefining damages component to include natural capital and ecosystem services
@defcomp green_damages begin
    DAMAGES    = Variable(index=[time])    #Damages (trillions 2005 USD per year)
    DAMAGES_NC = Variable(index=[time])    # GreenDICE: Natural capital Damages (in a measurement of NC)
    DAMAGES_ES = Variable(index=[time])    # GreenDICE: Ecosystem services Damages (in a measurement of NC)
    DAMFRAC    = Variable(index=[time])    #Damages (fraction of gross output)
    DAMFRAC_NC = Variable(index=[time])    # GreenDICE: Natural Capital Damages (fraction of NC)
    DAMFRAC_ES = Variable(index=[time])    # GreenDICE: Ecosystesm Services Damages (fraction of NC)

    TATM    = Parameter(index=[time])   #Increase temperature of atmosphere (degrees C from 1900)
    YGROSS  = Parameter(index=[time])   #Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
    NC      = Parameter(index=[time])   #GreenDICE: NC
    #ES      = Parameter(index=[time])   #GreenDICE: Ecosystem services
    a1      = Parameter()               #Damage coefficient
    a2      = Parameter()               #Damage quadratic term
    a3      = Parameter()               #Damage exponent
    a4      = Parameter()               #GreenDICE: Damage coefficient for NC
    a5      = Parameter()               #GreenDICE: Damage coefficient for ES
    damadj  = Parameter()               #Adjustment exponent in damage function
    usedamadj::Bool = Parameter()       # Only the Excel version uses the damadj parameter
    path    = Parameter()               #GreenDICE: PAth in which climate damage nature: 0 for none, 1 for NC, 2 for ES
    invNCfrac   = Parameter(index=[time])    #GreenDICE - investment control
    
     
    function run_timestep(p, v, d, t)
        #Define function for DAMFRAC
        v.DAMFRAC[t] = p.a1 * p.TATM[t] + p.a2 * p.TATM[t] ^ p.a3
        #Equation above is as in DICE Matlab 2016
        v.DAMFRAC_NC[t] = (1 / (1 + 0. * p.TATM[t] ^ 2))
        v.DAMFRAC_ES[t] = (1 / (1 + 0. * p.TATM[t] ^ 2))
        if p.path == 1
            if is_first(t)
                v.DAMFRAC_NC[t] = (1 / (1 + p.a4 * p.TATM[t] ^ 2))
            else
                v.DAMFRAC_NC[t] = (1 / (1 + (p.a4/(1+p.invNCfrac[t-1]*1875.60976)) * p.TATM[t] ^ 2)) #271.32 comes from McCarthy, investments in conservation
            end 
           
        elseif p.path == 2
           v.DAMFRAC_ES[t] = (1 / (1 + p.a5 * p.TATM[t] ^ 2))
        end
    

        v.DAMAGES[t] = p.YGROSS[t] * v.DAMFRAC[t]
        v.DAMAGES_NC[t] = v.DAMFRAC_NC[t]
        v.DAMAGES_ES[t] = v.DAMFRAC_ES[t]
    end
end
replace_comp!(GreenDICE,green_damages,:damages)
set_param!(GreenDICE,:damages,:invNCfrac,fill(0.,60))

#Defining the component of natural capital
@defcomp green_naturalcapital begin
    NC       = Variable(index=[time])   #Natural capital
    nonUV       = Variable(index=[time])   #non Use Value
    
    DAMAGES_NC = Parameter(index=[time])   #GreenDICE: total damages to NC 
    #DAMAGES_ES = Parameter(index=[time])   #GreenDICE: total damages to ES 
    k0         = Parameter()               #Initial capital value (trill 2005 USD)
    ratioNC    = Parameter()               #GreenDICE: ratio of NC to K0
    benefitsNC = Parameter(index=[time])   #GreenDICE: fraction of investment in NC
    ExtraN = Parameter(index=[time])          #GreenDICE: Year to add an extra asset of N, to compute investments 
   
   

    function run_timestep(p, v, d, t)
        #Evolution of NC and ES 

        if is_first(t)
            v.NC[t] = p.k0 / p.ratioNC + p.ExtraN[1] #Matching world Bank global estimation of the ratio of NC to K
            #v.ES[t] = 47 #changed from 80 to 47 because this is 2013 version instead of 2016. 47 is taken from original DICE 2013 parameters in the excel file within Mimi
            #v.NC[t] = 47
            v.nonUV[t] = 47 #47 is the initial consumption, following Sterner and Persson
        else
            #I can also include investment, growth, recovery, depreciation, inflation
            v.NC[t] = v.NC[t-1] * p.DAMAGES_NC[t] 
            #v.ES[t] = v.NC[t] * (47 / (p.k0 / p.ratioNC)) * p.DAMAGES_ES[t]
            v.nonUV[t] = v.NC[t] * (47 / (p.k0 / p.ratioNC)) 
            #v.ES[t] = v.NC[t] *  p.DAMAGES_ES[t]
        end

    end
end
#replace_comp!(GreenDICE,green_naturalcapital,:green_naturalcapital)
add_comp!(GreenDICE, green_naturalcapital,:green_naturalcapital, before=:neteconomy)
set_param!(GreenDICE,:green_naturalcapital,:k0,135.)
connect_param!(GreenDICE, :grosseconomy, :NC, :green_naturalcapital, :NC)
connect_param!(GreenDICE, :damages, :NC, :green_naturalcapital, :NC)
#connect_param!(GreenDICE, :damages, :ES, :green_naturalcapital, :ES)
connect_param!(GreenDICE, :green_naturalcapital, :DAMAGES_NC, :damages, :DAMAGES_NC)
#connect_param!(GreenDICE, :green_naturalcapital, :DAMAGES_ES, :damages, :DAMAGES_ES)
set_param!(GreenDICE,:green_naturalcapital,:ExtraN,fill(0.,60))


@defcomp green_welfare begin
    CEMUTOTPER      = Variable(index=[time])    #Period utility
    CUMCEMUTOTPER   = Variable(index=[time])    #Cumulative period utility
    PERIODU         = Variable(index=[time])    #One period utility function
    UTILITY         = Variable()                #Welfare Function
    Consumption_subs   = Variable(index=[time]) #GreenDICE: Internal variable for computing utility
    Environmental_subs = Variable(index=[time]) #GreenDICE: Internal variable for computing utility
    

    CPC             = Parameter(index=[time])   #Per capita consumption (thousands 2005 USD per year)
    l               = Parameter(index=[time])   #Level of population and labor
    rr              = Parameter(index=[time])   #Average utility social discount rate
    elasmu          = Parameter()               #Elasticity of marginal utility of consumption
    scale1          = Parameter()               #Multiplicative scaling coefficient
    scale2          = Parameter()               #Additive scaling coefficient
    utility_path    = Parameter()               #GreenDICE: defines whether ES are in the Utility function
    share           = Parameter()               #GreenDICE: share of utility produced by ES
    share2           = Parameter()               #GreenDICE: share of utility produced by nonUV
    theta      = Parameter()               #GreenDICE: Elasticity of substitution for ES and C
    theta2      = Parameter()               #GreenDICE: Elasticity of substitution for UV and nonUV
    sigma_subs      = Parameter()               #GreenDICE: Substitutability parameter
    #ES              = Parameter(index=[time])               
    ESPC              = Parameter(index=[time])               
    nonUV             = Parameter(index=[time])               



    function run_timestep(p, v, d, t)
        
        #if p.utility_path == 0
            # Define function for PERIODU
        #    v.PERIODU[t] = (p.CPC[t] ^ (1 - p.elasmu) - 1) / (1 - p.elasmu) - 1
        #else 
            #MATLAB: (1 - share)*1000/L(i)*Consumption^(1 - 1/elasticity_subs)
            #v.Consumption_subs[t] = ((1.0 - p.share) * p.CPC[t] ^ (1.0 - 1.0 / p.sigma_subs))
            #v.Environmental_subs[t] = ((p.share) * p.ES[t] ^ (1.0 - 1.0 / p.sigma_subs))
            #U = ((Consumption_subs + Environmental_subs)^((1 - p.elasmu)*(p.sigma_subs/(p.sigma_subs - 1)))-1) / (1 - p.elasmu) - 1;
            #v.PERIODU[t] = ((Consumption_subs[t] + Environmental_subs[t]) ^ ( (1. - p.elasmu) *  (p.sigma_subs / (p.sigma_subs - 1.))) -1.) / (1. - p.elasmu) - 1.
            ### PREVIOUS v.PERIODU[t] = (((1.0 - p.share) * p.CPC[t] ^ (1.0 - 1.0 / p.sigma_subs) + (p.share) * p.ES[t] ^ (1.0 - 1.0 / p.sigma_subs)) ^ ( (1. - p.elasmu) *  (p.sigma_subs / (p.sigma_subs - 1.))) -1.) / (1. - p.elasmu) - 1.          
            #using total consumption: v.PERIODU[t] = ((((1.0 - p.share) * (p.CPC[t] * p.l[t] / 1000) ^ (1.0 - 1.0 / p.sigma_subs)) + ((p.share) * p.ES[t] ^ (1.0 - 1.0 / p.sigma_subs)))^((1 - p.elasmu)*(p.sigma_subs/(p.sigma_subs - 1)))-1) / (1 - p.elasmu) - 1;
            #cpc and ES in trillion US$ v.PERIODU[t] = ((((1.0 - p.share) * (p.CPC[t]) ^ (1.0 - 1.0 / p.sigma_subs)) + ((p.share) * p.ES[t] ^ (1.0 - 1.0 / p.sigma_subs)))^((1 - p.elasmu)*(p.sigma_subs/(p.sigma_subs - 1)))-1) / (1 - p.elasmu) - 1; #edited 15Aug2019 using consumption per capita in thousands of dollars
            #cpc and ES in billion v.PERIODU[t] = ((((1.0 - p.share) * (p.CPC[t]) ^ (1.0 - 1.0 / p.sigma_subs)) + ((p.share) * 1000 * p.ES[t] ^ (1.0 - 1.0 / p.sigma_subs)))^((1 - p.elasmu)*(p.sigma_subs/(p.sigma_subs - 1)))-1) / (1 - p.elasmu) - 1; 
            #v.PERIODU[t] = ((((1.0 - p.share) * (p.CPC[t]) ^ (1.0 - 1.0 / p.sigma_subs)) + ((p.share) * (1000 / p.l[t] * p.ES[t]) ^ (1.0 - 1.0 / p.sigma_subs)))^((1 - p.elasmu)*(p.sigma_subs/(p.sigma_subs - 1)))-1) / (1 - p.elasmu) - 1; #edited 17Aug2019 using ES as thousands of $ per capita
            #v.PERIODU[t] = ((((1.0 - p.share) * (p.CPC[t]) ^ (1.0 - 1.0 / p.sigma_subs)) + ((p.share) * (p.ESPC[t]) ^ (1.0 - 1.0 / p.sigma_subs)))^((1 - p.elasmu)*(p.sigma_subs/(p.sigma_subs - 1)))-1) / (1 - p.elasmu) - 1; #edited 10Sept2019 using ES per capita
            #v.PERIODU[t] = ((((1-p.share2)*((1-p.share)*p.CPC[t]^p.theta)+(p.share*p.ESPC[t]^p.theta))^(p.theta2/p.theta)+p.share2*(p.nonUV[t])^p.theta2)^((1-p.elasmu)/p.theta2)-1)/(1-p.elasmu)-1 #nested utility function with constant elasticities of substiution
            v.PERIODU[t] = ((((1)*((1-p.share)*p.CPC[t]^p.theta)+(p.share*p.ESPC[t]^p.theta))^(p.theta2/p.theta)+p.share2*(p.nonUV[t])^p.theta2)^((1-p.elasmu)/p.theta2)-1)/(1-p.elasmu)-1 #nested utility function with constant elasticities of substiution


        # Define function for CEMUTOTPER
        v.CEMUTOTPER[t] = v.PERIODU[t] * p.l[t] * p.rr[t]

        # Define function for CUMCEMUTOTPER
        v.CUMCEMUTOTPER[t] = v.CEMUTOTPER[t] + (!is_first(t) ? v.CUMCEMUTOTPER[t-1] : 0)

        # Define function for UTILITY
        if t.t == 60
            v.UTILITY = v.CUMCEMUTOTPER[t] 
      end
    end
end
replace_comp!(GreenDICE,green_welfare,:welfare)
connect_param!(GreenDICE, :welfare, :nonUV, :green_naturalcapital, :nonUV)

@defcomp green_neteconomy begin

    ABATECOST   = Variable(index=[time])    #Cost of emissions reductions  (trillions 2005 USD per year)
    C           = Variable(index=[time])    #Consumption (trillions 2005 US dollars per year)
    CPC         = Variable(index=[time])    #Per capita consumption (thousands 2005 USD per year)
    ESPC         = Variable(index=[time])    #Per capita consumption of ES (thousands 2005 USD per year)
    CPRICE      = Variable(index=[time])    #Carbon price (2005$ per ton of CO2)
    I           = Variable(index=[time])    #Investment (trillions 2005 USD per year)
    InvNC       = Variable(index=[time])    #GreenDICE Investment in natural capital (trillions 2005 USD per year)
    benefitsNC  = Variable(index=[time])    #GreenDICE Benefits of the Investment in natural capital 
    MCABATE     = Variable(index=[time])    #Marginal cost of abatement (2005$ per ton CO2)
    Y           = Variable(index=[time])    #Gross world product net of abatement and damages (trillions 2005 USD per year)
    YNET        = Variable(index=[time])    #Output net of damages equation (trillions 2005 USD per year)
    YGreen      = Variable(index=[time])    #GreenDICE: a middle step from gross economy to consumption, to subtract investments in NC

    rr              = Parameter(index=[time])   #Average utility social discount rate
    cost1       = Parameter(index=[time])   #Abatement cost function coefficient
    DAMAGES     = Parameter(index=[time])   #Damages (Trillion $)
    l           = Parameter(index=[time])   #Level of population and labor
    MIU         = Parameter(index=[time])   #Emission control rate GHGs
    partfract   = Parameter(index=[time])   #Fraction of emissions in control regime
    pbacktime   = Parameter(index=[time])   #Backstop price
    S           = Parameter(index=[time])   #Gross savings rate as fraction of gross world product
    YGROSS      = Parameter(index=[time])   #Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
    ES      = Parameter(index=[time])   #Gross world ES
    expcost2    = Parameter()               #Exponent of control cost function
    invNCfrac   = Parameter(index=[time])    #GreenDICE - investment control
    ExtraC      = Parameter(index=[time])    #GreenDICE: Year to add an extra aggregated consumption, to compute SCC
    a4      = Parameter()               #GreenDICE: Damage coefficient for NC
    gama3   = Parameter()               #GreenDICE: natural Capital elasticity in production function
    NC      = Parameter(index=[time])   #GreenDICE: NC
    K      = Parameter(index=[time])   #Capital stock
    priceNC      = Parameter(index=[time])   #exogenous price of NC
    

    function run_timestep(p, v, d, t)
        #Define function for YNET
        v.YNET[t] = p.YGROSS[t] - p.DAMAGES[t]
    
        #Define function for ABATECOST
        v.ABATECOST[t] = p.YGROSS[t] * p.cost1[t] * (p.MIU[t]^p.expcost2) * (p.partfract[t]^(1 - p.expcost2))
    
        #Define function for MCABATE (equation from GAMS version)
        v.MCABATE[t] = p.pbacktime[t] * p.MIU[t]^(p.expcost2 - 1)
    
        #Define function for Y
        v.YGreen[t] = v.YNET[t] - v.ABATECOST[t]
    
        #GreenDICE
        #Define function for Investments in natural capital
        v.InvNC[t] = p.invNCfrac[t] * v.YGreen[t] 

        v.Y[t] = v.YGreen[t] - v.InvNC[t]

        #Define function for I
        v.I[t] = p.S[t] * v.Y[t]
        #v.I[t] = ( p.S[t] / ( 1 - v.InvNC[t])) * v.Y[t]
        
        #GreenDICE
        #Benefits of investment in Natural Capital
        #v.benefitsNC[t] = 1 + ((1-0.01)/(0.01*(0.076/47))) * p.invNCfrac[t] #According to 
        #v.benefitsNC[t] =  v.InvNC[t] / ((p.gama3/(0.3-p.gama3)) * (44760/15841)) #The investment is used to "buy" NC stock at price = ((p.gama3/(0.3-p.gama3)) * (151250/1836)), as given by  dY/dNC/dY/dK 
        #v.benefitsNC[t] =  v.InvNC[t] / ( ((p.gama3/(0.3-p.gama3)) * (p.K[t]/p.NC[t])) * ((0.1+0.015)/(0.015)) )#The investment is used to "buy" NC stock at price = ((p.gama3/(0.3-p.gama3)) * (151250/1836)), as given by  dY/dNC/dY/dK 
        #v.benefitsNC[t] =  v.InvNC[t] / ( ((p.gama3/(0.3-p.gama3)) * (44760/15841)) * ((0.1+0.015)/(0.015)) )#The investment is used to "buy" NC stock at price = ((p.gama3/(0.3-p.gama3)) * (151250/1836)), as given by  dY/dNC/dY/dK         #v.benefitsNC[t] =  v.InvNC[t] / ( ((p.gama3/(0.3-p.gama3)) * (44760/15841)) * (log(0.9/1.015)/(log(1/1.015))) )#The investment is used to "buy" NC stock at price = ((p.gama3/(0.3-p.gama3)) * (151250/1836)), as given by  dY/dNC/dY/dK 
        #v.benefitsNC[t] =  v.InvNC[t] / p.priceNC[t] #using exogenous price
        #((p.rr[2])^(1/5)-1) is the prtp
        #0.1 is the depreciation rate
        v.benefitsNC[t] =  v.InvNC[t] 


        #GreenDICE modification, include investments
        #Define function for C
        v.C[t] = v.Y[t] - v.I[t] + p.ExtraC[t]
    
        #Define function for CPC
        v.CPC[t] = 1000 * v.C[t] / p.l[t]
    
        #Define function for CPRICE (equation from GAMS version of DICE2013)
        v.CPRICE[t] = p.pbacktime[t] * (p.MIU[t] / p.partfract[t])^(p.expcost2 - 1)

        v.ESPC[t] = 1000 * p.ES[t] / p.l[t] #ES per capita
    end
end
replace_comp!(GreenDICE,green_neteconomy,:neteconomy)
connect_param!(GreenDICE, :green_naturalcapital, :benefitsNC, :neteconomy, :benefitsNC)
connect_param!(GreenDICE, :neteconomy, :ES, :grosseconomy, :ES)
connect_param!(GreenDICE, :welfare, :ESPC, :neteconomy, :ESPC)
connect_param!(GreenDICE, :neteconomy, :K, :grosseconomy, :K)
connect_param!(GreenDICE, :neteconomy, :NC, :green_naturalcapital, :NC)
connect_param!(GreenDICE, :damages, :invNCfrac, :neteconomy, :benefitsNC)
set_param!(GreenDICE,:neteconomy,:invNCfrac,fill(0.,60))
set_param!(GreenDICE,:neteconomy,:ExtraC,fill(0.,60))
set_param!(GreenDICE,:neteconomy,:priceNC,fill(10e100,60))
RR = [1/(1+0.015)^(t*5) for t in 0:59]
set_param!(GreenDICE,:neteconomy,:rr,RR)


@defcomp green_emissions begin
    CCA     = Variable(index=[time])    #Cumulative indiustrial emissions
    E       = Variable(index=[time])    #Total CO2 emissions (GtCO2 per year)
    EIND    = Variable(index=[time])    #Industrial emissions (GtCO2 per year)

    etree   = Parameter(index=[time])   #Emissions from deforestation
    MIU     = Parameter(index=[time])   #Emission control rate GHGs
    sigma   = Parameter(index=[time])   #CO2-equivalent-emissions output ratio
    YGROSS  = Parameter(index=[time])   #Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
    cca0    = Parameter()               #Initial cumulative industrial emissions
    ExtraCO2 = Parameter(index=[time])          #GreenDICE: Year to add an extra tonn of CO2, to compute SCC
   

    function run_timestep(p, v, d, t)
        #Define function for EIND
        v.EIND[t] = p.sigma[t] * p.YGROSS[t] * (1-p.MIU[t])
    
        #Define function for E
        
            v.E[t] = v.EIND[t] + p.etree[t] + p.ExtraCO2[t]
        
    
        #Define function for CCA
        if is_first(t)
            v.CCA[t] = p.cca0
        else
            v.CCA[t] = v.CCA[t-1] + v.EIND[t-1] * 5/3.666
        end

    end
end
replace_comp!(GreenDICE,green_emissions,:emissions)
set_param!(GreenDICE,:emissions,:ExtraCO2,fill(0.,60))