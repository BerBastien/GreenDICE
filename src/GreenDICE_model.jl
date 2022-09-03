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
    k0         = Parameter()               #Initial capital value (trill 2005 USD)
    ratioNC    = Parameter()               #GreenDICE: ratio of NC to K0
    share           = Parameter()               #GreenDICE: share of ES produced by Green Production Function
    ExtraK = Parameter(index=[time])          #GreenDICE: Year to add an extra asset of K, to compute investments 
    


    function run_timestep(p, v, d, t)
         if is_first(t)
            v.K[t] = p.k0 + p.ExtraK[TimestepIndex(1)] #before: TimestepIndex(1) p.ExtraK[1]
        else
            v.K[t] = (1 - p.dk)^5 * v.K[t-1] + 5 * p.I[t-1]  + p.ExtraK[t]
        end

        #Define function for YGROSS
        if is_first(t)
            v.YGROSS[t] = (p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama - p.gama3)) * ((p.k0 / p.ratioNC)^p.gama3)
            v.ES[t] = (p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama3)) * ((p.k0 / p.ratioNC)^(p.gama - p.gama3))
        else
            v.YGROSS[t] = ((p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama - p.gama3)) * (p.NC[t-1]^p.gama3))
            v.ES[t] = ((p.al[t] * (p.l[t]/1000)^(1-p.gama)) * (v.K[t]^(p.gama3)) * (p.NC[t-1]^(p.gama - p.gama3)))
        end
            

    end
end

#replace_comp!(GreenDICE,green_grosseconomy,:grosseconomy) deprecated
replace!(GreenDICE,:grosseconomy => green_grosseconomy)
#set_param!(GreenDICE,:grosseconomy,:k0,135.)
set_param!(GreenDICE,:k0,135.)
set_param!(GreenDICE,:grosseconomy,:ExtraK,fill(0.,60))


#Redefining damages component to include natural capital and ecosystem services
@defcomp green_damages begin
    DAMAGES    = Variable(index=[time])    #Damages (trillions 2005 USD per year)
    DAMAGES_NC = Variable(index=[time])    # GreenDICE: Natural capital Damages (in a measurement of NC)
    DAMFRAC    = Variable(index=[time])    #Damages (fraction of gross output)
    DAMFRAC_NC = Variable(index=[time])    # GreenDICE: Natural Capital Damages (fraction of NC)
    
    TATM    = Parameter(index=[time])   #Increase temperature of atmosphere (degrees C from 1900)
    YGROSS  = Parameter(index=[time])   #Gross world product GROSS of abatement and damages (trillions 2005 USD per year)
    NC      = Parameter(index=[time])   #GreenDICE: NC
    a1      = Parameter()               #Damage coefficient
    a2      = Parameter()               #Damage quadratic term
    a3      = Parameter()               #Damage exponent
    a4      = Parameter()               #GreenDICE: Damage coefficient for NC
    a5      = Parameter()               #GreenDICE: Damage coefficient for ES
    damadj  = Parameter()               #Adjustment exponent in damage function
    usedamadj=Parameter{Bool}()       # Only the Excel version uses the damadj parameter
    a_d = Parameter()                   #Total aggregated damages (meta paratameter, used to compute damage parameter a4, using damage_params function)
     
    function run_timestep(p, v, d, t)
        v.DAMFRAC[t] = p.a1 * p.TATM[t] + p.a2 * p.TATM[t] ^ p.a3
        
            if is_first(t)
                v.DAMFRAC_NC[t] = (1 / (1 + p.a4 * p.TATM[t] ^ 2))
            else
                v.DAMFRAC_NC[t] = 1 / (1 + (p.a4 * p.TATM[t] ^ 2))
            end 
        v.DAMAGES[t] = p.YGROSS[t] * v.DAMFRAC[t]
        v.DAMAGES_NC[t] = v.DAMFRAC_NC[t]
    end
end
replace!(GreenDICE,:damages=>green_damages)

#Defining the component of natural capital
@defcomp green_naturalcapital begin
    NC       = Variable(index=[time])   #Natural capital
    nonUV       = Variable(index=[time])   #non Use Value
    
    DAMAGES_NC = Parameter(index=[time])   #GreenDICE: total damages to NC 
    k0         = Parameter()               #Initial capital value (trill 2005 USD)
    ratioNC    = Parameter()               #GreenDICE: ratio of NC to K0
    benefitsNC = Parameter(index=[time])   #GreenDICE: fraction of investment in NC
    ExtraN = Parameter(index=[time])          #GreenDICE: Year to add an extra asset of N, to compute investments 
    g4 = Parameter()                        #gamma4, elasticity of natural capital
    invNCfrac   = Parameter(index=[time])    #GreenDICE - investment control
    w = Parameter()                  #damage reduction parameter
    
    function run_timestep(p, v, d, t)
    
        if is_first(t)
            v.NC[t] = p.k0 / p.ratioNC + p.ExtraN[TimestepIndex(1)]
            v.nonUV[t] = v.NC[t] ^ p.g4 
        else
            v.NC[t] = v.NC[t-1] - (v.NC[t-1] - v.NC[t-1] * p.DAMAGES_NC[t])* (1 - (p.benefitsNC[t-1])^(1/p.w))
            v.nonUV[t] = v.NC[t] ^ p.g4
        end
    end
end
add_comp!(GreenDICE, green_naturalcapital,:green_naturalcapital, before=:neteconomy)
update_param!(GreenDICE,:k0,135.)
connect_param!(GreenDICE, :grosseconomy, :NC, :green_naturalcapital, :NC)
connect_param!(GreenDICE, :damages, :NC, :green_naturalcapital, :NC)
connect_param!(GreenDICE, :green_naturalcapital, :DAMAGES_NC, :damages, :DAMAGES_NC)
set_param!(GreenDICE,:green_naturalcapital,:ExtraN,fill(0.,60))
set_param!(GreenDICE,:green_naturalcapital,:invNCfrac,fill(0.,60))
set_param!(GreenDICE,:green_naturalcapital,:w,2.8)


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
    share           = Parameter()               #GreenDICE: share of utility produced by ES
    share2           = Parameter()               #GreenDICE: share of utility produced by nonUV
    theta      = Parameter()               #GreenDICE: Elasticity of substitution for ES and C
    theta2      = Parameter()               #GreenDICE: Elasticity of substitution for UV and nonUV
    ESPC              = Parameter(index=[time])               
    nonUV             = Parameter(index=[time])
    ES              = Parameter(index=[time])                              



    function run_timestep(p, v, d, t)
        
        v.PERIODU[t] = ((((1)*((1-p.share)*p.CPC[t]^p.theta)+(p.share*p.ESPC[t]^p.theta))^(p.theta2/p.theta)+p.share2*(p.nonUV[t])^p.theta2)^((1-p.elasmu)/p.theta2)-1)/(1-p.elasmu)-1 #nested utility function with constant elasticities of substiution


        v.CEMUTOTPER[t] = v.PERIODU[t] * p.l[t] * p.rr[t]

        v.CUMCEMUTOTPER[t] = v.CEMUTOTPER[t] + (!is_first(t) ? v.CUMCEMUTOTPER[t-1] : 0)

        if t.t == 60
            v.UTILITY = v.CUMCEMUTOTPER[t] 
      end
    end
end
replace!(GreenDICE,:welfare=>green_welfare)
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
    NCcost     = Variable(index=[time])    #GreenDICE: a middle step from gross economy to consumption, to subtract investments in NC

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
    w      = Parameter()   #exogenous cost of abatement
    

    function run_timestep(p, v, d, t)
        #Define function for YNET
        v.YNET[t] = p.YGROSS[t] - p.DAMAGES[t]
    
        #Define function for ABATECOST
        v.ABATECOST[t] = p.YGROSS[t] * p.cost1[t] * (p.MIU[t]^p.expcost2) * (p.partfract[t]^(1 - p.expcost2))
    
        #Define function for MCABATE (equation from GAMS version)
        v.MCABATE[t] = p.pbacktime[t] * p.MIU[t]^(p.expcost2 - 1)
    
        #Define function for Y
        v.YGreen[t] = v.YNET[t] - v.ABATECOST[t]
    
        #Define function for Investments in natural capital
        v.InvNC[t] = v.YGreen[t] * p.invNCfrac[t]

        v.Y[t] = v.YGreen[t] - v.InvNC[t]

        #Define function for I
        v.I[t] = p.S[t] * v.Y[t]
        #v.I[t] = ( p.S[t] / ( 1 - v.InvNC[t])) * v.Y[t]
        
        v.benefitsNC[t] =  p.invNCfrac[t]


        #Define function for C
        v.C[t] = v.Y[t] - v.I[t] + p.ExtraC[t]
    
        #Define function for CPC
        v.CPC[t] = 1000 * v.C[t] / p.l[t]
    
        #Define function for CPRICE (equation from GAMS version of DICE2013)
        v.CPRICE[t] = p.pbacktime[t] * (p.MIU[t] / p.partfract[t])^(p.expcost2 - 1)

        v.ESPC[t] = 1000 * p.ES[t] / p.l[t] #ES per capita
    end
end
replace!(GreenDICE,:neteconomy=>green_neteconomy)
connect_param!(GreenDICE, :green_naturalcapital, :benefitsNC, :neteconomy, :benefitsNC)
connect_param!(GreenDICE, :neteconomy, :ES, :grosseconomy, :ES)
connect_param!(GreenDICE, :welfare, :ES, :grosseconomy, :ES)
connect_param!(GreenDICE, :welfare, :ESPC, :neteconomy, :ESPC)
connect_param!(GreenDICE, :neteconomy, :K, :grosseconomy, :K)
connect_param!(GreenDICE, :neteconomy, :NC, :green_naturalcapital, :NC)
connect_param!(GreenDICE, :green_naturalcapital, :invNCfrac, :neteconomy, :benefitsNC)
set_param!(GreenDICE,:invNCfrac,fill(0.,60))
set_param!(GreenDICE,:neteconomy,:ExtraC,fill(0.,60))
set_param!(GreenDICE,:neteconomy,:priceNC,fill(10e100,60))
RR = [1/(1+0.015)^(t*5) for t in 0:59]
set_param!(GreenDICE,:rr,RR)
set_param!(GreenDICE,:w,2.8)


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
replace!(GreenDICE,:emissions=>green_emissions)
set_param!(GreenDICE,:emissions,:ExtraCO2,fill(0.,60))
set_param!(GreenDICE,:k0,135.)