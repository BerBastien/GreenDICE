#####
#####
##### Packages (start)
#####
#####
    install.packages("tidyverse")
    install.packages("lubridate")
    install.packages("ggpubr")
    install.packages("ggridges")
    install.packages("wesanderson")
    install.packages("grid")
    install.packages("ggExtra")
    install.packages("gridExtra")
    install.packages("ggplot2")
    install.packages("ggalt")
    install.packages("randomForest")
    install.packages("rpart")
    library("ggExtra")
    library("gridExtra")
    library(grid)
    library(wesanderson)
    library(rpart)
    library("ggridges")
    library("ggplot2")
    library(lubridate)
    library("ggpubr")
    library("tidyverse")
    library(plyr)
    library(randomForest)
    library(randomForestExplainer)
    cbp1 <- c( "#0072B2","#009E73","#E69F00")
    mypath = "C:/Users/bastien/Documents/GitHub/GreenDICE/Results"
#####
#####
##### Packages (end)
#####
#####

  
#####
#####
##### Data preparation (start)
#####
#####





    #####
    #####
    ##### READING FILES (START)
    #####
    #####
        #Reading simple runs of GeenDICE (start)
            setwd("C:/Users/bastien/Documents/GitHub/GreenDICE/Results")
            files <- list.files(pattern = "\\.csv$")
            for (i in 1:length(files)){
            results_i <- read.csv(file=files[i], header=TRUE, sep=",")
            exp_code = substr(files[i],11,nchar(files[i])-4) #getting the name of the experiment
            names = colnames(results_i)
            colnames(results_i) = paste(names,exp_code, sep = "_")
            if (i ==1) {
            results = results_i
            } else {
                results = cbind(results,results_i)
            }
            }
            num_exp = length(files)
        #Reading simple runs of GeenDICE (end)

        #Reading spread of GeenDICE (start)

            setwd("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/spread")
            files <- list.files(pattern = "\\.csv$")
            Results_i <- read.csv(file=files[1], header=TRUE, sep=",")
            exp_code = substr(files[1],13,nchar(files[1])-4) #getting the name of the experiment
            names = colnames(Results_i)
            colnames(Results_i) = paste(names,exp_code, sep = "_")
            Results = Results_i

        #Reading spread of GeenDICE (end)


        #Reading NC_gama3 (start)
            setwd("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/NCtfp")
            files <- list.files(pattern = "\\.csv$")
            for (i in 1:length(files)){
            Results_nc_i <- read.csv(file=files[i], header=TRUE, sep=",")
            exp_code = substr(files[i],11,nchar(files[i])-4) #getting the name of the experiment
            names = colnames(Results_nc_i)
            colnames(Results_nc_i) = paste(names,exp_code, sep = "_")
            if (i ==1) {
            Results_nc = Results_nc_i
            num_cols = length(names)
            sens = substr(names[1],19,nchar(names))
            } else {
                Results_nc = cbind(Results_nc,Results_nc_i)
                num_cols = cbind(num_cols, length(names))
                sens = cbind(sens,substr(names[1],19,nchar(names)))
            }
            }
            num_exp_s = length(files)
        #Reading NC_gama3 (end)

        
            #Reading investment (start)
                setwd("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/investment/AssetInvestment")
                files <- list.files(pattern = "\\.csv$")
                for (i in 1:length(files)){
                Results_inv1_i <- read.csv(file=files[i], header=TRUE, sep=",")
                exp_code = substr(files[i],11,nchar(files[i])-4) #getting the name of the experiment
                names = colnames(Results_inv1_i)
                colnames(Results_inv1_i) = paste(names,exp_code, sep = "_")
                if (i ==1) {
                    Results_inv1 = Results_inv1_i
                } else {
                    Results_inv1 = cbind(Results_inv1,Results_inv1_i)
                }
                num_cols = length(names)
                sens = substr(names[1],19,nchar(names))
                }
                Results_AssetInv = Results_inv1


                setwd("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/investment/ReducedDamages")
                files <- list.files(pattern = "\\.csv$")
                for (i in 1:length(files)){
                Results_inv_i <- read.csv(file=files[i], header=TRUE, sep=",")
                exp_code = substr(files[i],11,nchar(files[i])-4) #getting the name of the experiment
                names = colnames(Results_inv_i)
                colnames(Results_inv_i) = paste(names,exp_code, sep = "_")
                Results_inv = Results_inv_i
                num_cols = length(names)
                sens = substr(names[1],19,nchar(names))
                }
                Results_ReducedD = Results_inv
            #Reading investment (end)

 
        
        # Read parametric sensitivity GreenDICE (START)
            setwd("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/sensitivity")
            files <- list.files(pattern = "\\.csv$")
            for (i in 1:length(files)){
            Results_s_i <- read.csv(file=files[i], header=TRUE, sep=",")
            exp_code = substr(files[i],11,nchar(files[i])-4) #getting the name of the experiment
            names = colnames(Results_s_i)
            colnames(Results_s_i) = paste(names,exp_code, sep = "_")
            if (i ==1) {
            Results_s = Results_s_i
            num_cols = length(names)
            sens = substr(names[1],19,nchar(names))
            } else {
                Results_s = cbind(Results_s,Results_s_i)
            num_cols = cbind(num_cols, length(names))
            sens = cbind(sens,substr(names[1],19,nchar(names)))
            }
            } 
            num_exp_s = length(files)
        #Read parametric sensitivity UVnonUV (END)
    #####
    #####
    ##### READING FILES (END)
    #####
    #####

    #####
    #####
    ##### ARRANGE DATASETS (START)
    #####
    #####



        #arrange plain runs of specifications (start)
            results[1]
            num_vars = 26
            num_exp = 4
            r_t  <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+2)]) %>%
            gather(key = "variable", value = "value_t", -1)
            names(r_t)[1] = "years"

            r_e  <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+3)]) %>%
            gather(key = "variable", value = "value_e", -1)
            names(r_e)[1] = "years"

            r_scc  <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+14)]) %>%
            gather(key = "variable", value = "value_scc", -1)
            names(r_scc)[1] = "years"

            r_nc  <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+5)]) %>%
            gather(key = "variable", value = "value_nc", -1)
            names(r_scc)[1] = "years"
            r_inv <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+9)]) %>%
            gather(key = "variable", value = "value_inv", -1)
            
            r_YGreen <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+25)]) %>%
            gather(key = "variable_YGreen", value = "value_YGreen", -1)

            r_YGross <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+7)]) %>%
            gather(key = "variable_YGross", value = "value_YGross", -1)
            
            r_NC <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+4)]) %>%
            gather(key = "variable_nc", value = "value_nc", -1)

            r_k <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+13)]) %>%
            gather(key = "variable_k", value = "value_k", -1)


            r_S <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+8)]) %>%
            gather(key = "variable_s", value = "value_s", -1)

            r_inv <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+9)]) %>%
            gather(key = "variable", value = "value_inv", -1)

            r_miu <- results %>%
            select(names(results)[c(1,num_vars*(0:(num_exp-1))+10)]) %>%
            gather(key = "variable_miu", value = "value_miu", -1)



            
            r_nc_perc <- r_NC[3] / (r_NC[3] + r_k[3])

            

            df_r <- cbind(r_t,r_e[3],r_scc[3]*1.18,r_YGreen[3],
            r_NC[3],r_S[3],r_inv[3],r_miu[3],r_k[3],r_YGross[3],r_nc_perc) #multiplied by 1.18 to pass from 2010 usd to 2019 usd
            names(df_r)[1] = "years"
            names(df_r)[dim(df_r)[2]] = "r_nc_perc"

        #arrange plain runs of specifications (end)


        ###Arrange sensitivity data (START)
            sens = c("cs","damage","elasmu","gama3","gama4","prtp","ratio","share1","share2","theta1","theta2")
            num_exp_ss = (num_cols-2)/25
            #num_exp_ss = c(10,10,10,10,3,10,3,3,16,4)
            for (num_exp_s_i in 1:11){
            initial = sum(num_cols[0:(num_exp_s_i-1)])
            final = sum(num_cols[0:num_exp_s_i])
            Results_s_i = Results_s[,(initial+1):final]
            num_vars = 25
            names(Results_s_i)
            #num_exp = (num_cols[num_exp_s_i]-2)/num_vars
            num_exp = num_exp_ss[num_exp_s_i]
            df_t <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+3)]) %>%
            gather(key = "variable", value = "value_t", -1)

            df_e <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+4)]) %>%
            gather(key = "variable_e", value = "value_e", -1)

            df_scc <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+15)]) %>%
            gather(key = "variable_scc", value = "value_scc", -1)
            df_scc[3] <- df_scc[3]*1.18  #multiplied by 1.18 to pass from 2010usd to 2019 usd


            df_gama <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+16)]) %>%
            gather(key = "variable_gama3", value = "value_gama3", -1)

            df_ratio <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+17)]) %>%
            gather(key = "variable_ratio", value = "value_ratio", -1)

            df_gama4 <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+18)]) %>%
            gather(key = "variable_gama4", value = "value_gama4", -1)




            df_theta1 <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+19)]) %>%
            gather(key = "variable_theta1", value = "value_theta1", -1)

            df_share <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+20)]) %>%
            gather(key = "variable_share1", value = "value_share1", -1)


            df_cs <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+21)]) %>%
            gather(key = "variable_cs", value = "value_cs", -1)

            df_prtp <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+22)]) %>%  
            gather(key = "variable_ratio", value = "value_prtp", -1)

            df_share2 <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+23)]) %>%
            gather(key = "variable_share2", value = "value_share2", -1)

            df_theta2 <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+24)]) %>%
            gather(key = "variable_theta2", value = "value_theta2", -1)


            df_damage <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+27)]) %>%
            gather(key = "variable_damage", value = "value_damage", -1)
            
            df_elasmu <- Results_s_i %>%
            select(names(Results_s_i)[c(1,num_vars*(0:(num_exp-1))+26)]) %>%
            gather(key = "variable_elasmu", value = "value_elasmu", -1)


            
            df_atfp <- df_ratio[3]^df_gama[3]
            names(df_atfp)= "value_gama3"

            id = rep((1+(num_exp_s_i-1)*num_exp):((num_exp_s_i)*num_exp), each = dim(df_t)[1]/num_exp)
            df_id = data.frame("id" = id)

            df <- cbind(df_t,df_e[3],df_scc[3],df_atfp,df_ratio[3],df_share[3],df_theta1[3],df_prtp[3],df_cs[3],df_share2[3],df_theta2[3],df_damage[3],df_gama4[3],df_elasmu[3],df_id)
            df$sensitivity = sens[num_exp_s_i]

            names(df)[1] = "years"
            if (num_exp_s_i  == 1) {
            DF = df } else {
            DF = rbind(DF,df)
            }
            }
        ####Arrange sensitivity data (END)


        # arrange combination NC_TFP of UVnonUV montecarlo simulation optimized (start)
            num_vars = 25
            num_exp = 50

            df_t <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+3)]) %>%
            gather(key = "variable", value = "value_t", -1)

            df_e <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+4)]) %>%
            gather(key = "variable_e", value = "value_e", -1)

            df_scc <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+15)]) %>%
            gather(key = "variable_scc", value = "value_scc", -1)
            df_scc[3] <- df_scc[3]*1.18  #multiplied by 1.18 to pass from 2010usd to 2019 usd



            df_gama <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+16)]) %>%
            gather(key = "variable_gama3", value = "value_gama3", -1)

            df_ratio <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+17)]) %>%
            gather(key = "variable_ratio", value = "value_ratio", -1)

            
            df_atfp <- df_ratio[3]^df_gama[3]
            names(df_atfp)[1] <- "value_atfp"

            df_subs <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+19)]) %>%
            gather(key = "variable_ratio", value = "value_subs", -1)

            df_share <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+20)]) %>%
            gather(key = "variable_ratio", value = "value_share", -1)


            df_cs <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+21)]) %>%
            gather(key = "variable_ratio", value = "value_cs", -1)

            df_prtp <- Results_nc %>%
            select(names(Results_nc)[c(1,num_vars*(0:(num_exp-1))+22)]) %>%
            gather(key = "variable_ratio", value = "value_prtp", -1)

            df <- cbind(df_t,df_e[3],df_scc[3],df_gama[3],df_ratio[3],df_atfp,df_share[3],df_subs[3],df_prtp[3],df_cs[3])
            df_nc <- df
        # arrange  combination NC_TFP of UVnonUV montecarlo simulation optimized (end)

        # arrange spread of UVnonUV montecarlo simulation optimized (start)
            num_vars = 25
            num_exp = (dim(Results)[2]-2)/(num_vars)
            df_t <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+3)]) %>%
            gather(key = "variable", value = "value_t", -1)

            df_e <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+4)]) %>%
            gather(key = "variable_e", value = "value_e", -1)

            df_cpc <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+12)]) %>%
            gather(key = "variable_cpc", value = "value_cpc", -1)

            df_scc <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+15)]) %>%
            gather(key = "variable_scc", value = "value_scc", -1)
            df_scc[3] <- df_scc[3]*1.18  #multiplied by 1.18 to pass from 2010usd to 2019 usd


            df_gama <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+16)]) %>%
            gather(key = "variable_gama3", value = "value_gama3", -1)

            df_gama4 <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+18)]) %>%
            gather(key = "variable_gama4", value = "value_gama4", -1)

            df_ratio <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+17)]) %>%
            gather(key = "variable_ratio", value = "value_ratio", -1)

            df_atfp <- df_ratio[3]^df_gama[3]
            names(df_atfp)[1] <- "value_atfp"

            df_theta1 <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+19)]) %>%
            gather(key = "variable_theta1", value = "value_theta1", -1)

            df_share1 <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+20)]) %>%
            gather(key = "variable_share1", value = "value_share1", -1)

            df_cs <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+21)]) %>%
            gather(key = "variable_cs", value = "value_cs", -1)

            df_prtp <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+22)]) %>%
            gather(key = "variable_prtp", value = "value_prtp", -1)

            df_share2 <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+23)]) %>%
            gather(key = "variable_ratio", value = "value_share2", -1)

            df_theta2 <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+24)]) %>%
            gather(key = "variable_ratio", value = "value_theta2", -1)


            df_damage <- Results %>%
            select(names(Results)[c(1,num_vars*(0:(num_exp-1))+27)]) %>%
            gather(key = "variable_damage", value = "value_damage", -1)

            df <- cbind(df_t,df_e[3],df_scc[3],df_atfp,df_ratio[3],df_damage[3],df_share1[3],df_share2[3],df_theta2[3],df_theta1[3],df_gama4[3],df_prtp[3],df_cs[3],df_cpc[3])
            df_mcs <- df
            glimpse(df_mcs)
        # arrange spread of UVnonUV montecarlo simulation optimized (end)


        # arrange iterations of GreenDICE + Asset investment (start)
            num_vars = 25
            num_expinv = (dim(Results_AssetInv)[2]-2)/(num_vars)
            num_exp = num_expinv


            df_t_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+3)]) %>%
            gather(key = "variable", value = "value_t", -1)
            
            df_e_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+4)]) %>%
            gather(key = "variable_e", value = "value_e", -1)
            df_scc_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+15)]) %>%
            gather(key = "variable_scc", value = "value_scc", -1)
            df_scc_inv[3] <- df_scc_inv[3]*1.18  #multiplied by 1.18 to pass from 2010usd to 2019 usd
            df_inv_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+10)]) %>%
            gather(key = "variable_inv", value = "value_inv", -1)
            
            df_YGreen <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+25)]) %>%
            gather(key = "variable_YGreen", value = "value_YGreen", -1)
            
            df_Y_inv = cbind(df_inv_inv[,1:2],df_inv_inv[3] * df_YGreen[3] *1.18*10^3)  #multiplied by 1.18 to pass from 2010usd to 2019 usd, multiplied by 10^3 to pass from trillion to billion
            names(df_Y_inv) = c("time_UVnonUV_inv","variable","value_Y")
            
            #df_price_inv <- Results_AssetInv %>%
            #select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+29)]) %>%
            #gather(key = "variable_price", value = "value_price", -1)

            df_inv_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+28)]) %>%
            gather(key = "variable_price", value = "value_price", -1)


            df_NC_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+5)]) %>%
            gather(key = "variable_nc", value = "value_nc", -1)


            
            df_YGross_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+8)]) %>%
            gather(key = "variable_YGross", value = "value_YGross", -1)

            df_k_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+14)]) %>%
            gather(key = "variable_k", value = "value_k", -1)


            df_S_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+9)]) %>%
            gather(key = "variable_s", value = "value_s", -1)

            df_inv_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+10)]) %>%
            gather(key = "variable_inv", value = "value_inv", -1)

            df_miu_inv <- Results_AssetInv %>%
            select(names(Results_AssetInv)[c(1,num_vars*(0:(num_exp-1))+10)]) %>%
            gather(key = "variable_miu", value = "value_miu", -1)

            id_var = rep(1:num_exp, each=60)

        
            df_nc_perc = df_NC_inv[3] / (df_NC_inv[3] + df_k_inv[3])
            df_inv <- cbind(df_t_inv,df_e_inv[3],df_scc_inv[3],df_YGreen[3],df_Y_inv[3],
                            #df_price_inv[3]*1.18,df_NC_inv[3],df_S_inv[3],df_inv_inv[3],df_miu_inv[3],df_k_inv[3],df_YGross_inv[3],df_nc_perc,id_var)
                            df_NC_inv[3],df_S_inv[3],df_inv_inv[3],df_miu_inv[3],df_k_inv[3],df_YGross_inv[3],df_nc_perc,id_var)

            names(df_inv)[1]="years"
            names(df_inv)[dim(df_inv)[2]-1] = "nc_perc"
            dim(df_inv)
            df_inv_last<- cbind(df_inv[1:60,2],aggregate(df_inv[(dim(df_inv)[1]-60*10+1):dim(df_inv)[1],3:15], by = list(df_inv[(dim(df_inv)[1]-60*10+1):dim(df_inv)[1],1]), mean))
            names(df_inv_last)[1:2] <- c("variable","years")
            
        # arrange iterations of GreenDICE + Asset investment (end)

        #arrange plain results with ReducedDamages (start)
            num_vars=25
            Results_ReducedDamages <- Results_ReducedD[,-2]
            num_exp = 3
            df_t_inv <- Results_ReducedDamages %>%
                select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+2)]) %>%
                gather(key = "variable", value = "value_t", -1)
            df_e_inv <- Results_ReducedDamages %>%
                select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+3)]) %>%
                gather(key = "variable_e", value = "value_e", -1)
            df_scc_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+14)]) %>%
            gather(key = "variable_scc", value = "value_scc", -1)
            df_scc_inv[3] <- df_scc_inv[3]*1.18  #multiplied by 1.18 to pass from 2010usd to 2019 usd
            df_inv_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+9)]) %>%
            gather(key = "variable_inv", value = "value_inv", -1)
            
            df_YGreen <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+24)]) %>%
            gather(key = "variable_YGreen", value = "value_YGreen", -1)
            
            df_Y_inv = cbind(df_inv_inv[,1:2],df_inv_inv[3] * df_YGreen[3] *1.18*10^3)  #multiplied by 1.18 to pass from 2010usd to 2019 usd, multiplied by 10^3 to pass from trillion to billion
            names(df_Y_inv) = c("time_UVnonUV_inv","variable","value_Y")
            

            df_NC_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+4)]) %>%
            gather(key = "variable_nc", value = "value_nc", -1)

                df_YGross <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+7)]) %>%
            gather(key = "variable_YGross", value = "value_YGross", -1)


            df_k_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+13)]) %>%
            gather(key = "variable_k", value = "value_k", -1)


            df_S_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+8)]) %>%
            gather(key = "variable_s", value = "value_s", -1)

            df_inv_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+9)]) %>%
            gather(key = "variable_inv", value = "value_inv", -1)

            df_miu_inv <- Results_ReducedDamages %>%
            select(names(Results_ReducedDamages)[c(1,num_vars*(0:(num_exp-1))+10)]) %>%
            gather(key = "variable_miu", value = "value_miu", -1)

            id_var = rep(1:num_exp, each=60)

            df_nc_perc = df_NC_inv[3] / (df_NC_inv[3] + df_k_inv[3])
            
            df_inv_reduceddamages <- cbind(df_t_inv,df_e_inv[3],df_scc_inv[3],df_YGreen[3],df_Y_inv[3],
                df_NC_inv[3],df_S_inv[3],df_inv_inv[3],df_miu_inv[3],df_k_inv[3],df_YGross[3],df_nc_perc,id_var)

                names(df_inv_reduceddamages)[1]="years"
                names(df_inv_reduceddamages)[dim(df_inv_reduceddamages)[2]-1]="nc_perc"

            #df_inv_reduceddamages <- cbind(df_inv_reduceddamages[1:60,2],aggregate(df_inv_reduceddamages[,3:14], list(df_inv_reduceddamages$years), mean))
            #names(df_inv_reduceddamages)[1:2] <- c("variable","years")
            df_inv_reduceddamages_sens <- df_inv_reduceddamages
            df_inv_reduceddamages <- df_inv_reduceddamages_sens[61:120,]
        #arrange plain results with ReducedDamages (end)


    #####
    #####
    ##### ARRANGE DATASETS (END)
    #####
    #####


    #####
    #####
    ##### GETTING QUANTILES (START)
    #####
        ##### quantiles spread
            df_t <- df_t[df_prtp$value==0.015,]
            #df_t <- df_t[df_damage$value==0.00480515,]
            qs_t = data.frame(
            do.call(
                rbind,
                tapply(
                    df_t$value, df_t[1] , function(i){quantile(i)})),
            t=1:60)
            qs_t$t = df_t[1:60,1]


            df_e <- df_e[df_prtp$value==0.015,]
            #df_e <- df_e[df_damage$value==0.00480515,]
            qs_e = data.frame(
            do.call(
                rbind,
                tapply(
                    df_e$value, df_e[1], function(i){quantile(i)})),
            t=1:60)
            qs_e$t = df_t[1:60,1]

            df_scc <- df_scc[df_prtp$value==0.015,]
            df_scc <- df_scc[!is.infinite(df_scc$value),]
            df_scc <- df_scc[!is.nan(df_scc$value),]
            #df_scc <- df_scc[df_damage$value==0.00480515,]
            qs_scc = data.frame(
            do.call(
                rbind,
                tapply(
                    df_scc$value, df_scc[1] , function(i){quantile(i)})),
            t=1:60)
            qs_scc$t = df_t[1:60,1]

            df_cpc <- df_cpc[df_prtp$value==0.015,]
            df_cpc <- df_cpc[!is.infinite(df_cpc$value),]
            df_cpc <- df_cpc[!is.nan(df_cpc$value),]
            #df_scc <- df_scc[df_damage$value==0.00480515,]
            qs_cpc = data.frame(
            do.call(
                rbind,
                tapply(
                    df_cpc$value, df_cpc[1] , function(i){quantile(i)})),
            t=1:60)
            qs_cpc$t = df_t[1:60,1]
        #
    #####
    #####
    ##### GETTING QUANTILES (END)
    #####
    #####


#####
#####
##### Data preparation (end)
#####
#####

#####
#####
##### FIGURE 2: SCC, T, E WITH QUANTILES (START)
#####
#####

  df_r$neworder <- df_r$variable
  df_r$neworder[df_r$variable=='TATM_standard'] <- 1
  df_r$neworder[df_r$variable=='TATM_UVmkt'] <- 2
  df_r$neworder[df_r$variable=='TATM_UV'] <- 3
  df_r$neworder[df_r$variable=='TATM_UVnonUV'] <- 4
  df_r <- df_r[order(df_r$neworder),]
    df_r$value_scc[df_r$years==2020]
    df_r$value_scc[df_r$years==2020]/df_r$value_scc[df_r$years==2020][1]
    df_r$value_t[df_r$years==2100]
    df_r$value_e[df_r$years==2100]
  theme_set(theme_classic())
  p_t <- ggplot(data = df_r, aes(years)) +
        geom_ribbon(data=qs_t, aes(x=t, ymin=X25., ymax=X75.),fill="gray30", alpha=0.2) + 
        geom_line(data = df_r, aes(x=years, y=value_t, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title=" ", y="Temperature (Degrees C)", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(0.8,3.25)) +
        #scale_linetype_manual("", values=c(4,3,2,2,1), labels=c("non-UV", "standard", "UV", "UV-mkt", "UV & non-UV")) +
        #scale_colour_manual("",values=c("violet","indianred","turquoise","blue","seagreen3"),labels=c("non-UV", "standard", "UV","UV-mkt","UV & non-UV")) 
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
  
  p_t

 
  p_e <- ggplot(data = df_r, aes(years)) +
        geom_ribbon(data=qs_e, aes(x=t, ymin=X25., ymax=X75.),fill="gray30", alpha=0.2) +
        geom_line(data = df_r, aes(x=years, y=value_e, group =neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title=" ", y="Emissions (GtCO2)", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(0,40)) +
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
    p_e

  p_scc <- ggplot(data = df_r, aes(years)) +
        geom_ribbon(data=qs_scc, aes(x=t, ymin=X25., ymax=X75.),fill="gray30", alpha=0.2) +
        geom_line(data = df_r, aes(x=years, y=value_scc, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title=" ", y="Opt. Carbon tax (USD/tCO2)", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(10,700)) +
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
    p_scc
  figure <- ggarrange(p_scc, ggarrange(p_e, p_t,
                      labels = c("b", "c"),
                      ncol = 2, legend = FALSE), labels = "a", nrow = 2, common.legend = TRUE, legend = "bottom")
  figure
  ggsave("F2_T_E_SCC.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)

  pdf(file = paste(mypath,"/Figures/F2_T_E_SCC.pdf", sep=""))
  figure
  dev.off()

#####
#####
##### FIGURE 2: SCC, T, E WITH QUANTILES (END)
#####
#####

#####
#####
##### FIGURE 3: Sensitivity (start)
#####
#####
  DF$value_ratio = 1/DF$value_ratio
  DF_2020 <- DF[which(DF[1]==2020),]
  DF_2100 <- DF[which(DF[1]==2100),]
    glimpse(DF_2020)
    glimpse(DF_2100)
  #order sensitivity (start)
    factor(DF_2020$sensitivity)
    DF_2020$value_gama3
    DF_2020$sensitivity <- factor(DF_2020$sensitivity, levels = c("elasmu","prtp","cs","damage","gama3","gama4","ratio","share2","theta2","share1","theta1"))
    type = cbind(c("reference","reference","reference","production","production","production","production","utility","utility","utility","utility"),c("elasmu","prtp","cs","damage","gama3","gama4","ratio","share2","theta2","share1","theta1"))
    type = data.frame(type)
    names(type)[2] = "sensitivity"
    DF_2020 <- merge(DF_2020, type, by="sensitivity")
    names(DF_2020)[dim(DF_2020)[2]] = "type"
    head(DF_2020)
    head(DF_2100)

    DF_2100$sensitivity <- factor(DF_2100$sensitivity, levels = c("elasmu","prtp","cs","damage","gama4","gama3","ratio","share2","theta2","share1","theta1"))
    DF_2100 <- merge(DF_2100, type, by="sensitivity")
    names(DF_2100)[dim(DF_2100)[2]] = "type"
  #order sensitivity (*End)
  

  cbp1 <- c( "#0072B2","#009E73","#E69F00")
  #new names
  DF_2020$names <- DF_2020$sensitivity
  
  DF_2020$names <- revalue(DF_2020$names, c("cs"="Climate sensitivity", "damage"="Damage to \n Natural Capital", "elasmu"="Relative risk \n aversion","gama4" = "Elasticity of non-use values \n to Natural Capital","gama3"="Natural Capital-adjusted \n total factor productivity", 
    "prtp"="Pure rate of \n time preference","ratio"="Natural Capital \n initial stock", "share1"="Ecosystem services \n initial value","share2"="Non-use value initial value", "theta1"="Substitutability between \n market and ES goods","theta2"="Substitutability between \n use and non-use values"))

  DF_2020$intensity <- numeric(dim(DF_2020)[1])  
  glimpse(DF_2020)
  variable = c("cs","damage","elasmu","gama3","gama4","prtp","ratio","share1","share2","theta1","theta2")
  for (i in 1:length(variable)[1]) { #calculating how high is each value
    value_of_interest = variable[i]
    column_of_interest = DF_2020[names(DF_2020)==paste('value_',variable[i],sep="")]
    spread = column_of_interest[DF_2020$sensitivity==value_of_interest,1]
    DF_2020$value_scc[DF_2020$sensitivity==value_of_interest]
    minval = min(spread)
    maxval = max(spread)
    intensity_i = (spread - minval) / (maxval-minval)
    length(intensity_i)
    if (i ==1) {
    intensity = intensity_i
    } else {
        intensity = c(intensity,intensity_i)
    }
  }
  DF_2020$intensity <- intensity
  library(viridis)
  a <- ggplot(DF_2020, aes(x = names, y = value_scc, group = names,fill = intensity)) +
      geom_jitter(width = .005, alpha = 0.7, size = 5, shape = 23) +
      theme_bw() +
      labs(
        y = "Opt. Carbon tax 2020 (USD/tonCO2)",
        x = ""
      ) +
      geom_hline(yintercept=df_r$value_scc[df_r$variable == 'TATM_UVnonUV' & df_r$years==2020] ) + #SCC value of GreenDICE
      geom_hline(aes(yintercept=df_r$value_scc[df_r$variable == 'TATM_standard' & df_r$years==2020]),linetype = "dashed" ) + #SCC value of Standard DICE in 2020
      coord_flip()   +
      #scale_color_viridis(option = "D")
      scale_color_gradientn(colours = cbp1) +
      scale_fill_gradientn(colours = cbp1,labels = c('low','','','','high'))
    #a <- a + scale_fill_manual(values = cbp1, labels=c("production","utilitys","reference"))
    a <- a + labs(fill="Relative value")
    a
    DF_2100$intensity <-  numeric(dim(DF_2100)[1]) 
  variable = c("cs","damage","elasmu","gama4","gama3","prtp","ratio","share1","share2","theta1","theta2")
  for (i in 1:length(variable)[1]) { #calculating how high is each value
    value_of_interest = variable[i]
    print(value_of_interest)
    column_of_interest = DF_2100[names(DF_2100)==paste('value_',variable[i],sep="")]
    spread = column_of_interest[DF_2100$sensitivity==value_of_interest,1]
    minval = min(spread)
    print(minval)
    maxval = max(spread)
    print(maxval)
    intensity_i = (spread - minval) / (maxval-minval)
    length(intensity_i)
    if (i ==1) {
    intensity = intensity_i
    } else {
        intensity = c(intensity,intensity_i)
    }
  }
  DF_2100$intensity <- intensity
    c <- ggplot(DF_2100, aes(x = sensitivity, y = value_t, group = sensitivity, fill = intensity)) +
      geom_jitter(width = .005, alpha = 0.7, size = 5, shape = 23) +
      guides(fill = "none") +
      theme_bw() +
      labs(
        y = "Temperature 2100 (Degrees C)",
        x = ""
      ) +
      geom_hline(yintercept=df_r$value_t[df_r$variable == 'TATM_UVnonUV' & df_r$years==2100] ) + #SCC value of GreenDICE
      geom_hline(aes(yintercept=df_r$value_t[df_r$variable == 'TATM_standard' & df_r$years==2100]),linetype = "dashed" ) + #SCC value of Standard DICE in 2020
        coord_flip()+
      theme(axis.text.y = element_blank()) +
      scale_color_gradientn(colours = cbp1) +
      scale_fill_gradientn(colours = cbp1,labels = c('low','','','','high'))
      
  c <- c + labs(fill="Relative value")
    c
  #c <- c + scale_fill_manual(values = cbp1)



                      #ggsave("sensitivity_UVnonUV_3.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
                      figure <- ggarrange(a,c, labels = c("","",""), ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", widths=c(1.5,1))
                      figure
                      #ggsave("F3_sensitivity_GreenDICE.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
                       pdf(file = paste(mypath,"/Figures/F3_sensitivity_GreenDICE.pdf", sep=""))
                        figure
                      dev.off()

  #PArametric sensitivity UVnonUV (end)
#####
#####
##### FIGURE 3: Sensitivity (end)
#####
#####

#####
#####
##### FIGURE 4: Minimal depth of random forests (start)
#####
#####

    #function
        grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {

        plots <- list(...)
        position <- match.arg(position)
        g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
        legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
        lheight <- sum(legend$height)
        lwidth <- sum(legend$width)
        gl <- lapply(plots, function(x) x + theme(legend.position="none"))
        gl <- c(gl, ncol = ncol, nrow = nrow)

        combined <- switch(position,
                        "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                                legend,
                                                ncol = 1,
                                                heights = unit.c(unit(1, "npc") - lheight, lheight)),
                        "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                                legend,
                                                ncol = 2,
                                                widths = unit.c(unit(1, "npc") - lwidth, lwidth)))

        grid.newpage()
        grid.draw(combined)

        # return gtable invisibly
        invisible(combined)

    }


    #RF for SCC in 2020 (start)
        df_mcs_2020 <- df_mcs[which(df_mcs[1]==2020),]
        glimpse(df_mcs_2020)
        names(df_mcs_2020) =c("years","model","value_t","value_e","value_scc","Natural Capital-adjusted \n total factor productivity","Natural Capital \n initial stock","Damage to \n Natural Capital","Ecosystem services \n initial value","Non-use values initial value","Substitutability between \n use and non-use values","Substitutability between \n market and ES goods","Elasticity of non-use \n values to Natural Capital","Pure rate of \n time preference","Climate sensitivity")
        df_mcs_2020 <- df_mcs_2020[,5:13]
        trf <- tuneRF(df_mcs_2020[,2:9], df_mcs_2020[,1])
        mt <- trf[which.min(trf[,2]), 1]
        Results_rf <- randomForest(df_mcs_2020[,2:9], df_mcs_2020[,1], importance = TRUE,tree = TRUE, mtry =mt)
        Results_rf1 <- Results_rf 
        min_depth_frame <- min_depth_distribution(Results_rf)
        md1 <- plot_min_depth_distribution(min_depth_frame, mean_sample = "relevant_trees", k = 15)
        #ggsave("RF_SCC2020_minDepth.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
        md1 <- md1 + labs(title="Opt. Carbon tax 2020", y="Number of trees", x = "", color="")
        md1 <- md1 + theme_minimal()
        
    #RF for SCC in 2020 (end)
    #RF for T in 2100 (start)
        df_mcs_2100 <- df_mcs[which(df_mcs[1]==2100),]
        head(df_mcs_2100)
        names(df_mcs_2100) = c("years","model","value_t","value_e","value_scc","Natural Capital-adjusted \n total factor productivity","Natural Capital \n initial stock","Damage to \n Natural Capital","Initial ecosystem \n services value","Initial non-use value","Substitutability between \n use and non-use values","Substitutability between \n market and ES goods","Non-use values elasicity \n to Natural Capital","Pure rate of \n time preference","Climate sensitivity")
        df_mcs_t2100 <- cbind(df_mcs_2100$value_t,df_mcs_2020[2:9])
        head(df_mcs_t2100)
        names(df_mcs_t2100)[1] = "value_t2100"
        trf <- tuneRF(df_mcs_t2100[,2:9], df_mcs_t2100[,1])
        mt <- trf[which.min(trf[,2]), 1]
        mt
        Results_rf <- randomForest(df_mcs_t2100[,2:9], df_mcs_t2100[,1], importance = TRUE,tree = TRUE, mtry =mt)
        Results_rf$importance
        plot(Results_rf)
        varImpPlot(Results_rf) #%IncMSE is the most robust and informative measure. It is the increase in mse of predictions(estimated with out-of-bag-CV) as a result of variable j being permuted(values randomly shuffled).
        varImpPlot(Results_rf, type = 1) 
        min_depth_frame <- min_depth_distribution(Results_rf)
        md2 <- plot_min_depth_distribution(min_depth_frame, mean_sample = "relevant_trees", k = 15)
        md2 <- md2 + labs(title="Temperature 2100", y="Number of trees", x = "", color="") + theme_minimal()
    #RF for T in 2100 (end)
        
        figure <- ggarrange(md1,md2, labels = c("","",""), ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
        annotate_figure(figure, top = NULL, bottom = NULL, left = NULL,
    right = NULL, fig.lab = "", fig.lab.pos = c("top.left"))
                        figure
                        
                        ggsave("F4_minimal_depth_panel.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
            pdf(file = paste(mypath,"/Figures/F4_minimal_depth_panel.pdf", sep=""))
                        figure
                      dev.off()
    #RF for SCC in 2100 (end)

#####
#####
##### FIGURE 4: Minimal depth of random forests (end)
#####
#####


#####
#####
##### FIGURE 5: Investments (start)
#####
#####


  ####################################
  ######################################

  ## Comparison between investments (start)


    df_r_inv_standard = df_r[which(df_r$variable=="TATM_standard" ),]
    df_r_inv_green = df_r[which(df_r$variable=="TATM_UVnonUV" ),]
    
    #comp_inv <- rbind(df_r_inv_standard,df_r_inv_green,df_inv_last,df_inv_reduceddamages)

    dif <- df_inv_last$value_YGross - df_r_inv_green$value_YGross
    df_inv_last <- cbind(df_inv_last,dif)

    dif <- df_inv_reduceddamages$value_YGross - df_r_inv_green$value_YGross
    df_inv_reduceddamage <- cbind(df_inv_reduceddamages,dif)


    
    plot_t2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          #geom_line(data = df_inv_last, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_reduceddamage, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          labs(title="Temperature", y="Degrees C", x = "years") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0.8,3.2)) +
          #scale_linetype_manual("", values=c(4,2,3,1), labels=c( "Asset Investment", "Damage Reduction", "Standard DICE", "GreenDICE")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan","indianred","seagreen3"),labels=c( "Asset Investment","Damage Reduction", "Standard DICE", "GreenDICE")) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Adaptive investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Adaptive investments", "Standard DICE", "GreenDICE")) 
         plot_t2

    plot_scc2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          #geom_line(data = df_inv_last, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_reduceddamage, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          labs(title="Opt. Carbon tax", y="USD/tonCO2", x = "years") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0,700)) +
          #scale_linetype_manual("", values=c(4,2,3,1), labels=c( "Asset Investment", "Damage Reduction", "Standard DICE", "GreenDICE")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan","indianred","seagreen3"),labels=c( "Asset Investment","Damage Reduction", "Standard DICE", "GreenDICE")) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Adaptive investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Adaptive investments", "Standard DICE", "GreenDICE")) 
         plot_scc2

     plot_inv2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_inv*100, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_inv*100, colour = variable, linetype=variable),size=1)  + 
          #geom_line(data = df_inv_last, aes(x=years, y=value_inv*value_YGross*1000, colour = variable, linetype=variable),size=1)  + 
          #geom_line(data = df_inv_reduceddamage, aes(x=years, y=value_inv*value_YGross*1000, colour = variable, linetype=variable),size=1)  + 
          #geom_line(data = df_inv_last, aes(x=years, y=value_inv, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_reduceddamage, aes(x=years, y=value_inv*100, colour = variable, linetype=variable),size=1)  + 
          labs(title="Adaptive investments", y="Percent of GWP", x = "years") +
          #coord_cartesian(xlim = c(2010, 2100),ylim=c(0,600)) +
          #scale_linetype_manual("", values=c(4,2), labels=c( "Asset Investment", "Damage Reduction")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan"),labels=c( "Asset Investment","Damage Reduction")) +
          #scale_linetype_manual("", values=c(4,2,3,1), labels=c( "Asset Investment", "Damage Reduction", "Standard DICE", "GreenDICE")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan","indianred","seagreen3"),labels=c( "Asset Investment","Damage Reduction", "Standard DICE", "GreenDICE")) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Adaptive investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Adaptive investments", "Standard DICE", "GreenDICE")) +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0.01,1))
    plot_inv2

    
    plot_nc2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          #geom_line(data = df_inv_last, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_reduceddamage, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          labs(title="Natural capital stock", y="Value (USD)", x = "years") +
          #scale_linetype_manual("", values=c(4,2,3,1), labels=c( "Asset Investment", "Damage Reduction", "Standard DICE", "GreenDICE")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan","indianred","seagreen3"),labels=c( "Asset Investment","Damage Reduction", "Standard DICE", "GreenDICE")) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Adaptive investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Adaptive investments", "Standard DICE", "GreenDICE")) +
         coord_cartesian(xlim = c(2010, 2100),ylim=c(0,40)) 
    plot_nc2 
    
    comparison_investments <- ggarrange(plot_inv2,plot_nc2,plot_scc2,plot_t2,ncol=2,nrow=2, 
      labels = c("A", "B", "C","D"),legend="bottom", common.legend=TRUE)
    comparison_investments

    pdf(file = paste(mypath,"/Figures/F5_adaptiveinvestments.pdf", sep=""))
                       figure
                   dev.off()

    ggsave("F5_investments.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
  ## Comparison between investments (end)

#####
#####
##### FIGURE 5: Investments (END)
#####
#####



#####
#####
##### FIGURE S.1.: SCC in 2020 by NC ratio (START)
#####
#####
  #Generate new dataframes (start)
    newdf <- df_nc[which(df_nc$time_UVnonUV_combination_NC_tfp==2020),]
    newdf2100 <- df_nc[which(df_nc$time_UVnonUV_combination_NC_tfp==2100),]
    head(newdf2100)
    names(newdf2100)[3] = "value_t2100"
    newdf <- cbind(newdf,newdf2100$value_t2100)
    names(newdf)[13] = "value_t2100"

    newdf[14] <- 1 / newdf$value_ratio
    names(newdf)[14] = "NC_K"
    newdf_default = newdf[which(newdf$value_gama3==newdf$value_gama3[1]),]
    newdf_default = newdf_default[which(newdf_default$value_ratio==newdf_default$value_ratio[1]),]
  #Generate new dataframes (end)

    gscc <- ggplot(newdf, aes(x = value_atfp,y = value_scc))
    gscc <- gscc + geom_jitter(aes( size = NC_K, fill = NC_K), shape = 21, alpha = 0.7,colour = "transparent")  + 
    scale_fill_gradientn(colours = cbp1, guide = "legend", name = expression(paste( over(NC[0],K[0])))) +
    scale_size_continuous(range = c(1, 8), name = expression(paste(over(NC[0],K[0])))) +
    geom_smooth(method="lm", se=T) +
    labs(title=" ", y="Temperature in 2100", x = "Natural capital-adjusted total factor productivity", color="") +
    guides(color = FALSE) + 
    geom_point(aes(x =  1.014353465 , y = df_r$value_scc[which(df_r$variable=="TATM_UVnonUV")[3]]), shape = "*", color = "red", size = 15)
    gscc <- gscc + theme_bw() 
    gscc

    gT <- ggplot(newdf, aes(x = value_atfp,y = value_t2100))
        gT <- gT + geom_jitter(aes( size = NC_K, fill = NC_K), shape = 21, alpha = 0.7,colour = "transparent")  + 
        scale_fill_gradientn(colours = cbp1, guide = "legend", name = expression(paste( over(NC[0],K[0])))) +
        scale_size_continuous(range = c(1, 8), name = expression(paste(over(NC[0],K[0])))) +
        geom_smooth(method="lm", se=T) +
        labs(title=" ", y="Temperature in 2100", x = "Natural capital-adjusted total factor productivity", color="") +
        guides(color = FALSE) + 
        geom_point(aes(x =  1.014353465 , y = df_r$value_t[which(df_r$variable=="TATM_UVnonUV")[19]]), shape = "*", color = "red", size = 15)
    gT <- gT + theme_bw() 
        gT

    gT <- gT +labs(title=" ", y="Temperature in 2100 (Degrees C)", x = "Natural capital-adjusted total factor productivity", color="Natural Capital \n initial stock") +
        theme(legend.position="top") 
    gT
    legend <- get_legend(gT)
    
    gT <- gT + theme(legend.position="none",plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"))
    gT
    gT <- ggMarginal(gT, type = "histogram", margins = "y", fill = "gray", alpha = 0.7)
    gT

    gS <- gscc + labs(title=" ", y="Opt C tax in 2020 (USD/tonCO2)", x = "", color="") +
        theme(legend.position="none", plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines")) + scale_color_gradientn(colours = cbp1)

    gS
    gS <- ggMarginal(gS, type = "histogram", margins = "y", fill = "gray", alpha = 0.7)
    gS

    myplotSCC <- arrangeGrob(gS, top = textGrob("A", x = unit(0, "npc")
        , y   = unit(1, "npc"), just=c("left","top"),
        gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))

    myplotT <- arrangeGrob(gT, top = textGrob("B", x = unit(0, "npc")
        , y = unit(1, "npc"), just=c("left","top"),
        gp=gpar(col="black", fontsize=18, fontfamily="Times Roman")))

    png("C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures/FS1_Temp2100_SCC2020_NC_gama_2.png", width = 6, height = 6, units = 'in', res = 600)
    
    fig <- grid.arrange(myplotSCC, myplotT, legend, nrow = 3,heights = c(2.5, 2.5,0.5))
    dev.off()

    pdf(file = paste(mypath,"/Figures/FS1_Temp2100_SCC2020_NC_gama_2.pdf", sep=""))
                       fig
                   dev.off()
    ##Revision (Figure NCTFP)
    ##### TEMP2100 sensitivity to NC and Gama (end)
#####
#####
##### FIGURE S.1.: SCC in 2020 by NC ratio (END)
#####
#####

#####
#####
##### FIGURE S.2.: Cost of damages reduction
#####
#####
    theme_set(theme_classic())


    plot_inv_nc <- ggplot(data = df_inv_reduceddamages_sens, aes(years)) +
          geom_line(data = df_inv_reduceddamages_sens, aes(x=years, y=value_nc, group = variable, color = id_var),size=1)  + 
          #geom_line(data = df_inv, aes(x=years, y=value_k, group = variable, colour = 'red'),size=1)  + 
          labs(title="Natural Capital stocks", y="USD", x = "years", color="") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0,40)) +
          scale_color_gradientn(colours = cbp1, labels = c('low','','','','high')) + labs(color="cost of damage reduction") 
          #scale_linetype_manual("", values=c(3,1,2), labels=c( "standard DICE", "GreenDICE","GreenDICE + investment")) +
          #scale_colour_manual("",values=c("indianred","seagreen3","blue"),labels=c("standard DICE", "GreenDICE", "GreenDICE + investment")) + 
          #theme(legend.position="")
    plot_inv_nc
    df_r_inv_standard = df_r[which(df_r$variable=="TATM_standard" ),]
    df_r_inv_green = df_r[which(df_r$variable=="TATM_UVnonUV" ),]
    plot_inv_nc2 <- plot_inv_nc +
    #geom_line(data = df_inv_reduceddamages, aes(x=years, y=value_nc, group = variable),size=1,color="red",linetype=2) +
    geom_line(data = df_r_inv_standard, aes(x=years, y=value_nc, group = variable),size=1,color="darkgray",linetype=3) +
    geom_line(data = df_r_inv_green, aes(x=years, y=value_nc, group = variable),size=1,color="darkgray",linetype=4) 
    #scale_colour_manual("",values=c("indianred","seagreen3","blue"),labels=c("standard DICE", "standard GreenDICE")) 
    plot_inv_nc2


    plot_inv <- ggplot(data = df_inv_reduceddamages_sens, aes(years)) +
          geom_line(data = df_inv_reduceddamages_sens, aes(x=years, y=value_s*100, group = variable, colour = id_var),size=1, linetype=1)  + 
          geom_line(data = df_inv_reduceddamages_sens, aes(x=years, y=value_inv*100, group = variable, colour = id_var),size=1, linetype=2)  + 
          #geom_line(data = r_t_inv, aes(x=names(r_t_inv)[1], y=value_t, group = variable, colour = variable, linetype=variable),size=1)  + 
          labs(title="Investments in damage reduction", y="Percentage of GWP", x = "years", color="") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0.01,0.5)) +
          scale_color_gradientn(colours = cbp1, labels = c('low','','','','high')) + labs(color="cost of damage reduction") 
          #scale_colour_manual("",values=c("indianred","seagreen3","blue"),labels=c("standard DICE", "GreenDICE", "GreenDICE + investment")) + 
          theme(legend.position="")
    plot_inv + scale_y_log10()
    plot_inv2 <- plot_inv +
    #geom_line(data = df_inv_reduceddamages, aes(x=years, y=value_inv*100, group = variable),size=1,color="red",linetype=2) 
    geom_line(data = df_r_inv_standard, aes(x=years, y=value_s*100, group = variable),size=1,color="darkgray",linetype=3) +
    geom_line(data = df_r_inv_green, aes(x=years, y=value_s*100, group = variable),size=2,color="darkgray",linetype=4) 
    plot_inv2 <- plot_inv2
          plot_inv2

    plot_scc <- ggplot(data = df_inv_reduceddamages_sens, aes(years)) +
          geom_line(data = df_inv_reduceddamages_sens, aes(x=years, y=value_scc, group = variable, colour = id_var),size=1)  + 
          #geom_line(data = df_inv, aes(x=years, y=value_inv, group = variable, colour = id_var),size=1)  + 
          #geom_line(data = r_t_inv, aes(x=names(r_t_inv)[1], y=value_t, group = variable, colour = variable, linetype=variable),size=1)  + 
          labs(title="Optimal Carbon tax", y="USD/tonCO2", x = "years", color="") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(50,750)) +
          scale_color_gradientn(colours = cbp1, labels = c('low','','','','high')) + labs(color="cost of damage reduction") 
           #scale_linetype_manual("", values=c(3,1,2), labels=c( "standard DICE", "GreenDICE","GreenDICE + investment")) +
          #scale_colour_manual("",values=c("indianred","seagreen3","blue"),labels=c("standard DICE", "GreenDICE", "GreenDICE + investment")) + 
          theme(legend.position="")
    plot_scc
    plot_scc2 <- plot_scc +
    geom_line(data = df_r_inv_standard, aes(x=years, y=value_scc, group = variable),size=1,color="darkgray",linetype=3) +
    geom_line(data = df_r_inv_green, aes(x=years, y=value_scc, group = variable),size=1,color="darkgray",linetype=4) 
    plot_scc2


    plot_t <- ggplot(data = df_inv_reduceddamages_sens, aes(years)) +
          geom_line(data = df_inv_reduceddamages_sens, aes(x=years, y=value_t, group = variable, colour = id_var),size=1)  + 
          #geom_line(data = df_r_inv_standard, aes(x=years, y=value_t, group = variable),size=1,color="indianred",linetype=3) +
            #geom_line(data = df_r_inv_green, aes(x=years, y=value_t, group = variable),size=1,color="seagreen3",linetype=4) +
            labs(title="Temperature", y="degrees C", x = "years", color="",linetype="") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0.8,3)) +
          scale_color_gradientn(colours = cbp1, labels = c('low','','','','high')) + labs(color="cost of damage reduction") 
          #scale_linetype_manual(values=c(3,1,2,1,1), labels=c( "standard DICE", "GreenDICE","GreenDICE + investment","a","b")) 
          #scale_colour_manual("",values=c("indianred","seagreen3","blue"),labels=c("standard DICE", "GreenDICE", "GreenDICE + investment")) + 
          #labels = c('low','','','','high')
    plot_t
    plot_t2 <- plot_t +
    geom_line(data = df_r_inv_standard, aes(x=years, y=value_t, group = variable),size=1,color="darkgray",linetype=3) +
    geom_line(data = df_r_inv_green, aes(x=years, y=value_t, group = variable),size=1,color="darkgray",linetype=4) 
          plot_t2

    fig <- ggarrange(plot_t2,plot_scc2,plot_inv,plot_inv_nc2,ncol=2,nrow=2, common.legend=TRUE)
    fig
    ggsave("FS2_ReducedDam_NC_sens.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
    pdf(file = paste(mypath,"/Figures/FS2_ReducedDam_NC_sens.pdf", sep=""))
                       fig
                   dev.off()
  #figures of iterations (end)
#####
#####
##### FIGURE S.2.: Cost of damages reduction
#####
#####


#####
#####
##### FIGURE S.3.: Iterations in investments (START)
#####
#####



     plot_t2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_last, aes(x=years, y=value_t, colour = variable, linetype=variable),size=1)  + 
          labs(title="Temperature", y="Degrees C", x = "years") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0.8,3.2)) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Asset investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Asset investments", "Standard DICE", "GreenDICE")) 
         plot_t2

    plot_scc2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_last, aes(x=years, y=value_scc, colour = variable, linetype=variable),size=1)  + 
          labs(title="Opt. Carbon tax", y="USD/tonCO2", x = "years") +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0,700)) +
          #scale_linetype_manual("", values=c(4,2,3,1), labels=c( "Asset Investment", "Damage Reduction", "Standard DICE", "GreenDICE")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan","indianred","seagreen3"),labels=c( "Asset Investment","Damage Reduction", "Standard DICE", "GreenDICE")) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Asset investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Asset investments", "Standard DICE", "GreenDICE")) 
         plot_scc2

     plot_inv2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_inv*100, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_inv*100, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_inv_last, aes(x=years, y=value_inv*100, colour = variable, linetype=variable),size=1)  + 
          labs(title="Asset investments", y="Percent of GWP", x = "years") +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Asset investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Asset investments", "Standard DICE", "GreenDICE")) +
          coord_cartesian(xlim = c(2010, 2100),ylim=c(0.01,50))
    plot_inv2

    
    plot_nc2 <- ggplot(data =  df_r_inv_standard, aes(years)) +
          geom_line(data = df_inv_last, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_standard, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          geom_line(data = df_r_inv_green, aes(x=years, y=value_nc, colour = variable, linetype=variable),size=1)  + 
          labs(title="Natural capital stock", y="Value (USD)", x = "years") +
          #scale_linetype_manual("", values=c(4,2,3,1), labels=c( "Asset Investment", "Damage Reduction", "Standard DICE", "GreenDICE")) +
          #scale_colour_manual("",values=c("firebrick2","darkcyan","indianred","seagreen3"),labels=c( "Asset Investment","Damage Reduction", "Standard DICE", "GreenDICE")) +
          scale_linetype_manual("", values=c(2,3,1), labels=c("Asset investments", "Standard DICE", "GreenDICE")) +
          scale_colour_manual("",values=c("darkcyan","indianred","seagreen3"),labels=c( "Asset investments", "Standard DICE", "GreenDICE")) +
         coord_cartesian(xlim = c(2010, 2100),ylim=c(0,250)) 
    plot_nc2 
    fig <- ggarrange(plot_t2,plot_scc2,plot_inv2,plot_nc2,ncol=2,nrow=2,legend="bottom", common.legend=TRUE)
    fig    
    ggsave("FS3_investment_NC_iterations.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
    pdf(file = paste(mypath,"/Figures/FS3_investment_NC_iterations.pdf", sep=""))
                       fig
                   dev.off()
  
  #figures of iterations (end)
#####
#####
##### FIGURE S.3.: Iterations in investments (END)
#####
#####

#####
#####
##### FIGURE S.4.: Other variables(START)
#####
#####
 p_gwp <- ggplot(data = df_r, aes(years)) +
        geom_line(data = df_r, aes(x=years, y=value_YGross, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title="Economic Output", y="trill USD", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(80,500)) +
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
  
  p_s <- ggplot(data = df_r, aes(years)) +
        geom_line(data = df_r, aes(x=years, y=value_s, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title="Capital savings", y="Fraction of Gross Output", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2200),ylim=c(0.225,0.275)) +
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
  
  
  p_nc <- ggplot(data = df_r, aes(years)) +
        geom_line(data = df_r, aes(x=years, y=value_nc, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title="Natural Capital Stock", y="trill USD", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(0,40)) +
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
    p_nc
    p_k <- ggplot(data = df_r, aes(years)) +
        geom_line(data = df_r, aes(x=years, y=value_k, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title="Manufactured Capital Stock", y="trill USD", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(200,1000)) +
        scale_linetype_manual("", values=c(3,4,2,1), labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) +
        scale_colour_manual("",values=c("indianred","violet","blue","seagreen3"),labels=c("DICE", "Market Only", "All Use Values", "GreenDICE")) 
    p_k

     figure <- ggarrange(p_s,p_gwp,p_nc,p_k, labels = c("a","b","c","d"), common.legend = TRUE, legend = "bottom")
  figure
  #ggsave("FS4_Savings_GWP.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
  pdf(file = paste(mypath,"/Figures/FS4_Savings_GWP.pdf", sep=""))
                       fig
                   dev.off()
  
  p_cpc <- ggplot(data = df_cpc, aes(years)) +
        geom_ribbon(data=qs_cpc, aes(x=t, ymin=X25., ymax=X75.),fill="gray30", alpha=0.2)  +
        geom_ribbon(data=qs_cpc, aes(x=t, ymin=X0., ymax=X100.),fill="gray30", alpha=0.2) +
        #geom_line(data = df_r, aes(x=years, y=value_cpc, group = neworder, colour = neworder, linetype=neworder),size=1)  + 
        labs(title="Consumption per Capita", y="thousand USD/ person", x = "years", color="") +
        coord_cartesian(xlim = c(2020, 2100),ylim=c(5,40)) 
  p_cpc        
  
  min_cpc <- ddply(df_cpc, "time_nonUV_spread_mcs_opt", summarise, minimum_cpc = min(value_cpc))
  names(min_cpc)[1] = "year"
  min_cpc <- min_cpc[1:19,]
  
  #fi
#####
#####
##### FIGURE S.4.: Other variables(END)
#####
#####


#####
#####
##### FIGURE S.5.: Minimal Depth interactions
#####
#####
  df_mcs_2020 <- df_mcs[which(df_mcs[1]==2020),]
  glimpse(df_mcs_2020)
  names(df_mcs_2020) = c("years","model","value_t","value_e","value_scc","atfp","NC0/K0","damage","share1","share2","theta2","theta1","gama4","prtp","cs")
  df_mcs_2020 <- df_mcs_2020[,5:13]
  trf <- tuneRF(df_mcs_2020[,2:9], df_mcs_2020[,1])
  mt <- trf[which.min(trf[,2]), 1]
  Results_rf <- randomForest(df_mcs_2020[,2:9], df_mcs_2020[,1], importance = TRUE,tree = TRUE, mtry =mt)
  importance_frame <- measure_importance(Results_rf)
  (vars <- important_variables(importance_frame, k = 5, measures = c("mean_min_depth", "no_of_trees")))
  interactions_frame <- min_depth_interactions(Results_rf, vars, mean_sample = "relevant_trees", uncond_mean_sample = "relevant_trees")
  interactions_frame <- interactions_frame[-which(interactions_frame$interaction %in% c("damage:damage")),]
  plot_min_depth_interactions(interactions_frame)
  #ggsave("FS5_interaction.png", path="C:/Users/bastien/Documents/GitHub/GreenDICE/Results/Figures", dpi=600)
    pdf(file = paste(mypath,"/Figures/FS5_interaction.pdf", sep=""))
                       fig
                   dev.off()
  
#####
#####
##### FIGURE S.5.: Minimal Depth interactions
#####
#####


