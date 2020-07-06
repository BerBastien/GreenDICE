#Prefered parameters
share1_param = 0.1
share2_param = 0.1
k_nc_param = 3.87
atfp_param = 1.014353465
gamma4_param = 0.5
elasmu_param = 1.45


#Vectors of estimations (start)
elasmus = [1.08, elasmu_param, 1.82] 
share1_i = [0.05, share1_param, 0.15] 
share2_i = [0.05, share2_param, 0.15] 
theta1_i = [0.74, 0.86, 0.68, 0.69, 0.32, 0.58, -0.16, 0.41, -0.1, 0.73, 0.79, 0.76, 0.80]
theta2_i = [0.63,0.78,0.62,0.27]
theta1_param = mean(theta1_i)
theta2_param = mean(theta2_i)
damagek_i = [1,2,3]
damaged_i = [0.0026686075, (0.0026686075/1.25)*2]
prtp_i = [0.001,0.015,0.03]
damagek3 = 0.624 #From Nordhaus 1994 Survey
#Vectors of estimations (end)

#Distributions (start)
cs_lnd = LogNormal(log(3.2), 0.12) #based on Roe and Baker 2007, using the mean of Nordhaus
prtp_ud = Uniform(0.001, 0.03) 
tfp_param = Normal(atfp_param,0.028783366) #based on wighted adjusted TFP by country GDP (see gamma3_computation.R)
k_nc_nd = Normal(3.87,2.11)
damagek1_nd = Normal(0.151, 0.147) #From regression using Howard and Sterner (2017) data
atfp_nd =  Normal(atfp_param,0.028783366) #From Brandt et al (2017) weighted mean and standard deviation
gamma4_ud = Uniform(0.2,0.8)
#Distributions (end) 
