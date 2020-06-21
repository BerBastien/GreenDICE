#Vectors of estimations (start)
Adjusted_tfp_i = [1.01144,0.9565217391,1.006993007,1,1.1052631579,1.0090909091,0.96,1.0555555556,1.0070422535,1.0526315789,0.8644067797,0.9705882353,1.0194174757,1.2222222222,1.1237113402,1.1312217195,1.049382716]; #removing repeated 1's, 17 elasticities
k_nc_i = [44760/15841, 1967/6421, 6531/6949, 28527/18960, 59096/80104, 195929/19525] #Ref: Changing Wealth of Naitons, World Bank 2018, Appendix B
latfp = log.(Adjusted_tfp_i)
lknc = log.(k_nc_i)
lknc_inv = 1 ./ lknc
gama3_iM = latfp * lknc_inv'
gama3_i = [gama3_iM[1:end,1];gama3_iM[1:end,2];gama3_iM[1:end,3];gama3_iM[1:end,4];gama3_iM[1:end,5];gama3_iM[1:end,6]]
gama3_i = gama3_i[gama3_i.>0]
gama3_i = gama3_i[gama3_i.<0.2999999]
share1_i = [0.05, 0.1, 0.15] #invented 
share2_i = [0.05, 0.1, 0.15] #invented 
theta1_i = [0.74,0.86,0.68,0.69,0.78,0.32,0.62,0.27,0.58,-0.16,0.41,-0.1,0.73,0.79,0.76,0.80]
theta2_i = [0.63,0.78,0.62,0.27]
damagek_i = [1,2,3]
damaged_i = [0.0026686075,0.0026686075, (0.0026686075/1.25)*2]
prtp_i = [0.001,0.015,0.03]
damagek3 = 0.624 #From Nordhaus 1994 Survey
#Vectors of estimations (end)

#Distributions (start)
cs_lnd = LogNormal(log(3.2), 0.12) #based on Roe and Baker 2007, using the mean of Nordhaus
prtp_ud = Uniform(0.001, 0.03) 
tfp_param = Normal(1.01144,0.000192385) #based on wighted adjusted TFP by country GDP (see gamma3_computation.R)
gama3_nd = Normal(log(1.01144)/log(44760/15841), std(gama3_i)) #normal distributions
k_nc_nd = Normal(3.87,2.11)
damagek1_nd = Normal(0.151, 0.147) #From regression using Howard and Sterner (2017) data
damagek2_ud = Uniform(0.3,0.5) #From Howard and Sylvan (2015)
atfp_nd = Normal(1.0144,0.0324) #From Brandt et al (2017) weighted mean and standard deviation
gamma4_ud = Uniform(0.6,1)

#Distributions (end) 
