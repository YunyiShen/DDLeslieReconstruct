## Model predict percision
nage = 8
res_link = "E:/UW Lab jobs/3. Lislie Reconstruction/Lislie_Reconstruct_ChicagoDeer/figs/temp"
res_link = "C:/Users/yshen99/Documents/GitHub/Lislie_Reconstruct_ChicagoDeer/figs/temp" 
### MSE & MPA
n_data = 0
MSE = 0
MPA = 0
MAD = 0
MAD_vec = numeric()

for(i in 1:nage){
  file_temp = paste0(res_link,"/age",i,".csv")
  temp = read.csv(file_temp)
  realdata = temp$point=="data"
  realvalue = temp[realdata,]
  model_predict = temp[realdata==0,]
  
  n_sub_data = sum(realdata)
  
  MSE = MSE + sum((model_predict$mean-realvalue$mean)^2)
  MAD = MAD + sum(abs(model_predict$mean-realvalue$mean))
  MAD_vec = c(MAD_vec,abs(model_predict$mean-realvalue$mean))
  
  MPA1 = apply(matrix(1:n_sub_data),1,function(k,realvalue,prediction){
    (realvalue$mean[k]>=prediction$low[k])&(realvalue$mean[k]<=prediction$high[k])
  }
              ,realvalue,model_predict)
  
  MPA = MPA + sum(MPA1)
  
  n_data = n_data + n_sub_data
}
MSE = MSE/n_data
MPA = MPA/n_data
MAD = MAD/n_data # mean absolute difference
SE_MAD = sd(MAD_vec)/sqrt(n_data)

### MSD # mean standard diviation
MSD = mean(apply(Chicago_RES$lx.mcmc,2,sd))
SE_MSD = sd(apply(Chicago_RES$lx.mcmc,2,sd))/sqrt(nrow(Chicago_RES$lx.mcmc))
