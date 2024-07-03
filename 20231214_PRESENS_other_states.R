Data_MA <- read.csv("20231214_PRESENS_MA.csv")
Data_NJ <- read.csv("20231214_PRESENS_NJ.csv")
Data_UT <- read.csv("20231214_PRESENS_UT.csv")

library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(rstan)
library(brms)
library(loo)

data_2 <- Data_UT ##Data名を変更
data_2 <- select(data_2, Date, IAV_0, IAV_nov_0, IAV_case)
colnames(data_2) <- c("collection_date", "IAV", "IAV_normalized", "IAV_day")
data_2 <- mutate(data_2, IAV_normalized = IAV_normalized + 9.05) #PMMoV; 9.05
data_2 <- mutate(data_2, IAV_day = IAV_day * 38.9 / 3.41) #NJ; 9.01, Ut; 3.41, MA: 6.97

#data名変える(移動平均)
data_lead <- mutate(data_2, IAV_1 = lag(IAV, 1), IAV_2 = lag(IAV, 2), IAV_3 = lag(IAV, 3), IAV_4 = lag(IAV, 4), IAV_5 = lag(IAV, 5), IAV_6 = lag(IAV, 6), IAV_7 = lag(IAV, 7), IAV_8 = lag(IAV, 8))
data_behind <- mutate(data_2, IAV_1 = lead(IAV, 1), IAV_2 = lead(IAV, 2),IAV_3 = lead(IAV, 3), IAV_4 = lead(IAV, 4), IAV_5 = lead(IAV, 5))

#Normalized version
data_lead <- mutate(data_2, IAV = IAV_normalized, IAV_1 = lag(IAV_normalized, 1), IAV_2 = lag(IAV_normalized, 2), IAV_3 = lag(IAV_normalized, 3), IAV_4 = lag(IAV_normalized, 4), IAV_5 = lag(IAV_normalized, 5), IAV_6 = lag(IAV_normalized, 6), IAV_7 = lag(IAV_normalized, 7), IAV_8 = lag(IAV_normalized, 8))
data_behind <- mutate(data_2, IAV_1 = lead(IAV_normalized, 1), IAV_2 = lead(IAV_normalized, 2),IAV_3 = lead(IAV_normalized, 3), IAV_4 = lead(IAV_normalized, 4), IAV_5 = lead(IAV_normalized, 5))

sample_size <- nrow(data_lead) #ここのdataファイル

data_WBE <- select(data_lead, collection_date, IAV_2) #ここのdataファイル名と数字を変える
colnames(data_WBE) <- c("collection_date", "IAV")
data_WBE <- data_WBE %>% drop_na()


#Pick up row that contains IAV concentration data
IAV_data_row <- data.frame(true = which(!is.na(data_lead$IAV_2)))#ここのdataファイル名と数字を変える
n1 <- nrow(IAV_data_row)


#artificially dataframe for shedding dynamics
Data_shedding <- read.csv("20231214_PRESENS_shedding.csv")
n2 <- nrow(Data_shedding)

data_list <- list(n = sample_size, IAV = data_2$IAV_day, na = n1, WBE = data_WBE$IAV,  pick_colum = IAV_data_row$true, nb = n2, feces = Data_shedding$Shedding)

mcmc_CA <- stan(
  file = "20231214_PRESENS_other_states.stan",
  data = data_list,
  seed = 1,
  chain = 4,
  iter = 15000,
  warmup = 4000,
  thin = 1
)

print(mcmc_CA, pars = c("k", "sigma", "lp__"), probe = c(0.025, 0.50, 0.975))


log_lik <- extract_log_lik(mcmc_CA)
waic(log_lik)
loo(log_lik)




#prediction

mcmc_sample <- rstan::extract(mcmc_CA)
k <- mcmc_sample[["k"]]
sigma <- mcmc_sample[["sigma"]]

IAV_data_row <- data.frame(true = which(!is.na(data_lead$IAV_2)))
nrow <- nrow(Data_UT) #Data名変更

summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

#31 10/24 3.631 Log copies/L
for(i in 26:nrow){
  WBE = data_2$IAV_day[i]*Data_shedding$Shedding[1] + data_2$IAV_day[i - 1]*Data_shedding$Shedding[2] + data_2$IAV_day[i - 2]*Data_shedding$Shedding[3] + data_2$IAV_day[i - 3]*Data_shedding$Shedding[4] + data_2$IAV_day[i - 4]*Data_shedding$Shedding[5] +
    data_2$IAV_day[i - 5]*Data_shedding$Shedding[6] + data_2$IAV_day[i - 6]*Data_shedding$Shedding[7] + data_2$IAV_day[i - 7]*Data_shedding$Shedding[8] + data_2$IAV_day[i - 8]*Data_shedding$Shedding[9] + data_2$IAV_day[i - 9]*Data_shedding$Shedding[10] +
    data_2$IAV_day[i - 10]*Data_shedding$Shedding[11] + data_2$IAV_day[i - 11]*Data_shedding$Shedding[12] + data_2$IAV_day[i - 12]*Data_shedding$Shedding[13] + data_2$IAV_day[i - 13]*Data_shedding$Shedding[14] + data_2$IAV_day[i - 14]*Data_shedding$Shedding[15] +
    data_2$IAV_day[i - 15]*Data_shedding$Shedding[16] + data_2$IAV_day[i - 16]*Data_shedding$Shedding[17] + data_2$IAV_day[i - 17]*Data_shedding$Shedding[18] + data_2$IAV_day[i - 18]*Data_shedding$Shedding[19] + data_2$IAV_day[i - 19]*Data_shedding$Shedding[20] +
    data_2$IAV_day[i - 10]*Data_shedding$Shedding[21] + data_2$IAV_day[i - 21]*Data_shedding$Shedding[22] + data_2$IAV_day[i - 22]*Data_shedding$Shedding[23] + data_2$IAV_day[i - 23]*Data_shedding$Shedding[24] + data_2$IAV_day[i - 24]*Data_shedding$Shedding[25] + data_2$IAV_day[i - 25]*Data_shedding$Shedding[26]
  WBE_2 = k * WBE 
  p_predict <- do.call(rnorm, c(10000, list(mean = log10(WBE_2), sd = sigma)))
  
  a <- quantile(p_predict, probs = c(0.025)) 
  b <- quantile(p_predict, probs = c(0.50))
  c <- quantile(p_predict, probs = c(0.975))
  
  summary_2_A <- c(summary_2_A, a)
  summary_2_B <- c(summary_2_B, b)
  summary_2_C <- c(summary_2_C, c)
}

data_result <- data.frame(A= summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("2.5%", "50%", "97.5%")
write.csv(x = data_result, file = "C:/2023_R/PRESENS_prediction_WBE_UT.csv")









##
## Determine all parametes in shedding models
Data_MA <- read.csv("20231214_PRESENS_MA.csv")
Data_NJ <- read.csv("20231214_PRESENS_NJ.csv")
Data_UT <- read.csv("20231214_PRESENS_UT.csv")


data_2 <- Data_NJ ##Data名を変更
data_2 <- select(data_2, Date, IAV_0, IAV_nov_0, IAV_case)
colnames(data_2) <- c("collection_date", "IAV", "IAV_normalized", "IAV_day")
data_2 <- mutate(data_2, IAV_normalized = IAV_normalized + 9.05) #PMMoV; 9.05
data_2 <- mutate(data_2, IAV_day = IAV_day * 38.9 / 9.01) #NJ; 9.01, Ut; 3.41, MA: 6.97


#Normalized version
data_lead <- mutate(data_2, IAV = IAV_normalized, IAV_1 = lag(IAV_normalized, 1), IAV_2 = lag(IAV_normalized, 2), IAV_3 = lag(IAV_normalized, 3), IAV_4 = lag(IAV_normalized, 4), IAV_5 = lag(IAV_normalized, 5), IAV_6 = lag(IAV_normalized, 6), IAV_7 = lag(IAV_normalized, 7), IAV_8 = lag(IAV_normalized, 8))
data_behind <- mutate(data_2, IAV_1 = lead(IAV_normalized, 1), IAV_2 = lead(IAV_normalized, 2),IAV_3 = lead(IAV_normalized, 3), IAV_4 = lead(IAV_normalized, 4), IAV_5 = lead(IAV_normalized, 5))


sample_size <- nrow(data_behind) #ここのdataファイル

data_WBE <- select(data_behind, collection_date, IAV_4) #ここのdataファイル名と数字を変える
colnames(data_WBE) <- c("collection_date", "IAV")
data_WBE <- data_WBE %>% drop_na()


#Pick up row that contains IAV concentration data
IAV_data_row <- data.frame(true = which(!is.na(data_behind$IAV_4)))#ここのdataファイル名と数字を変える
n1 <- nrow(IAV_data_row)


#artificially dataframe for shedding dynamics
shedding_time <- data.frame(values = 0:25)
n2 <- nrow(shedding_time)

data_list_CA <- list(n = sample_size, IAV = data_2$IAV_day, na = n1, WBE = data_WBE$IAV,  pick_colum = IAV_data_row$true, nb = n2, t = shedding_time$values)

mcmc_CA <- stan(
  file = "20231204_PRESENS.stan",
  data = data_list_CA,
  seed = 1,
  chain = 4,
  iter = 15000,
  warmup = 4000,
  thin = 1
)

print(mcmc_CA, pars = c("A", "a", "b", "lp__"), probe = c(0.025, 0.50, 0.975))

#WAIC
log_lik <- extract_log_lik(mcmc_CA)
waic(log_lik)
loo(log_lik)


#Traceplot
library(bayesplot)
mcmc_combo(mcmc_CA, pars = c("A", "a", "b", "sigma"))

traceplot(mcmc_CA, inc_warmup = T)


#fecal shedding
mcmc_sample <- rstan::extract(mcmc_CA)
name <- "c1" #stanの予測変数名
A <- mcmc_sample[["A"]]
a <- mcmc_sample[["a"]]
b <- mcmc_sample[["b"]]
parameter <- data.frame(A = A, a = a, b = b)


t_peak <- (log(a/b))/(a-b)
C_peak <- (a*A/(b-a))*((a/b)^(a/(b-a)))*(1-a/b)

quantile(t_peak, probs = c(0.025, 0.5, 0.975))
quantile(C_peak, probs = c(0.025, 0.5, 0.975))


##Maximum Shedding dynamics
c <- seq(0.1, 25, by = 0.1)
A <- mcmc_sample[["A"]]
a <- mcmc_sample[["a"]]
b <- mcmc_sample[["b"]]
parameter <- data.frame(A = A, a = a, b = b)
date <- data.frame(date = c)

summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

for (i in 1: nrow(date)){
  date_a <- date$date[i]
  p_fitted <- a*A/(b-a) * exp(-a*date_a)*(1-exp((a-b)*date_a))
  m <- quantile(p_fitted, probs = c(0.025)) 
  l <- quantile(p_fitted, probs = c(0.50))
  n <- quantile(p_fitted, probs = c(0.975))
  
  summary_2_A <- c(summary_2_A, m)
  summary_2_B <- c(summary_2_B, l)
  summary_2_C <- c(summary_2_C, n)
}

nrow(summary_2_C)

data_result <- data.frame(A = summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("2.5%", "50%", "97.5%")
write.csv(x = data_result, file = "C:/2023_R/sheeding_prediction.csv")


#prediction

mcmc_sample <- rstan::extract(mcmc_CA)
A <- mcmc_sample[["A"]]
a <- mcmc_sample[["a"]]
b <- mcmc_sample[["b"]]
sigma <- mcmc_sample[["sigma"]]

feces_1 <-  A * a / (b - a) * exp(-a * shedding_time$values[1]) * (1 - exp((a - b) * shedding_time$values[1]))
feces_2 <-  A * a / (b - a) * exp(-a * shedding_time$values[2]) * (1 - exp((a - b) * shedding_time$values[2]))
feces_3 <-  A * a / (b - a) * exp(-a * shedding_time$values[3]) * (1 - exp((a - b) * shedding_time$values[3]))
feces_4 <-  A * a / (b - a) * exp(-a * shedding_time$values[4]) * (1 - exp((a - b) * shedding_time$values[4]))
feces_5 <-  A * a / (b - a) * exp(-a * shedding_time$values[5]) * (1 - exp((a - b) * shedding_time$values[5]))
feces_6 <-  A * a / (b - a) * exp(-a * shedding_time$values[6]) * (1 - exp((a - b) * shedding_time$values[6]))
feces_7 <-  A * a / (b - a) * exp(-a * shedding_time$values[7]) * (1 - exp((a - b) * shedding_time$values[7]))
feces_8 <-  A * a / (b - a) * exp(-a * shedding_time$values[8]) * (1 - exp((a - b) * shedding_time$values[8]))
feces_9 <-  A * a / (b - a) * exp(-a * shedding_time$values[9]) * (1 - exp((a - b) * shedding_time$values[9]))
feces_10 <-  A * a / (b - a) * exp(-a * shedding_time$values[10]) * (1 - exp((a - b) * shedding_time$values[10]))
feces_11 <-  A * a / (b - a) * exp(-a * shedding_time$values[11]) * (1 - exp((a - b) * shedding_time$values[11]))
feces_12 <-  A * a / (b - a) * exp(-a * shedding_time$values[12]) * (1 - exp((a - b) * shedding_time$values[12]))
feces_13 <-  A * a / (b - a) * exp(-a * shedding_time$values[13]) * (1 - exp((a - b) * shedding_time$values[13]))
feces_14 <-  A * a / (b - a) * exp(-a * shedding_time$values[14]) * (1 - exp((a - b) * shedding_time$values[14]))
feces_15 <-  A * a / (b - a) * exp(-a * shedding_time$values[15]) * (1 - exp((a - b) * shedding_time$values[15]))
feces_16 <-  A * a / (b - a) * exp(-a * shedding_time$values[16]) * (1 - exp((a - b) * shedding_time$values[16]))
feces_17 <-  A * a / (b - a) * exp(-a * shedding_time$values[17]) * (1 - exp((a - b) * shedding_time$values[17]))
feces_18 <-  A * a / (b - a) * exp(-a * shedding_time$values[18]) * (1 - exp((a - b) * shedding_time$values[18]))
feces_19 <-  A * a / (b - a) * exp(-a * shedding_time$values[19]) * (1 - exp((a - b) * shedding_time$values[19]))
feces_20 <-  A * a / (b - a) * exp(-a * shedding_time$values[20]) * (1 - exp((a - b) * shedding_time$values[20]))
feces_21 <-  A * a / (b - a) * exp(-a * shedding_time$values[21]) * (1 - exp((a - b) * shedding_time$values[21]))
feces_22 <-  A * a / (b - a) * exp(-a * shedding_time$values[22]) * (1 - exp((a - b) * shedding_time$values[22]))
feces_23 <-  A * a / (b - a) * exp(-a * shedding_time$values[23]) * (1 - exp((a - b) * shedding_time$values[23]))
feces_24 <-  A * a / (b - a) * exp(-a * shedding_time$values[24]) * (1 - exp((a - b) * shedding_time$values[24]))
feces_25 <-  A * a / (b - a) * exp(-a * shedding_time$values[25]) * (1 - exp((a - b) * shedding_time$values[25]))
feces_26 <-  A * a / (b - a) * exp(-a * shedding_time$values[26]) * (1 - exp((a - b) * shedding_time$values[26]))

IAV_data_row <- data.frame(true = which(!is.na(data_lead$IAV_1)))
nrow(Data_CA)

summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

#31 10/23 3.631 Log copies/L
for(i in 30:239){
  WBE = data_2$IAV_day[i]*0 + data_2$IAV_day[i - 1]*feces_2 + data_2$IAV_day[i - 2]*feces_3 + data_2$IAV_day[i - 3]*feces_4 + data_2$IAV_day[i - 4]*feces_5 +
    data_2$IAV_day[i - 5]*feces_6 + data_2$IAV_day[i - 6]*feces_7 + data_2$IAV_day[i - 7]*feces_8 + data_2$IAV_day[i - 8]*feces_9 + data_2$IAV_day[i - 9]*feces_10 +
    data_2$IAV_day[i - 10]*feces_11 + data_2$IAV_day[i - 11]*feces_12 + data_2$IAV_day[i - 12]*feces_13 + data_2$IAV_day[i - 13]*feces_14 + data_2$IAV_day[i - 14]*feces_15 +
    data_2$IAV_day[i - 15]*feces_16 + data_2$IAV_day[i - 16]*feces_17 + data_2$IAV_day[i - 17]*feces_18 + data_2$IAV_day[i - 18]*feces_19 + data_2$IAV_day[i - 19]*feces_20 +
    data_2$IAV_day[i - 10]*feces_21 + data_2$IAV_day[i - 21]*feces_22 + data_2$IAV_day[i - 22]*feces_23 + data_2$IAV_day[i - 23]*feces_24 + data_2$IAV_day[i - 24]*feces_25 + data_2$IAV_day[i - 25]*feces_26
  
  p_predict <- do.call(rnorm, c(10000, list(mean = log10(WBE), sd = sigma)))
  
  a <- quantile(p_predict, probs = c(0.025)) 
  b <- quantile(p_predict, probs = c(0.50))
  c <- quantile(p_predict, probs = c(0.975))
  
  summary_2_A <- c(summary_2_A, a)
  summary_2_B <- c(summary_2_B, b)
  summary_2_C <- c(summary_2_C, c)
}

data_result <- data.frame(A = summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("2.5%", "50%", "97.5%")
write.csv(x = data_result, file = "C:/2023_R/PRESENS_prediction_WBE_0_final_day.csv")


