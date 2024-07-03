Data <- read.csv("20231129_IAV_logistic.csv")

library(plyr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(rstan)
library(brms)

data <- select(Data, Probability, IAV_2, Sample)
colnames(data) <- c("Probability", "IAV", "Sample")
data <- data %>% drop_na()

mean(data$Sample)

sample_size <- nrow(data)
p <- data$Probability/100
IAV <- log10(data$IAV/389)

data_list_CA <- list(p = p, IAV = IAV, n = sample_size, nos = data$Sample)

mcmc_CA <- stan(
  file = "20231203_logit.stan",#, 20231129_IAV_logistic.stan
  data = data_list_CA,
  seed = 1,
  chain = 4,
  iter = 10000,
  warmup = 2000,
  thin = 1
)

print(mcmc_CA, probe = c(0.025, 0.50, 0.975))

log_lik <- extract_log_lik(mcmc_CA)
waic(log_lik)

loo(log_lik)

mcmc_sample <- rstan::extract(mcmc_CA)
name <- "c1" #stanの予測変数名
c1 <- mcmc_sample[["c1"]]
intercept <- mcmc_sample[["intercept"]]
sigma <- mcmc_sample[["sigma"]]
parameter <- data.frame(c1 = c1, intercept = intercept, sigma = sigma)

parameter <- data.frame(c1 = c1, intercept = intercept)

summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

for (i in 1: nrow(data)){
  IAV <- data$IAV[i]/389
  p_fitted <- exp(intercept + c1*log10(IAV)) / (1 + exp(intercept + c1*log10(IAV)))
  p_predict <- do.call(rnorm, c(10000, list(mean = p_fitted, sd = sigma)))
  a <- quantile(p_predict, probs = c(0.025)) 
  b <- quantile(p_predict, probs = c(0.50))
  c <- quantile(p_predict, probs = c(0.975))

  summary_2_A <- c(summary_2_A, a)
  summary_2_B <- c(summary_2_B, b)
  summary_2_C <- c(summary_2_C, c)
  }

data_result <- data.frame(IAV_Case = data$IAV, A = summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("IAV_cases", "2.5%", "50%", "97.5%")
write.csv(x = data_result, file = "C:/2023_R/logistic_regression_prediction.csv")

#Back calculation: 1 per 10000 cases
summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

p_fitted <- exp(intercept + c1*log10(1)) / (1 + exp(intercept + c1*log10(1)))
p_predict <- do.call(rnorm, c(10000, list(mean = p_fitted, sd = sigma)))
a <- quantile(p_predict, probs = c(0.025)) 
b <- quantile(p_predict, probs = c(0.50))
c <- quantile(p_predict, probs = c(0.975))

summary_2_A <- c(summary_2_A, a)
summary_2_B <- c(summary_2_B, b)
summary_2_C <- c(summary_2_C, c)


data_result <- data.frame(IAV_Case = data$IAV, A = summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("IAV_cases", "2.5%", "50%", "97.5%")



##Bionmial_fitting(Approximation, stan_file_name: 20231203_logit.stan)
mcmc_sample <- rstan::extract(mcmc_CA)

c1 <- mcmc_sample[["c1"]]
intercept <- mcmc_sample[["intercept"]]

summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

for (i in 1: nrow(data)){
  IAV <- data$IAV[i]/389
  p_fitted <- exp(intercept + c1*log10(IAV)) / (1 + exp(intercept + c1*log10(IAV)))
  p_predict <- do.call(rnorm, c(10000, list(mean = p_fitted, sd = sqrt(p_fitted*(1-p_fitted)/data$Sample[i]))))
  a <- quantile(p_predict, probs = c(0.025)) 
  b <- quantile(p_predict, probs = c(0.50))
  c <- quantile(p_predict, probs = c(0.975))
  
  summary_2_A <- c(summary_2_A, a)
  summary_2_B <- c(summary_2_B, b)
  summary_2_C <- c(summary_2_C, c)
}

data_result <- data.frame(IAV_Case = data$IAV, Data = data$Probability, A = summary_2_A*100, B = summary_2_B*100, C = summary_2_C*100)
colnames(data_result) <- c("IAV_cases", "Data", "Min", "Median", "Max")

#California population: 38.9 million ここは修正する必要あり。
plot <- ggplot(data_result, aes(x = log10(IAV_cases/38.9*0.1))) +
  geom_ribbon(aes(ymin = Min, ymax = Max), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = Median), color = "blue", size = 1) +
  geom_point(aes(y = Data), color = "blue", size = 3) +
  labs(title = "IAVに対するグラフ", x = "IAV", y = "Values")
plot <- plot + theme_classic()
plot <- plot + theme(
  axis.line = element_line(linewidth = 1.0, lineend = "square"),
  text = element_text(colour ="black", size = 14),
  legend.position = "none",
  axis.ticks = element_line(linewidth = 1.5),
  axis.ticks.length = unit(-2, "mm"))
plot


#Back calculation: 1 per 10000 cases
summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

  p_fitted <- exp(intercept + c1*log10(1)) / (1 + exp(intercept + c1*log10(1)))
  p_predict <- do.call(rnorm, c(10000, list(mean = p_fitted, sd = sqrt(p_fitted*(1-p_fitted)/100))))
  a <- quantile(p_predict, probs = c(0.025)) 
  b <- quantile(p_predict, probs = c(0.50))
  c <- quantile(p_predict, probs = c(0.975))
  
  summary_2_A <- c(summary_2_A, a)
  summary_2_B <- c(summary_2_B, b)
  summary_2_C <- c(summary_2_C, c)


data_result <- data.frame(IAV_Case = data$IAV, A = summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("IAV_cases", "2.5%", "50%", "97.5%")

write.csv(x = data_result, file = "C:/2023_R/logistic_regression_prediction.csv")


#Back calculation: estimation
summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL

p_fitted <- exp(intercept + c1*log10(0.78)) / (1 + exp(intercept + c1*log10(0.78)))
p_predict <- do.call(rnorm, c(10000, list(mean = p_fitted, sd = sqrt(p_fitted*(1-p_fitted)/100))))
quantile(p_predict, probs = c(0.025, 0.50, 0.975)) 

Median <- 10^(-median(intercept)/(median(c1)))





##memo
sample_size_2 <- nrow(parameter)
summary_2_A <- NULL
summary_2_B <- NULL
summary_2_C <- NULL
for (i in 1:sample_size){
  IAV <- log10(data$IAV[i])
  p_estimate_2 <- NULL
  
  for(m in 1:sample_size_2){
    c1 <- parameter$c1[m]
    intercept <- parameter$intercept[m]
    p_estimate <- exp(intercept + c1*IAV)/(1 + exp(intercept + c1*IAV))
    p_estimate_2 <- c(p_estimate_2, p_estimate)
  }
  summary <- quantile(p_estimate_2, probs = c(0.025, 0.5, 0.975))
  summary_A <- quantile(p_estimate_2, probs = c(0.025))
  summary_B <- quantile(p_estimate_2, probs = c(0.5))
  summary_C <- quantile(p_estimate_2, probs = c(0.975))
  
  summary_2_A <- c(summary_2_A, summary_A)
  summary_2_B <- c(summary_2_B, summary_B)
  summary_2_C <- c(summary_2_C, summary_C)
}
data_result <- data.frame(A = summary_2_A, B = summary_2_B, C = summary_2_C)
colnames(data_result) <- c("2.5%", "50%", "97.5%")
data_final <- cbind(data, data_result)
write.csv(x = data_final, file = "C:/2023_R/regression_model.csv")


#prediction
