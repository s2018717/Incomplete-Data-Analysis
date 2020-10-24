# 3(a) a = 2 and b = 0

require(MASS)
n = 500

#defining the covariance matrix
Sigma = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3, byrow = T)

#simulating the data
set.seed(123)
Z = mvrnorm(n = n, mu = c(0,0,0), Sigma = Sigma)
#storing each column of Z in a separate variable
Z1 = Z[,1]; Z2 <- Z[,2]; Z3 <- Z[,3]
Y1 = 1 + Z1
Y2 = 5 + 2*Z1 + Z2

#extracting the index that meets the condition
ind_mar = which(2*(Y1-1) + Z3 < 0)
#storing the observed values in new variables
Y2_MAR_obs = Y2[-ind_mar]

#plotting the densities
plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y[2]),
     main = "MAR", ylim = c(0, 0.3))
lines(density(Y2_MAR_obs), lwd = 2, col = "red")
legend(7.5, 0.3, legend = c("Complete data", "Observed data"),col = c("blue", "red"),
       lty = c(1,1), lwd = c(2,2), bty ="n")


# 3(b) stochastic regression imputation

Y2_MAR = Y2
Y2_MAR[ind_mar] = NA
df_Y = data.frame(Y1, Y2, Y2_MAR)
fit = lm(Y2_MAR ~ Y1)
summary(fit)
set.seed(1234)
predicted_sri = predict(fit, newdata = df_Y) + rnorm(nrow(df_Y), 0, sigma(fit))
Y_sri = ifelse(is.na(df_Y$Y2_MAR), predicted_sri, df_Y$Y2_MAR)

plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y[2]),
     main = "Stochastic Regression Imputation", ylim = c(0, 0.25))
lines(density(Y_sri), lwd = 2, col = "red")
legend(2, 0.25, legend = c("Complete data", "Completed (after imputation) data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")

# 3(c) a = 0 and b = 2

#extracting the index that meets the condition
ind_mnar = which(2*(Y2-5) + Z3 < 0)
#storing the observed values in new variables
Y2_MNAR_obs = Y2[-ind_mnar]

#plotting the densities
plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y[2]),
     main = "MNAR", ylim = c(0, 0.3))
lines(density(Y2_MNAR_obs), lwd = 2, col = "red")
legend(7.5, 0.3, legend = c("Complete data", "Observed data"),col = c("blue", "red"),
       lty = c(1,1), lwd = c(2,2), bty ="n")

# 3(d) stochastic regression imputation
Y2_MNAR = Y2
Y2_MNAR[ind_mnar] = NA
df_Y = data.frame(Y1, Y2, Y2_MAR, Y2_MNAR)
fit = lm(Y2_MNAR ~ Y1)
summary(fit)
set.seed(111)
predicted_sri2 = predict(fit, newdata = df_Y) + rnorm(nrow(df_Y), 0, sigma(fit))
Y_sri2 = ifelse(is.na(df_Y$Y2_MNAR), predicted_sri2, df_Y$Y2_MNAR)

plot(density(Y2), lwd = 2, col = "blue", xlab = expression(Y[2]),
     main = "Stochastic Regression Imputation", ylim = c(0, 0.35))
lines(density(Y_sri2), lwd = 2, col = "red")
legend(2, 0.35, legend = c("Complete data", "Completed (after imputation) data"),
       col = c("blue", "red"), lty = c(1,1), lwd = c(2,2), bty ="n")


# 4(a) complete case analysis

load("databp.Rdata")
ind = which(is.na(databp$recovtime) == FALSE)
mccoverall = mean(databp$recovtime, na.rm = TRUE)
seccoverall = sd(databp$recovtime, na.rm = TRUE)/sqrt(length(ind))
mccoverall; seccoverall
corrd = cor(databp$recovtime[ind], databp$logdose[ind])
corrb = cor(databp$recovtime[ind], databp$bloodp[ind])
corrd; corrb

# 4(b) mean imputation

recmi = ifelse(is.na(databp$recovtime) == TRUE, mean(databp$recovtime, na.rm = TRUE), databp$recovtime)
mmi = mean(recmi)
n = nrow(databp)
semi = sd(recmi)/sqrt(n)
mmi; semi
corrdmi = cor(recmi, databp$logdose)
corrbmi = cor(recmi, databp$bloodp)
corrdmi; corrbmi

# 4(c) mean regression imputation

fit = lm(recovtime ~ logdose + bloodp, data = databp)
summary(fit)
predicted_ri = predict(fit, newdata = databp)
recri = ifelse(is.na(databp$recovtime), predicted_ri, databp$recovtime)
mri = mean(recri)
sdri = sd(recri)/sqrt(n)
mri; sdri
corrdri = cor(recri, databp$logdose)
corrbri = cor(recri, databp$bloodp)
corrdri; corrbri

# 4(d) stochastic regression imputation

set.seed(123)
predicted_sri = predict(fit, newdata = databp) + rnorm(nrow(databp), 0, sigma(fit))
recsri = ifelse(is.na(databp$recovtime), predicted_sri, databp$recovtime)
msri = mean(recsri)
sesri = sd(recsri)/sqrt(n)
msri; sesri
corrdsri = cor(recsri, databp$logdose)
corrbsri = cor(recsri, databp$bloodp)
corrdsri; corrbsri

# 4(e) predictive mean matching

ind_mis = which(is.na(databp$recovtime))
fit = lm(recovtime ~ logdose + bloodp, data = databp)
summary(fit)
predicted = predict(fit, newdata = databp)
predicted

# calculate squared difference

sd_4 = (predicted[-ind_mis]-predicted[4])^2 
sd_10 = (predicted[-ind_mis]-predicted[10])^2 
sd_22 = (predicted[-ind_mis]-predicted[22])^2 

# find the index of the donor
which.min(sd_4)
which.min(sd_10)
which.min(sd_22)

# replace every missing value
rechd = databp$recovtime
rechd[4] = predicted[6]
rechd[10] = predicted[2]
rechd[22] = predicted[17]

mhd <- mean(rechd); sehd <- sd(rechd)/sqrt(n)
mhd; sehd
corrdhd = cor(rechd, databp$logdose)
corrbhd = cor(rechd, databp$bloodp)
corrdhd; corrbhd