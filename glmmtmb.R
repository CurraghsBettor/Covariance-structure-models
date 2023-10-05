rm(list=ls())

pkg <- c("glmmTMB", "MASS", "reshape2", "ez", "Hmisc", "car")
for(package in pkg) {
  library(package, character.only = T)
}

set.seed(21450)

## parameters for the bivariate Normal distribution
# Number of cells per group
n <- 30
# mu
mu <- c(0, 1, 2, 3, 3.8)
mu2 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
# Variance-Covariance matrix
Sigma <- matrix(0.8, nrow = 5, ncol = 5); diag(Sigma) <- 1
Sigma2 <- matrix(0.1, nrow = 5, ncol = 5); diag(Sigma2) <- 1
Sigma3 <- matrix(c(1, .8, .6, .4, .2,
                 .8, 1, .8, .6, .4,
                 .6, .8, 1, .8, .6,
                 .4, .6, .8, 1, .8,
                 .2, .4, .6, .8, 1), nrow = 5, ncol = 5); print(Sigma3)

## Multivariate Normal distributions
Data <- mvrnorm(n = n, mu = mu, Sigma = Sigma)
Datab <- mvrnorm(n = n, mu = mu2, Sigma = Sigma)
Data2 <- mvrnorm(n = n, mu = mu, Sigma = Sigma2)
Data2b <- mvrnorm(n = n, mu = mu2, Sigma = Sigma2)
Data3 <- mvrnorm(n = n, mu = mu, Sigma = Sigma3)
Data3b <- mvrnorm(n = n, mu = mu2, Sigma = Sigma3)

## Pearson
cor <- sapply(list(Data, Data2, Data3), 
              function(data) rcorr(data))
cor[, 1]
cor[, 2]
cor[, 3]

# Merge
Data <- rbind(Data, Datab)
Data2 <- rbind(Data2, Data2b)
Data3 <- rbind(Data3, Data3b)

rm(Datab, Data2b, Data3b)

# Number of observations 
ID <- 1:sum(n,n)

# merge              
Data <- cbind(ID, Data)
Data2 <- cbind(ID, Data2)
Data3 <- cbind(ID, Data3)

Data <- as.data.frame(Data)
Data2 <- as.data.frame(Data2)
Data3 <- as.data.frame(Data3)

Data$Group <- c(rep(1, n), rep(2, n))
Data2$Group <- c(rep(1, n), rep(2, n))
Data3$Group <- c(rep(1, n), rep(2, n))

colnames(Data) <- c("ID", "Time1", "Time2", "Time3", "Time4", "Time5", "Group")
colnames(Data2) <- c("ID", "Time1", "Time2", "Time3", "Time4", "Time5", "Group")
colnames(Data3) <- c("ID", "Time1", "Time2", "Time3", "Time4", "Time5", "Group")

Data <- melt(data = Data, id.vars = c("ID", "Group"), 
             measure.vars = c("Time1", "Time2", "Time3", "Time4", "Time5"),
             variable.name = "Time", value.name = "DV")

Data2 <- melt(data = Data2, id.vars = c("ID", "Group"), 
             measure.vars = c("Time1", "Time2", "Time3", "Time4", "Time5"),
             variable.name = "Time", value.name = "DV")

Data3 <- melt(data = Data3, id.vars = c("ID", "Group"), 
              measure.vars = c("Time1", "Time2", "Time3", "Time4", "Time5"),
              variable.name = "Time", value.name = "DV")
## Factors
Data$Time <- as.factor(Data$Time)
Data$Group <- as.factor(Data$Group)
Data$ID <- as.factor(Data$ID)

Data2$Time <- as.factor(Data2$Time)
Data2$Group <- as.factor(Data2$Group)
Data2$ID <- as.factor(Data2$ID)

Data3$Time <- as.factor(Data3$Time)
Data3$Group <- as.factor(Data3$Group)
Data3$ID <- as.factor(Data3$ID)

qqPlot(Data$DV)
ezPlot(data = Data, dv =. (DV), wid =.(ID),
       within =. (Time), between =. (Group), x=. (Time), 
       split =. (Group), do_lines = T)

qqPlot(Data2$DV)
ezPlot(data = Data2, dv =. (DV), wid =.(ID),
       within =. (Time), between =. (Group), x=. (Time), 
       split =. (Group), do_lines = T)

qqPlot(Data3$DV)
ezPlot(data = Data3, dv =. (DV), wid =.(ID),
       within =. (Time), between =. (Group), x=. (Time),
       split =. (Group), do_lines = T)

### Data 1
## ML
# covariance Structure: Unstructured 
Model1 <- glmmTMB(DV ~ Group + Time + Group:Time + us(0 + Time|ID), data = Data,
                  control = glmmTMBControl(optimizer=nlminb),
                  dispformula = ~1, REML = F)
# Heterogeneous Compound Symmetry
Model2 <- glmmTMB(DV ~ Group + Time + Group:Time + cs(0 + Time|ID), data = Data,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")), 
                  dispformula = ~1, REML = F)
# Variance Components/Diagonal Heterogeneous Variance
Model3 <- glmmTMB(DV ~ Group + Time + Group:Time + diag(0 + Time|ID), data = Data,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                  dispformula = ~1, REML = F)
# Variance Components/Diagonal Homogeneous Variance
Model4 <- glmmTMB(DV ~ Group + Time + Group:Time + homdiag(0 + Time|ID), data = Data,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                  dispformula = ~1, REML = F)

AIC(Model1, Model2, Model3, Model4)
BIC(Model1, Model2, Model3, Model4)

VarC <- sapply(list(Model1, Model2, Model3, Model4),
               function(model) VarCorr(model))

Vcov <- sapply(list(Model1, Model2, Model3, Model4),
            function(model) vcov(model))
VarC[[1]]
Vcov[[1]]

rm(Model1, Model2, Model3, Model4, VarC, Vcov)

# fit the best model with REML
Model1b <- glmmTMB(DV ~ Group + Time + Group:Time + us(0 + Time|ID), data = Data,
                   control = glmmTMBControl(optimizer=nlminb), 
                   dispformula = ~1, REML = T)
summary(Model1b)
Anova(Model1b, type = 'III', test.statistic = "Chisq")

Model1c <- glmmTMB(DV ~ Group + Time + Group:Time + us(0 + Time|ID), data = Data,
                   control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                   dispformula = ~I(Time), REML = T) ## does not converge
# summary(Model1c)
# Anova(Model1c, type = 'III', test.statistic = "Chisq")

rm(Model1b, Model1c)

### Data 2
## ML
Model1 <- glmmTMB(DV ~ Group + Time + Group:Time + us(0 + Time|ID), data = Data2,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")),
                  dispformula = ~1, REML = F)
Model2 <- glmmTMB(DV ~ Group + Time + Group:Time + cs(0 + Time|ID), data = Data2,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")), 
                  dispformula = ~1, REML = F)
Model3 <- glmmTMB(DV ~ Group + Time + Group:Time + diag(0 + Time|ID), data = Data2,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                  dispformula = ~1, REML = F)
Model4 <- glmmTMB(DV ~ Group + Time + Group:Time + homdiag(0 + Time|ID), data = Data2,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                  dispformula = ~1, REML = F)

AIC(Model1, Model2, Model3, Model4)
BIC(Model1, Model2, Model3, Model4)

VarC <- sapply(list(Model1, Model2, Model3, Model4),
               function(model) VarCorr(model))

Vcov <- sapply(list(Model1, Model2, Model3, Model4),
            function(model) vcov(model))

# fit the best model with REML
Model2b <- glmmTMB(DV ~ Group + Time + Group:Time + cs(0 + Time|ID), data = Data2,
                   control = glmmTMBControl(optimizer=optim, optArgs=list(method="L-BFGS-B")), 
                   dispformula = ~1, REML = T)
summary(Model2b)
Anova(Model2b, type = 'III', test.statistic = "Chisq")

Model2c <- glmmTMB(DV ~ Group + Time + Group:Time + cs(0 + Time|ID), data = Data2,
                   control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                   dispformula = ~I(Time), REML = T)
summary(Model2c)
Anova(Model2c, type = 'III', test.statistic = "Chisq")

rm(Model1, Model2, Model3, Model4, VarC, Vcov)

Model1 <- glmmTMB(DV ~ Group + Time + Group:Time + us(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer=nlminb),
                  dispformula = ~1, REML = F)
Model2 <- glmmTMB(DV ~ Group + Time + Group:Time + cs(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")), 
                  dispformula = ~1, REML = F)
Model3 <- glmmTMB(DV ~ Group + Time + Group:Time + diag(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="CG")), 
                  dispformula = ~1, REML = F)
Model4 <- glmmTMB(DV ~ Group + Time + Group:Time + homdiag(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")), 
                  dispformula = ~1, REML = F)
# First-Order Auto-Regressive AR1
Model5 <- glmmTMB(DV ~ Group + Time + Group:Time + ar1(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")),
                  dispformula = ~1, REML = F)
# Heterogeneous First-Order Auto-Regressive ARH1
Model6 <- glmmTMB(DV ~ Group + Time + Group:Time + ar1(0 + Time|ID) + diag(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer = optim, optArgs = list(method="CG")),
                  dispformula = ~1, REML = F)

AIC(Model1, Model2, Model3, Model4, Model5, Model6)
BIC(Model1, Model2, Model3, Model4, Model5, Model6)

VarC <- sapply(list(Model1, Model2, Model3, Model4, Model5, Model6),
               function(model) VarCorr(model))

# the worst model: Model3
VarC[[3]]
# the best model: Model5
VarC[[5]]

Vcov <- sapply(list(Model1, Model2, Model3, Model4, Model5, Model6),
            function(model) vcov(model))

# the worst model: Model3
Vcov[[3]]
# the best model: Model5
Vcov[[5]]

rm(Model1, Model2, Model3, Model4, Model5, Model6, Vcov, VarC)

## fit the best model with REML
Model5b <- glmmTMB(DV ~ Group + Time + Group:Time + ar1(0 + Time|ID), data = Data3,
                  control = glmmTMBControl(optimizer = optim, optArgs = list(method="BFGS")),
                  dispformula = ~1, REML = T)
summary(Model5b)
Anova(Model5b, type = "III")
