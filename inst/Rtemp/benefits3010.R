
## -----------------------------------------------------------------------
# Calculating the benefit of continuing endocrine therapy
# assuming survival to 5 years
start <- 6
m.oth <- m.oth[start:15]
s.cum.oth <- 1-cumsum(m.oth)

# Generate the annual breast cancer specific mortality rate
m.br <- base.m.br[start:15]*exp(rx[start:15,])

# Calculate the cumulative breast cancer mortality rate
m.cum.br <- apply(m.br, 2, cumsum)

# Calculate the cumulative breast cancer survival
s.cum.br <- exp(- m.cum.br)
m.cum.br <- 1- s.cum.br

m.br <- m.cum.br
for (j in 1:cols) {
  for (i in 2:10) {
    m.br[i,j] <- m.cum.br[i,j] - m.cum.br[i-1,j]
  }
}

# Cumulative all cause mortality conditional on surviving breast and all cause mortality
m.cum.all <- 1 - s.cum.oth*s.cum.br
s.cum.all <- 100-100*m.cum.all

# Annual all cause mortality
m.all <- m.cum.all
for (j in 1:cols) {
  for (i in 2:10) {
    m.all[i,j] <- m.cum.all[i,j] - m.cum.all[i-1,j]
  }
}

# Proportion of all cause mortality that is breast cancer
prop.br      <- m.br/(m.br + m.oth)
pred.m.br    <- prop.br*m.all
pred.cum.br  <- apply(pred.m.br, 2, cumsum)
pred.m.oth   <- m.all - pred.m.br
pred.cum.oth <- apply(pred.m.oth, 2, cumsum)
pred.cum.all <- pred.cum.br + pred.cum.oth

surv_conditional <- 100*(1-pred.cum.all)[c(10), c(1:3)]

# rm(list=setdiff(ls(), c("benefits", "surv", "surv_conditional")))

