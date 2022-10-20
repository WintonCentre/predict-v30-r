# Version as implemented in PREDICT v2.0 on web
# Error in estimating benefit of treatment fixed (line 178)
#   benefit measure as diff in all cause mortality and 
#   not breast cnacer specific mortality
# Addition of survival to 15 years
# Addition of additional hormone therapy for 5 year survivors
# Addition of bisphosphonates
# Addition of radiotherapyhelp rm
# Additions of DFS

library(tidyverse)

#ydirectory <- "/Users/paul/Google Drive/Breast/Predict/Predict-v2"
#setwd(directory)

# Set input paramenters
age.start  <- 65
screen     <- 0     # Clinically detected = 0, screen detected = 1
size       <- 19    # Tumour size mm
grade      <- 2     # Tumour grade
nodes      <- 1     # Number positive nodes
er         <- 1     # ER+ = 1, ER- = 0
her2       <- 0     # HER2+ = 1, HER2- = 0, missing = 9
ki67       <- 0     # KI67+ = 1, KI67- = 0, missing = 9
generation <- 0     # Chemo generation 0, 2 or 3 only
horm       <- 1     # Hormone therapy Yes = 1, no = 0
traz       <- 0     # Trastuzumab therapy Yes = 1, no = 0
bis        <- 0     # Bisphosphonate therapy Yes = 1, no = 0
radio      <- 1     # Radiotherapy Yes = 1, no = 0

# Derived variables
time      <- c(1:15)
age       <- age.start - 1 + time
grade.val <- ifelse(er==1, grade, ifelse(grade==1, 0, 1)) # Grade variable for ER neg

# Generate coefficients
age.mfp.1   <- ifelse(er==1, (age.start/10)^-2-.0287449295, age.start-56.3254902)
age.beta.1  <- ifelse(er==1, 34.53642, 0.0089827)
age.mfp.2   <- ifelse(er==1, (age.start/10)^-2*log(age.start/10)-.0510121013, 0)
age.beta.2  <- ifelse(er==1, -34.20342, 0)
size.mfp    <- ifelse(er==1, log(size/100)+1.545233938, (size/100)^.5-.5090456276)
size.beta   <- ifelse(er==1, 0.7530729, 2.093446)
nodes.mfp   <- ifelse(er==1, log((nodes+1)/10)+1.387566896, log((nodes+1)/10)+1.086916249)
nodes.beta  <- ifelse(er==1, 0.7060723, .6260541)
grade.beta  <- ifelse(er==1, 0.746655, 1.129091)
screen.beta <- ifelse(er==1, -0.22763366, 0)
her2.beta   <- ifelse(her2==1, 0.2413,
                    ifelse(her2==0, -0.0762,0 ))
ki67.beta   <- ifelse(ki67==1 & er==1, 0.14904,
                    ifelse(ki67==0 & er==1, -0.1133,0 ))
r.base.br   <- 0.133   # To set baseline of no radioRX breast mortality HR=0.82 and
r.base.oth  <- -0.047  # mortality HR=1.07 and 65% getting radiotherapy
dfs         <- 0.223   # Relative risk relapse compared to breast mortality

# Calculate the other mortality prognostic index (mi)
mi <- 0.0698252*((age.start/10)^2-34.23391957) + r.base.oth

# Calculate the breast cancer mortality prognostic index (pi)
pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
  nodes.beta*nodes.mfp + grade.beta*grade.val + screen.beta*screen + 
  her2.beta + ki67.beta + r.base.br

# Generate the therapy related coefficients
c     <- ifelse(generation == 0, 0, ifelse(generation == 2, -.248, -.446))
h     <- ifelse(horm==1 & er==1, -0.3857, 0)
t     <- ifelse(her2==1 & traz==1, -.3567, 0)
b     <- ifelse(bis==1, -0.198, 0) # Only applicable to menopausal women.
#o    <- ifelse(ov==1, -0.2, 0) # Ready to add oophorectomy option
r.br  <- ifelse(radio==1, -0.198, 0) 
r.oth <- ifelse(radio==1, 0.068, 0) 
hr    <- h + r.br
hrc   <- h + r.br + c
hrct  <- h + r.br + c + t
hrctb <- h + r.br + c + t + b

# Generate cumulative baseline other mortality
base.m.cum.oth <- exp(-6.052919 + (1.079863*log(time)) + (.3255321*time^.5))

# Generate cumulative survival non-breast mortality
s.cum.oth <- exp(-exp(mi)*base.m.cum.oth)

# Generate annual baseline breast mortality
base.m.oth <- base.m.cum.oth
for (i in 2:15) {
  base.m.oth[i] <- base.m.cum.oth[i] - base.m.cum.oth[i-1] }

# Loop for different treatment options
rx.oth <- c(b = 0, h = 0, r = r.oth, c = 0, 
        t = 0, hr = r.oth, hrc = r.oth, hrct = r.oth, hrctb = r.oth)
cols <- length(rx.oth) # Number of RX categories

# Generate the annual breast cancer specific mortality rate
m.oth.rx <- sapply(rx.oth, function(rx.oth, x.vector = base.m.oth) {
  output <-  x.vector*exp(mi + rx.oth)
  return(output)
}
)

# Generate annual survival from cumulative survival
m.cum.oth <- 1 - s.cum.oth 
m.oth <- 1 - s.cum.oth 
for (i in 2:15) {
  m.oth[i] <- m.cum.oth[i] - m.cum.oth[i-1]
}

# Calculate the cumulative other mortality rate
m.cum.oth.rx <- apply(m.oth.rx, 2, cumsum) 

# Calculate the cumulative breast cancer survival
s.cum.oth.rx <- exp(-m.cum.oth.rx)  

# Convert cumulative mortality rate into cumulative risk
m.cum.oth.rx <- 1- s.cum.oth.rx

m.oth.rx <- m.cum.oth.rx
for (j in 1:cols) {
  for (i in 2:15) {
    m.oth.rx[i,j] <- m.cum.oth.rx[i,j] - m.cum.oth.rx[i-1,j]
  }
}

# Generate cumulative baseline breast mortality
if (er==1) {
  base.m.cum.br <- exp(0.7424402 - 7.527762/time^.5 - 1.812513*log(time)/time^.5) 
} else { base.m.cum.br <- exp(-1.156036 + 0.4707332/time^2 - 3.51355/time)
}

# Generate annual baseline breast mortality
base.m.br <- base.m.cum.br
for (i in 2:15) {
  base.m.br[i] <- base.m.cum.br[i] - base.m.cum.br[i-1] }

# Loop for different treatment options
rx <- c(b = 0, h = h, r = r.br, c = c, 
           t = t, hr = hr, hrc = hrc, hrct = hrct, hrctb = hrctb)
cols <- length(rx) # Number of RX categories

# Generate the annual breast cancer specific mortality rate
m.br.rx <- sapply(rx, function(rx, x.vector = base.m.br) {
  output <-  x.vector*exp(pi + rx)
  return(output)
}
)

# Calculate the cumulative breast cancer mortality rate
m.cum.br.rx <- apply(m.br.rx, 2, cumsum) 

# Calculate the cumulative breast cancer survival
s.cum.br.rx <- exp(- m.cum.br.rx)  

# Convert cumulative mortality rate into cumulative risk
m.cum.br.rx <- 1- s.cum.br.rx

m.br.rx <- m.cum.br.rx
for (j in 1:cols) {
  for (i in 2:15) {
    m.br.rx[i,j] <- m.cum.br.rx[i,j] - m.cum.br.rx[i-1,j]
  }
}

# Cumulative all cause mortality conditional on surviving breast and all cause mortality
m.cum.all.rx <- 1 - s.cum.oth.rx*s.cum.br.rx
s.cum.all.rx <- 100-100*m.cum.all.rx

# Annual all cause mortality
m.all.rx <- m.cum.all.rx
for (j in 1:cols) {
  for (i in 2:15) {
    m.all.rx[i,j] <- m.cum.all.rx[i,j] - m.cum.all.rx[i-1,j]
  }
}

# Proportion of all cause mortality that is breast cancer
prop.br.rx      <- m.br.rx/(m.br.rx + m.oth)
pred.m.br.rx    <- prop.br.rx*m.all.rx
pred.cum.br.rx  <- apply(pred.m.br.rx, 2, cumsum) 
pred.m.oth.rx   <- m.all.rx - pred.m.br.rx
pred.cum.oth.rx <- apply(pred.m.oth.rx, 2, cumsum)
pred.cum.all.rx <- pred.cum.br.rx + pred.cum.oth.rx

# rx benefits
benefits2.1 <- 100*(pred.cum.all.rx[,1] - pred.cum.all.rx)

################################
# DFS as the outcome measure   #
################################

# RH distance recurrence in first five years based on EBCTCC (2011) 
#   where there are multiple recurrence/BCSM absolute risk comparisons.
#   Estimated mean HR for 0-4, 5-14 years
#   as function of baseline breast mortality rate
dfs <- c(rep(0.617,5), rep(-0.0603,10))
base.rel <- base.m.br*exp(dfs)

# Generate the baseline annual relapse rate
base.rel.rx <- sapply(rx, function(rx, x.vector = base.rel) {
  output <-  x.vector*exp(pi + rx)
  return(output)
}
)

# Calculate the cumulative relapse rate
rel.cum.rx <- apply(base.rel.rx, 2, cumsum) 

# Calculate the cumulative relapse survival
s.cum.rel.rx <- exp(- rel.cum.rx)

# Convert cumulative relapse rate into cumulative risk
rel.cum.rx <- 1- s.cum.rel.rx

# Calculate annual relapse risk
rel.rx <- rel.cum.rx
for (j in 1:cols) {
  for (i in 2:15) {
    rel.rx[i,j] <- rel.cum.rx[i,j] - rel.cum.rx[i-1,j]
  }
}

# Cumulative all cause survival conditional on surviving relapse and all cause mortality
m.cum.all.rx <- 1 - s.cum.oth.rx*s.cum.rel.rx
s.cum.all.rx <- 100-100*m.cum.all.rx

# Annual all cause event rate
m.all.rx <- m.cum.all.rx
for (j in 1:cols) {
  for (i in 2:15) {
    m.all.rx[i,j] <- m.cum.all.rx[i,j] - m.cum.all.rx[i-1,j]
  }
}

# Proportion of all events that is relapse
prop.rel.rx <- rel.rx/(rel.rx + m.oth)

# Calulate the predicted relapse and other mortality
pred.m.rel.rx       <- prop.rel.rx*m.all.rx
pred.cum.rel.rx     <- apply(pred.m.rel.rx, 2, cumsum) 
pred.m.oth.rx       <- m.all.rx - pred.m.rel.rx
pred.cum.oth.rx     <- apply(pred.m.oth.rx, 2, cumsum)
pred.cum.all.rel.rx <- pred.cum.rel.rx + pred.cum.oth.rx

# rx benefits

benefits2.1.2 <- 100*(pred.cum.all.rx[,1] - pred.cum.all.rx)

################################
# 10 years hormone from year 5 #
################################

# Conditional on survival to five years (entered as delay)
# Only variable is duration hormone therapy
# Other treatment options fixed
delay <- 5   # Set delay to 0 to check predictions with main model
rows  <- 15 - delay
time  <- c(1:rows)
age   <- age.start + delay - 1 + time

# Generate cumulative survival non-breast mortality
s.cum.oth10 <- exp(-exp(mi + r.oth)*base.m.cum.oth)

# Generate annual survival from cumulative survival
m.oth.10 <- m.oth[(1+delay):15]
m.cum.oth.10 <- cumsum(m.oth.10)
s.cum.oth.10 <- 1 - m.cum.oth.10

# Generate the time specific treatment coefficients
h5   <- c(rep(h,rows))
h10  <- c(rep(h, rows - 5), rep(-.2+h, 5))
pi5  <- pi + h5 + r.br + c + t + b
pi10 <- pi + h10 + r.br + c + t + b
rx10 <- cbind(pi5, pi10)

# Generate the annual breast cancer specific mortality rate
m.br.rx.10 <- base.m.br[(1+delay):15]*exp(rx10)

# Calculate the cumulative breast cancer mortality rate
m.cum.br.rx.10 <- apply(m.br.rx.10, 2, cumsum)

# Calculate the cumulative breast cancer survival
s.cum.br.rx.10 <- exp(- m.cum.br.rx.10)  
m.cum.br.rx.10 <- 1- s.cum.br.rx.10

m.br.rx.10 <- m.cum.br.rx.10
for (j in 1:2) {
  for (i in 2:rows) {
    m.br.rx.10[i,j] <- m.cum.br.rx.10[i,j] - m.cum.br.rx.10[i-1,j]
  }
}

# Cumulative all cause mortality conditional on surviving breast and all cause mortality
m.cum.all.rx.10 <- 1 - s.cum.oth.10*s.cum.br.rx.10
s.cum.all.rx.10 <- 100-100*m.cum.all.rx.10

# Annual all cause mortality
m.all.rx.10 <- m.cum.all.rx.10
for (j in 1:2) {
  for (i in 2:rows) {
    m.all.rx.10[i,j] <- m.cum.all.rx.10[i,j] - m.cum.all.rx.10[i-1,j]
  }
}

# Proportion of all cause mortality that is breast cancer
prop.br.rx.10      <- m.br.rx.10/(m.br.rx.10 + m.oth.10)
pred.m.br.rx.10    <- prop.br.rx.10*m.all.rx.10
pred.cum.br.rx.10  <- apply(pred.m.br.rx.10, 2, cumsum) 
pred.m.oth.rx.10   <- m.all.rx.10 - pred.m.br.rx.10
pred.cum.oth.rx.10 <- apply(pred.m.oth.rx.10, 2, cumsum)
pred.cum.all.rx.10 <- pred.cum.br.rx.10 + pred.cum.oth.rx.10

# rx benefits
benefits2.1.10 <- 100*(pred.cum.all.rx.10[,1] - pred.cum.all.rx.10)

rm("m.all.rx",
   "m.all.rx.10",
   "m.br.rx",
   "m.br.rx.10",
   "m.cum.all.rx",
   "m.cum.all.rx.10",
   "m.cum.br.rx" ,
   "m.cum.br.rx.10",
   "m.cum.oth", 
   "m.cum.oth.10",
   "m.cum.oth.rx",
   "rel.cum.rx",
   "m.oth",
   "m.oth.10",
   "m.oth.rx",
#   "m.rel.rx",
   "prop.br.rx",
   "prop.br.rx.10",
   "prop.rel.rx",
   "r.base.br",
   "r.base.oth", "r.br", "r.oth", "s.cum.all.rx", "s.cum.all.rx.10",
   "s.cum.br.rx", "s.cum.br.rx.10", "s.cum.oth", "s.cum.oth.10",
   "s.cum.oth.rx", "s.cum.oth10", "s.cum.rel.rx")
                                                                                                                                 


