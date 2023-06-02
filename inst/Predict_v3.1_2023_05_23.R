library(tidyverse)

# Version with size  parameterised as  1 -  exp(-size/20)

coeffs <- read_csv("coefficients_v3.csv")

ag1_er0 <- unlist(coeffs[1,2])
ag2_er0 <- unlist(coeffs[2,2])
sz1_er0 <- unlist(coeffs[3,2])
nd1_er0 <- unlist(coeffs[4,2])
gr1_er0 <- unlist(coeffs[5,2])
sc1_er0 <- unlist(coeffs[6,2])
yr1_er0 <- unlist(coeffs[7,2])
ag1_er1 <- unlist(coeffs[8,2])
ag2_er1 <- unlist(coeffs[9,2])
sz1_er1 <- unlist(coeffs[10,2])
nd1_er1 <- unlist(coeffs[11,2])
gr1_er1 <- unlist(coeffs[12,2])
sc1_er1 <- unlist(coeffs[13,2])
yr1_er1 <- unlist(coeffs[14,2])
ag_ot_ea.beta.1 <- unlist(coeffs[15,2])
ag_ot_ea.beta.2 <- unlist(coeffs[16,2])
yr_ot_ea.beta <- unlist(coeffs[17,2])
h0_br_i <- unlist(coeffs[18,2])
h0_br_t1 <- unlist(coeffs[19,2])
h0_br_t2 <- unlist(coeffs[20,2])
h1_br_i <- unlist(coeffs[21,2])
h1_br_t1 <- unlist(coeffs[22,2])
h1_br_t2 <- unlist(coeffs[23,2])
h_ot_i <- unlist(coeffs[24,2])
h_ot_t1 <- unlist(coeffs[25,2])
h_ot_t2 <- unlist(coeffs[26,2])

## ----input, echo=FALSE---------------------------------------------------
year       <- 2017  # For web tool this is not an optional input - fixed at 2017
age.start  <- 50
screen     <- 0     # Clinically detected = 0, screen detected = 1
size       <- 20    # Tumour size mm
grade      <- 2     # Tumour grade
nodes      <- 2     # Number positive nodes
er         <- 1     # ER+ = 1, ER- = 0
her2       <- 1     # HER2+ = 1, HER2- = 0, missing = 9
ki67       <- 0     # KI67+ = 1, KI67- = 0, missing = 9
pr         <- 0     # PR+ = 1, PR- = 0, missing = 9

## --- treatment
radio      <- 1     # Radiotherapy Yes = 1, no = 0
horm       <- 1     # Hormone therapy Yes = 1, no = 0
generation <- 3     # Chemo generation 0, 2 or 3 only
traz       <- 1     # Trastuzumab therapy Yes = 1, no = 0
bis        <- 1     # Bisphosphonate therapy Yes = 1, no = 0
heart.gy   <- 4     # No grays radiotherapy to heart

## --- lifestyle
smoker     <- 1     # Never/ex = 0, current = 1

##----------------------------------------------------------------
##[WINTON FIX] Fix inputs for missing data
screen    <- ifelse(screen == 2, 0.204, screen)
grade     <- ifelse(grade == 9, 2.13, grade)

## ------------------------------------------------------------------------
time      <- c(1:15)
age       <- age.start - 24

##[WINTON FIX] - Input changed to include grade = 9
# Grade variable for ER neg now takes on 3 values 1, 2, 3

## ------------------------------------------------------------------------
age.mfp.1   <- ifelse(er==1, ((age/100)^-0.5), (age/100))
age.beta.1  <- ifelse(er==1, ag1_er1, ag1_er0)
age.mfp.2   <- ifelse(er==1, ((age/100)^2), (age/100) * log(age/100))
age.beta.2  <- ifelse(er==1, ag2_er1, ag2_er0)
size.mfp  <- ifelse(er==1, 1-exp(-size/20), log(size))
size.beta <- ifelse(er==1, sz1_er1, sz1_er0)
nodes.mfp   <- ifelse(er==1, log(nodes + 1), log(nodes + 1))
nodes.beta  <- ifelse(er==1, nd1_er1, nd1_er0)
grade.beta  <- ifelse(er==1, gr1_er1, gr1_er0)
screen.beta <- ifelse(er==1, sc1_er1, sc1_er0)
year.br.beta  <- ifelse(er==1, yr1_er1, yr1_er0)
## HER2 time varying effect for ER positive disease
her2.beta   <- case_when(her2==1 & er==1 ~ c(0.608, .532, .457, .382, .307, .231, .156, .081, .006, 0, 0, 0, 0, 0, 0),
                         her2==0 & er==1 ~ c(-.053, -.046, -.040, -.033, -.027, -.020, -.014, -.007, 0, 0, 0, 0, 0, 0, 0),
                         her2==1 & er==0 ~ 0.2316,
                         her2==0 & er==0 ~ -0.08589,
                         TRUE ~ 0)
ki67.beta   <- case_when(ki67==1 & er==1 ~ 0.14904,
                         ki67==0 & er==1 ~ -0.11333,
                         TRUE ~0 )
pr.beta    <- case_when(pr==1 & er==1 ~ -0.0619,
                        pr==0 & er==1 ~ 0.2624,
                        pr==1 & er==0 ~ -0.2231,
                        pr==0 & er==0 ~ 0.0296,
                        TRUE ~0)

## --- smoking adjustment ------------------------------------------------
smoker.prop  <- 0.1  # Proportion of population that are current smokers
smoker.rr    <- 2    # Relative risk non-br mortality in smokers
cvd.prop     <- 0.25 # Proportion of non-br mortality due to smoking related disease
                     # Proportion of I2*, I6*, C15 and C34 in SEARCH breast deaths
smoker.rr.acm <- cvd.prop*smoker.rr + 1 - cvd.prop

smoker.beta <- ifelse(smoker==0, log(1/(1 - smoker.prop + smoker.rr.acm*smoker.prop)),
                      log(smoker.rr.acm/(1 - smoker.prop + smoker.rr.acm*smoker.prop))) # Assume smoking RR of 2 and 10% population smokers

## ----baseline_adjust-----------------------------------------------------
c.other  <- 1.2  # RH non-breast mortality from Kerr et al (2022)
r.prop   <- 0.69 # Proportion of population receiving radiotherapy
r.breast <- 0.82 # Relative hazard breast specific mortality from Darby et al
r.other  <- 1.04 # Relative hazard other mortality per Gy heart dose from Taylor et al (2017)

r.base.br  <- log(1/((1-r.prop) + r.prop*r.breast))
r.base.oth <- log(1/((1-r.prop) + r.prop*(r.other^2))) # Assume 2 Gy average heart dose

## ------------------------------------------------------------------------
# Other mortality prognostic index (mi)
c.oth <- log(c.other)
r.oth <- log(r.other)*heart.gy
mi <-  ag_ot_ea.beta.1*((age/100)^3) + ag_ot_ea.beta.2*((age/100)^3*log(age/100))  + yr_ot_ea.beta*(year-2000) +
  r.base.oth + smoker.beta

# Breast cancer mortality prognostic index (pi)
pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
  nodes.beta*nodes.mfp + grade.beta*grade +
  screen.beta*screen + year.br.beta*(year-2000) +
  her2.beta + ki67.beta + pr.beta + r.base.br

c     <- ifelse(generation == 0, 0, ifelse(generation == 2, -.248, -.446))
h     <- ifelse(horm==1 & er==1, -0.3857, 0)
h10  <- if(h==0) {0
  } else {c(rep(h, 10), rep(-.26+h, 5)) } #including both ATLAS and aTTom trials
t     <- ifelse(her2==1 & traz==1, -.3567, 0)
b     <- ifelse(bis==1, -0.198, 0) # Only applicable to menopausal women.
r  <- ifelse(radio==1, log(r.breast), 0)

rx <- tibble(s = rep(0, 15),
             r = r,
             rh = r + h,
             rc = r + c,
             rt = r + t ,
             rb = r + b,
             rhc = r + h + c,
             rht = r + h + t,
             rhb = r + h + b,
             rct = r + c + t,
             rcb = r + c + b,
             rtb = r + t + b,
             rhct = r + h + c + t,
             rhcb = r + h + c + b,
             rhtb = r + h + t + b,
             rctb = r + c + t + b,
             rhctb = r + h + c + t + b,
             h = h,
             hc = h + c,
             ht = h + t,
             hb = h + b,
             hct = h + c + t,
             hcb= h + c + b,
             htb = h + t + b,
             hctb = h + c + t + b,
             c = c,
             ct =  c + t,
             cb = c + b,
             ctb = c + t + b,
             t = t,
             tb = t + b,
             rh10 = r + h10,
             rh10c = r + h10 + c,
             rh10t = r + h10 + t,
             rh10b = r + h10 + b,
             rh10ct = r + h10 + c + t,
             rh10cb = r + h10 + c + b,
             rh10tb = r + h10 + t + b,
             rh10ctb = r + h10 + c + t + b,
             h10 = h10,
             h10c = h10 + c,
             h10t = h10 + t,
             h10b = h10 + b,
             h10ct = h10 + c + t,
             h10cb= h10 + c + b,
             h10tb = h10 + t + b,
             h10ctb = h10 + c + t + b
             )

pi.rx <- rx + pi

cols <- ncol(rx)

mi.rx <- tibble(s = rep(mi, 15),
                r = mi + r.oth,
                rh = mi + r.oth,
                rc = mi + r.oth + c.oth,
                rt = mi + r.oth ,
                rb = mi + r.oth,
                rhc = mi + r.oth + c.oth,
                rht = mi + r.oth,
                rhb = mi + r.oth,
                rct = mi + r.oth + c.oth,
                rcb = mi + r.oth + c.oth,
                rtb = mi + r.oth,
                rhct = mi + r.oth + c.oth,
                rhcb = mi + r.oth + c.oth,
                rhtb = mi + r.oth,
                rctb = mi + r.oth + c.oth,
                rhctb = mi + r.oth + c.oth,
                h = mi,
                hc = mi + c.oth,
                ht = mi,
                hb = mi,
                hct = mi + c.oth,
                hcb= mi + c.oth,
                htb = mi,
                hctb = mi + c.oth,
                c = mi + c.oth,
                ct =  mi + c.oth,
                cb = mi + c.oth,
                ctb = mi + c.oth,
                t = mi,
                tb = mi,
                rh10 = mi + r.oth,
                rh10c = mi + r.oth + c.oth,
                rh10t = mi + r.oth,
                rh10b = mi + r.oth,
                rh10ct = mi + r.oth + c.oth,
                rh10cb = mi + r.oth + c.oth,
                rh10tb = mi + r.oth,
                rh10ctb = mi + r.oth + c.oth,
                h10 = mi,
                h10c = mi + c.oth,
                h10t = mi,
                h10b = mi,
                h10ct = mi + c.oth,
                h10cb= mi + c.oth,
                h10tb = mi,
                h10ctb = mi + c.oth
)

## ------------------------------------------------------------------------
# Generate cumulative baseline other mortality
base.m.cum.oth <- exp(h_ot_i+ h_ot_t1*log(time/10) + h_ot_t2*(time/10))

# Generate annual baseline other mortality
base.m.oth <- base.m.cum.oth
for (i in 2:15) {
  base.m.oth[i] <- base.m.cum.oth[i] - base.m.cum.oth[i-1] }

# Generate the annual other cancer specific mortality rate
m.oth <- base.m.oth*exp(mi.rx)

# Calculate the cumulative other cancer mortality rate
m.cum.oth <- apply(m.oth, 2, cumsum)

# Calculate the cumulative breast cancer survival
s.cum.oth <- exp(- m.cum.oth)
m.cum.oth <- 1 - s.cum.oth

m.oth <- m.cum.oth
for (j in 1:cols) {
  for (i in 2:15) {
    m.oth[i,j] <- m.cum.oth[i,j] - m.cum.oth[i-1,j]
  }
}

## ------------------------------------------------------------------------
# Generate cumulative baseline breast mortality
if (er==1) {
  base.m.cum.br <- exp(h1_br_i + h1_br_t1*((time/10)^-0.5) + h1_br_t2*((time/10)^-0.5*log(time/10)))
} else {exp(h0_br_i + h0_br_t1*((time/10)^-1) + h0_br_t2*((time/10)^-1*log(time/10)))
}

# Generate annual baseline breast mortality
base.m.br <- base.m.cum.br
for (i in 2:15) {
  base.m.br[i] <- base.m.cum.br[i] - base.m.cum.br[i-1] }

# Generate the annual breast cancer specific mortality rate
m.br <- base.m.br*exp(pi.rx)

# Calculate the cumulative breast cancer mortality rate
m.cum.br <- apply(m.br, 2, cumsum)

# Calculate the cumulative breast cancer survival
s.cum.br <- exp(- m.cum.br)
m.cum.br <- 1- s.cum.br

m.br <- m.cum.br
for (j in 1:cols) {
  for (i in 2:15) {
    m.br[i,j] <- m.cum.br[i,j] - m.cum.br[i-1,j]
  }
}

# Cumulative all cause mortality conditional on surviving breast and all cause mortality
m.cum.all <- 1 - s.cum.oth*s.cum.br
s.cum.all <- 100 - 100*m.cum.all

# Annual all cause mortality
m.all <- m.cum.all
for (j in 1:cols) {
  for (i in 2:15) {
    m.all[i,j] <- m.cum.all[i,j] - m.cum.all[i-1,j]
  }
}

## ------------------------------------------------------------------------
# Proportion of all cause mortality that is breast cancer
prop.br      <- m.br/(m.br + m.oth)
pred.m.br    <- prop.br*m.all
pred.cum.br  <- apply(pred.m.br, 2, cumsum)
pred.m.oth   <- m.all - pred.m.br
pred.cum.oth <- apply(pred.m.oth, 2, cumsum)
pred.cum.all <- pred.cum.br + pred.cum.oth

surv <- as_tibble(100*(1-pred.cum.all))

benefits <- benefits <- tibble(r_benefit = round(100*(pred.cum.all[,1] - pred.cum.all[,2]), 1),
                               h_benefit = round(100*(pred.cum.all[,2] - pred.cum.all[,3]), 1),
                               c_benefit = round(100*(pred.cum.all[,3] - pred.cum.all[,7]), 1),
                               t_benefit = round(100*(pred.cum.all[,7] - pred.cum.all[,13]), 1),
                               b_benefit = round(100*(pred.cum.all[,13] - pred.cum.all[,17]), 1),
                               total_benefit = round(100*(pred.cum.all[,1] - pred.cum.all[,17]), 1))

## -----------------------------------------------------------------------
# Calculating the benefit of continuing endocrine therapy
# assuming sruvival to 5 years
start <- 6

# Generate the annual other cancer specific mortality rate
m.oth <- base.m.oth[start:15]*exp(mi.rx[start:15,])

# Calculate the cumulative other cancer mortality rate
m.cum.oth <- apply(m.oth, 2, cumsum)

# Calculate the cumulative oth cancer survival
s.cum.oth <- exp(- m.cum.oth)
m.cum.oth <- 1- s.cum.oth

m.oth <- m.cum.oth
for (j in 1:cols) {
  for (i in 2:10) {
    m.oth[i,j] <- m.cum.oth[i,j] - m.cum.oth[i-1,j]
  }
}

# Generate the annual breast cancer specific mortality rate
m.br <- base.m.br[start:15]*exp(pi.rx[start:15,])

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

surv_conditional <- as_tibble(100*(1-pred.cum.all)) #[c(10), c(1:3)])

#rm(list=setdiff(ls(), c("benefits", "surv", "surv_conditional")))

