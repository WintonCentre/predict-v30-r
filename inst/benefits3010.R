## ----setup, include=FALSE------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

#' Title Calculate benefits of continuing endocrine therapy assuming survival to 5 years
#'
#' @param age.start Patient age in years
#' @param screen Clinically detected = 0, screen detected = 1
#' @param size Tumour size mm
#' @param grade Tumour grade
#' @param nodes Number positive nodes
#' @param er ER+ = 1, ER- = 0
#' @param her2 HER2+ = 1, HER2- = 0, missing = 9
#' @param ki67 KI67+ = 1, KI67- = 0, missing = 9
#' @param pr ptogeserpne satus PR+ = 1 PR- = 0
#' @param generation Chemo generation 0, 2 or 3 only
#' @param horm Hormone therapy Yes = 1, no = 0
#' @param traz Trastuzumab therapy Yes = 1, no = 0
#' @param bis Bisphosphonate therapy Yes = 1, no = 0
#' @param radio Radiotherapy Yes = 1, no = 0
#' @param heart.gy Number of grays radiation to the heart
#' @param smoker Never/ex = 0, current = 1
#'
#' @return table of treatment benefits
#' @export benefits3010
#'
#' @examples benefits3010()

benefits3010 <- function() {
  locals <- benefits30(  age.start  = 65,
                         screen     = 0,     # Clinically detected = 0, screen detected = 1
                         size       = 25,    # Tumour size mm
                         grade      = 2,     # Tumour grade
                         nodes      = 2,     # Number positive nodes
                         er         = 1,     # ER+ = 1, ER- = 0
                         her2       = 0,     # HER2+ = 1, HER2- = 0, missing = 9
                         ki67       = 1,     # KI67+ = 1, KI67- = 0, missing = 9
                         pr         = 1,     # PR+ = 1, PR- = 0, missing = 9

                         ## --- treatment
                         generation = 2,     # Chemo generation 0, 2 or 3 only
                         horm       = 1,     # Hormone therapy Yes = 1, no = 0
                         traz       = 0,     # Trastuzumab therapy Yes = 1, no = 0
                         bis        = 1,     # Bisphosphonate therapy Yes = 1, no = 0
                         radio      = 1,     # Radiotherapy Yes = 1, no = 0
                         heart.gy   = 1,     # No grays radiotherapy to heart

                         ## --- lifestyle
                         smoker     = 1     # Never/ex = 0, current = 1
  );
  return(locals)
}

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

