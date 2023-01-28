## ----setup, include=FALSE------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

#' Title Calculate benefits of continuing endocrine therapy assuming survival to 5 years
#'
#' @param age.start.in Patient age in years
#' @param screen.in Clinically detected = 0, screen detected = 1
#' @param size.in Tumour size mm
#' @param grade.in Tumour grade
#' @param nodes.in Number positive nodes
#' @param er.in ER+ = 1, ER- = 0
#' @param her2.in HER2+ = 1, HER2- = 0, missing = 9
#' @param ki67.in KI67+ = 1, KI67- = 0, missing = 9
#' @param pr.in ptogesterone satus PR+ = 1 PR- = 0
#' @param generation.in Chemo generation 0, 2 or 3 only
#' @param horm.in Hormone therapy Yes = 1, no = 0
#' @param traz.in Trastuzumab therapy Yes = 1, no = 0
#' @param bis.in Bisphosphonate therapy Yes = 1, no = 0
#' @param radio.in Radiotherapy Yes = 1, no = 0
#' @param heart.gy.in Number of grays radiation to the heart
#' @param smoker.in Never/ex = 0, current = 1
#'
#' @return list of environment variables
#' @export benefits3010
#'
#' @examples benefits3010()

benefits3010 <- function(age.start.in  = 65,
                         screen.in     = 0,     # Clinically detected = 0, screen detected = 1
                         size.in       = 25,    # Tumour size mm
                         grade.in      = 2,     # Tumour grade
                         nodes.in      = 2,     # Number positive nodes
                         er.in         = 1,     # ER+ = 1, ER- = 0
                         her2.in       = 0,     # HER2+ = 1, HER2- = 0, missing = 9
                         ki67.in       = 1,     # KI67+ = 1, KI67- = 0, missing = 9
                         pr.in         = 1,     # PR+ = 1, PR- = 0, missing = 9

                         ## --- treatment
                         generation.in = 2,     # Chemo generation 0, 2 or 3 only
                         horm.in       = 1,     # Hormone therapy Yes = 1, no = 0
                         traz.in       = 0,     # Trastuzumab therapy Yes = 1, no = 0
                         bis.in        = 1,     # Bisphosphonate therapy Yes = 1, no = 0
                         radio.in      = 1,     # Radiotherapy Yes = 1, no = 0
                         heart.gy.in   = 1,     # No grays radiotherapy to heart

                         ## --- lifestyle
                         smoker.in     = 1)     # Never/ex = 0, current = 1
{
  locals3010 <- list2env(benefits30(
    age.start.in,
    screen.in,
    size.in,
    grade.in,
    nodes.in,           # Number positive nodes
    er.in,
    her2.in,
    ki67.in,
    pr.in,

    ## --- treatment
    generation.in,
    horm.in,
    traz.in,
    bis.in,
    radio.in,
    heart.gy.in,

    ## --- lifestyle
    smoker.in))
  with(locals3010, {
    ## -----------------------------------------------------------------------
    # Calculating the benefit of continuing endocrine therapy
    # assuming survival to 5 years
    start <- 6
    m.oth <- m.oth[start:15]
    s.cum.oth.ten <- 1-cumsum(m.oth)

    # Generate the annbreast cancer specific mortality rate
    m.br0 <- base.m.br[start:15]*exp(rx[start:15,])

    # Calculate the cumulative breast cancer mortality rate
    m.cum.br <- apply(m.br0, 2, cumsum)

    # Calculate the cumulative breast cancer survival
    s.cum.br <- exp(- m.cum.br)
    m.cum.br <- 1- s.cum.br

    m.br <- m.cum.br
    for (j in 1:cols) {
      for (i in 2:10) {       # In this R version, years 6:15 correspond to  i in 2:10 (in R inclusive notation)
                              # But shouldn't this be i in 2:11 as 2:10 only covers years 6:14 ??
        m.br[i,j] <- m.cum.br[i,j] - m.cum.br[i-1,j]
      }
    }

    # Cumulative all cause mortality conditional on surviving breast and all cause mortality
    m.cum.all <- 1 - s.cum.oth.ten*s.cum.br
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

    surv_conditional <- 100*(1-pred.cum.all)[c(1:10), c(1:3)]
    return(as.list.environment(locals3010))
  })
}
