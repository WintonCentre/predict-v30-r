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
#' @export benefits3110
#'
#' @examples benefits3110()

benefits3110 <- function(age.start.in  = 65,
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
  locals3110 <- list2env(benefits31(
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

  options(pillar.sigfig = 4)

  # Version with size  parameterised as  1 -  exp(-size/20)

  #coeffs probably doesn't need to be retuened as is in locals
  csv_path <- system.file("extdata", "coefficients_v3.csv", package = "PREDICT3")
  coeffs <- read_csv(csv_path, show_col_types = FALSE)

  # Move into locals3110 environment if you need to see these...
  # NB. We could inherit the coeffs from benefits31() with a bit of refactoring
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

  with(locals3110, {
    ## -----------------------------------------------------------------------
    # Calculating the benefit of continuing endocrine therapy
    # assuming sruvival to 5 years
    start <- 6

    # New version of m.oth, based on benefits31 directly
    # m.oth <- m.oth[start:15,]

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

    # For debug purposes
    s.cum.oth.ten <- m.oth
      for (j in 1:cols) {
        s.cum.oth.ten[,j] <- 1-cumsum(m.oth[,j])
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
    return(as.list.environment(locals3110))
  })
}
