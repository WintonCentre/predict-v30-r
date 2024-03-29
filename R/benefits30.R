## ----setup, include=FALSE------------------------------------------------
#knitr::opts_chunk$set(echo = TRUE)

#' Title Calculate Predict3.0 treatment benefits given patient details
#'
#' @param age.start.in Patient age in years
#' @param screen.in Clinically detected = 0, screen detected = 1
#' @param size.in Tumour size mm
#' @param grade.in Tumour grade
#' @param nodes.in Number positive nodes
#' @param er.in ER+ = 1, ER- = 0
#' @param her2.in HER2+ = 1, HER2- = 0, missing = 9
#' @param ki67.in KI67+ = 1, KI67- = 0, missing = 9
#' @param pr.in ptogeserpne satus PR+ = 1 PR- = 0
#' @param generation.in Chemo generation 0, 2 or 3 only
#' @param horm.in Hormone therapy Yes = 1, no = 0
#' @param traz.in Trastuzumab therapy Yes = 1, no = 0
#' @param bis.in Bisphosphonate therapy Yes = 1, no = 0
#' @param radio.in Radiotherapy Yes = 1, no = 0
#' @param heart.gy.in Number of grays radiation to the heart
#' @param smoker.in Never/ex = 0, current = 1
#'
#' @return table of treatment benefits
#' @export benefits30
#'
#' @examples benefits30()


benefits30 <- function(
    age.start.in  = 65,
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
    smoker.in     = 1     # Never/ex = 0, current = 1
) {
  #set tibble print precision
  options(pillar.sigfig = 4)


  # Create an environment for local vars
  locals30 <- new.env()
  with(locals30, {

    age.start  <- age.start.in;
    screen     <- screen.in;     # Clinically detected <- 0, screen detected = 1
    size       <- size.in;    # Tumour size mm
    grade      <- grade.in;     # Tumour grade
    nodes      <- nodes.in;     # Number positive nodes
    er         <- er.in;     # ER+ = 1, ER- = 0
    her2       <- her2.in;     # HER2+ = 1, HER2- = 0, missing = 9
    ki67       <- ki67.in;     # KI67+ = 1, KI67- = 0, missing = 9
    pr         <- pr.in;     # PR+ = 1, PR- = 0, missing = 9

    ## --- treatment
    generation <- generation.in;     # Chemo generation 0, 2 or 3 only
    horm       <- horm.in;     # Hormone therapy Yes = 1, no = 0
    traz       <- traz.in;     # Trastuzumab therapy Yes = 1, no = 0
    bis        <- bis.in;     # Bisphosphonate therapy Yes = 1, no = 0
    radio      <- radio.in;     # Radiotherapy Yes = 1, no = 0
    heart.gy   <- heart.gy.in;     # No grays radiotherapy to heart

    smoker <- smoker.in;

    ##----------------------------------------------------------------
    ##[WINTON FIX] Fix inputs for missing data
    screen    <- ifelse(screen == 2, 0.204, screen)
    grade     <- ifelse(grade == 9, 2.13, grade)

    ## ------------------------------------------------------------------------
    time      <- c(1:15)
    age       <- age.start - 1 + time
    ##[WINTON FIX] - Input changed to include grade = 9
    grade.val <- ifelse(er==1, grade, ifelse(grade>=2, 1, 0)) # Grade variable for ER neg

    ## ------------------------------------------------------------------------
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
                          ifelse(ki67==0 & er==1, -0.11333,0 ))
    pr.beta    <- ifelse(pr==1 & er==1, -0.0619,
                         ifelse(pr==0 & er==1, 0.2624,
                                ifelse(pr==1 & er==0, -0.2231,
                                       ifelse(pr==0 & er==0, 0.0296, 0))))
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
    c.oth <- ifelse(generation==0, 0, log(c.other))
    r.oth <- ifelse(radio==0, 0, log(r.other)*heart.gy)
    zagebase <- 0.0698252*((age.start/10)^2-34.23391957)
    mi <- zagebase + r.base.oth + smoker.beta + c.oth + r.oth
    #mi <- 0.0698252*((age.start/10)^2-34.23391957) + r.base.oth + smoker.beta + c.oth + r.oth

    # Breast cancer mortality prognostic index (pi)
    pi <- age.beta.1*age.mfp.1 + age.beta.2*age.mfp.2 + size.beta*size.mfp +
      nodes.beta*nodes.mfp + grade.beta*grade.val + screen.beta*screen +
      her2.beta + ki67.beta + pr.beta + r.base.br

    c     <- ifelse(generation == 0, 0, ifelse(generation == 2, -.248, -.446))
    h     <- ifelse(horm==1 & er==1, -0.3857, 0)

    # h10 betas are the same as h betas for 10 years, and then gain another -.26 benefit in years 11 to 15.
    # Note that this means the benefit of taking extended hormone from year 6 only shows up 5 years later.
    h10  <- ifelse(h==0, 0, c(rep(h, 10), rep(-.26+h, 5))) #including both ATLAS and aTTom trials
    t     <- ifelse(her2==1 & traz==1, -.3567, 0)
    b     <- ifelse(bis==1, -0.198, 0) # Only applicable to menopausal women.
    r.br  <- ifelse(radio==1, log(r.breast), 0)

    rxbase <- tibble(surg = rep(0, 15),

                 h = h,
                 hr = h + r.br,
                 hrc = h + r.br + c,
                 hrct = h + r.br + c + t,
                 hrctb = h + r.br + c + t + b,


                 h10 = h10,
                 h10r = h10 + r.br,
                 h10rc = h10 + r.br + c,
                 h10rct = h10 + r.br + c + t,
                 h10rctb = h10 + r.br + c + t + b
                 )

    rx <- rxbase + pi

    # export treatment combination names for clojurescript test purposes. Otherwise this is unused.
    treatments <- colnames(rx)

    cols <- ncol(rxbase)

    ## ------------------------------------------------------------------------
    # Generate cumulative baseline other mortality
    base.m.cum.oth <- exp(-6.052919 + (1.079863*log(time)) + (.3255321*time^.5))

    # Generate cumulative survival non-breast mortality
    # Incorporates the increased mortality associated with chemo and radiotherapy
    # *** WINTON Fix: c.oth and r.oth have already been included, so don't include again here. ***
    #s.cum.oth <- exp(-exp(mi + c.oth + r.oth)*base.m.cum.oth)
    s.cum.oth <- exp(-exp(mi)*base.m.cum.oth)

    # Convert cumulative mortality rate into cumulative risk
    m.cum.oth <- 1- s.cum.oth

    # Annual other mortality rate
    m.oth <- m.cum.oth
    for (i in 2:15) {
        m.oth[i] <- m.cum.oth[i] - m.cum.oth[i-1]
    }

    #
    # NB This m.oth is carried forward into the benefits3010 calculation
    #

    ## ------------------------------------------------------------------------
    # Generate cumulative baseline breast mortality
    if (er==1) {
      base.m.cum.br <- exp(0.7424402 - 7.527762/time^.5 - 1.812513*log(time)/time^.5)
    } else { base.m.cum.br <- exp(-1.156036 + 0.4707332/time^2 - 3.51355/time)
    }

    # Generate annual baseline breast mortality
    base.m.br <- base.m.cum.br
    for (i in 2:15) {
      base.m.br[i] <- base.m.cum.br[i] - base.m.cum.br[i-1] }

    #
    # This base.m.br is carried forward into the benefits3010 calculation
    #

    # Generate the annual breast cancer specific mortality rate.
    # benefits3010 repeats this calculation, but with base.m.br and rx resampled over years 6:15
    # See benefits3010.R from line 79.
    m.br.1 <- base.m.br*exp(rx)

    # Calculate the cumulative breast cancer mortality rate
    m.cum.br.1 <- apply(m.br.1, 2, cumsum)

    # Calculate the cumulative breast cancer survival
    s.cum.br <- exp(- m.cum.br.1)
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

    # surv <- 100*(1-pred.cum.all)[c(5,10,15), ]

    ## ------------------------------------------------------------------------
    # rx benefits
    # version implemented on web has benefit as difference in breast specific mortality
    benefits30 <- 100*(pred.cum.all[,1] - pred.cum.all)
    benefits30[,1] <- 100*(1 - pred.cum.all[,1]) # patch in baseline in surgery column


    # Original Paul Pharoah code commented out below. Not used in tests
    # benefits <- as_tibble(signif(100*(pred.cum.all[,1] - pred.cum.all), 3)[c(5,10,15),]) %>%
    #     mutate(year = c(5, 10, 15)) %>%
    #     pivot_longer(cols=1:3, names_to = "rx", values_to = "benefit")

    return(as.list.environment(locals30))
  })
}

# Call benefits30 with parameters corresponding to seed=100 in the cljs parameter generator
# See predict.models.v3_opencpu_test/generate-r-call.
s100 <- function() {
  benefits30(  age.start.in = 33,
                     bis.in =  1,
                     er.in =  1,
                     generation.in =  0,
                     grade.in = 1,
                     heart.gy.in = 7,
                     her2.in = 1,
                     horm.in = 0,
                     ki67.in = 1,
                     nodes.in =  3,
                     pr.in =  1,
                     radio.in = 1,
                     screen.in =  1,
                     size.in =  21,
                     smoker.in = 1,
                     traz.in = 1)
}
s200 <- function() {
  benefits30(bis.in=1,age.start.in=40,size.in=3,screen.in=1,heart.gy.in=12,nodes.in=0.5,generation.in=0,pr.in=1,smoker.in=1,horm.in=0,grade.in=1,er.in=1,ki67.in=1,traz.in=1,her2.in=1,radio.in=1)
}
s251 <- function() {
  benefits30(bis.in=0,age.start.in=55,size.in=40,screen.in=0,heart.gy.in=25,nodes.in=4,generation.in=2,pr.in=0,smoker.in=1,horm.in=1,grade.in=2,er.in=1,ki67.in=1,traz.in=1,her2.in=1,radio.in=1)
}
