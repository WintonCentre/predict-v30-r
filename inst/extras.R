prbeta <- function(erx,prx) {
  ifelse(prx==1 & erx==1, -0.0619,
         ifelse(prx==0 & erx==1, 0.2624,
                ifelse(prx==1 & erx==0, -0.2231,
                       ifelse(prx==0 & erx==0, 0.0296, 0))))
}
pr.df <- tibble(
  x = 0:3,
  er = x%%2,
  pr = as.integer(x/2),
  beta = prbeta(erx=er,prx=pr)
)

# prbeta <- function(erx,prx) {
#   ifelse(prx==1 & erx==1, -0.0619,
#          ifelse(prx==0 & erx==1, 0.2624,
#                 ifelse(prx==1 & erx==0, -0.2231,
#                        ifelse(prx==0 & erx==0, 0.0296, 0))))

# library(jsonify)
test_data <- function() {
  params <- new.env(parent = environment());
  with(params, {
  age.start  <- 25;
  screen     <- 0;     # Clinically detected = 0, screen detected = 1
  size       <- 25;    # Tumour size mm
  grade      <- 2;     # Tumour grade
  nodes      <- 2;     # Number positive nodes
  er         <- 1;     # ER+ = 1, ER- = 0
  her2       <- 0;     # HER2+ = 1, HER2- = 0, missing = 9
  ki67       <- 1;     # KI67+ = 1, KI67- = 0, missing = 9
  pr         <- 1;     # PR+ = 1, PR- = 0, missing = 9

  ## --- treatment
  generation <- 2;     # Chemo generation 0, 2 or 3 only
  horm       <- 1;     # Hormone therapy Yes = 1, no = 0
  traz       <- 0;     # Trastuzumab therapy Yes = 1, no = 0
  bis        <- 1;     # Bisphosphonate therapy Yes = 1, no = 0
  radio      <- 1;     # Radiotherapy Yes = 1, no = 0
  heart.gy   <- 1;     # No grays radiotherapy to heart

  ## --- lifestyle
  smoker     <- 1;     # Never/ex = 0, current = 1
  })
  return(params)
}

# dump <- function(params, file = "example1.json") {
#   x <- to_json(params, unbox=TRUE);
#   write(x, file)
# }
#
# dump(test_data())
