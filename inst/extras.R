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
