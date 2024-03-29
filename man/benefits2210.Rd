% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benefits2210.R
\name{benefits2210}
\alias{benefits2210}
\title{v2.2 predict survival benefits years 1 to 15}
\usage{
benefits2210(
  age.start = 57,
  screen = 1,
  size = 20,
  grade = 2,
  nodes = 10,
  er = 1,
  her2 = 1,
  ki67 = 1,
  generation = 3,
  horm = 1,
  traz = 0,
  bis = 0,
  radio = 0,
  delay = 5
)
}
\arguments{
\item{age.start}{Age at time of surgery}

\item{screen}{Clinically detected = 0, screen detected = 1}

\item{size}{Tumour size mm}

\item{grade}{Tumour grade}

\item{nodes}{Number positive nodes}

\item{er}{ER+ = 1, ER- = 0}

\item{her2}{HER2+ = 1, HER2- = 0, missing = 9}

\item{ki67}{KI67+ = 1, KI67- = 0, missing = 9}

\item{generation}{Chemo generation 0, 2 or 3 only}

\item{horm}{Hormone therapy Yes = 1, no = 0}

\item{traz}{Trastuzumab therapy Yes = 1, no = 0}

\item{bis}{Bisphosphonate therapy Yes = 1, no = 0}

\item{radio}{Radiotherapy Yes = 1, no = 0}

\item{delay}{0 or 5 years since surgery for h10 benefit}
}
\description{
h10 over h5 benefits with delay 0 or 5.
}
\examples{
benefits2210()
}
