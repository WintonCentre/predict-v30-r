#
# Run this script to setup rstudio to develop this package.
# Then use BUILD>check and BUILD>install to install locally as a package
#
library(devtools)
use_package("tidyr")
use_package("utils")
use_package("dplyr")
use_package("readr")
use_tibble()
use_import_from("tibble", "as_tibble")
use_import_from("dplyr", "mutate")
use_import_from("readr", "read_csv")
use_import_from("dplyr", 'case_when')
use_pipe(export=TRUE)
use_import_from("tidyr", "pivot_longer")
devtools::document()
