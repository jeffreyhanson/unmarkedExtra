# compile c++ attributes
library(Rcpp)
setwd(file.path(Sys.getenv('HOME'), 'github', 'unmarkedExtra'))
compileAttributes()

# document code
setwd(file.path(Sys.getenv('HOME'), 'github', 'unmarkedExtra'))
library(devtools)
library(roxygen2)
document()

# find obvious errors
setwd(file.path(Sys.getenv('HOME'), 'github', 'unmarkedExtra'))
library(devtools)
library(roxygen2)
load_all() 

# formal package tests
setwd(file.path(Sys.getenv('HOME'), 'github', 'unmarkedExtra'))
library(devtools)
library(roxygen2)
test()

# local install
library(devtools)
setwd(file.path(Sys.getenv('HOME'), 'github', 'unmarkedExtra'))

# install from github
library(devtools)
install_github('paleo13/marxan')

# make vignettes
library(devtools)
setwd(file.path(Sys.getenv('HOME'), 'github', 'unmarkedExtra'))
build_vignettes()

