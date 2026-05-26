# To source file location
setwd("MixTwice")
list.files()
library(devtools)
rm(list = c("mixtwice")) # remove any from current R session
document() # generate package documentation file
library(MixTwice)
packageVersion("MixTwice")
args(mixtwice)
devtools::check(cran = TRUE)
mixtwice_v3 <- MixTwice::mixtwice


devtools::document()
devtools::test()
devtools::run_examples()
devtools::check(cran = TRUE)



detach("package:MixTwice", unload = TRUE, character.only = TRUE)
unloadNamespace("MixTwice")


# load version 2.0
library(remotes)
lib_v2 <- path.expand("~/R/lib_v2")
dir.create(lib_v2, recursive = TRUE, showWarnings = FALSE)
remotes::install_version(
  "MixTwice",
  version = "2.0",
  lib = lib_v2,
  upgrade = "never",
  force = TRUE
)
library(MixTwice,lib.loc = lib_v2)
packageVersion("MixTwice")
mixtwice_v2 <- MixTwice::mixtwice






