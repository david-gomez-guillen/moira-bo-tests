# The bimodal problem, whose code is in its folder. Its files are read relative
# to that folder, both while it is sourced and later, when THALASSA asks for the
# overview or the code sample, so the working directory stays there.
source('1-bimodal/thalassa_interface.R', chdir = TRUE)
setwd('1-bimodal')
