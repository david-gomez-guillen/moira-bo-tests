# The valley problem, whose code is in its folder. Its files are read relative
# to that folder, both while it is sourced and later, when THALASSA asks for the
# overview or the code sample, so the working directory stays there.
source('7-valley/thalassa_interface.R', chdir = TRUE)
setwd('7-valley')
