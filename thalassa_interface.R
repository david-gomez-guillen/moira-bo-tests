# The noisy problem, whose code is in its folder. Its files are read relative
# to that folder, both while it is sourced and later, when THALASSA asks for the
# overview or the code sample, so the working directory stays there.
source('4-noisy/thalassa_interface.R', chdir = TRUE)
setwd('4-noisy')
