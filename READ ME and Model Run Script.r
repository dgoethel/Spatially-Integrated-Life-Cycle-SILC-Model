rm(list=(ls()))

# following are packages that need to be installed to run the various r scripts (not all of them are used, but I have forgotten which ones are needed at this point)
# you also need to have some version of ADMB installed (http://www.admb-project.org/downloads/admb-12.2/)
# (I like to use AD Studio which comes packaged with the EMACS text editor and syntax highlighting)

load_libraries<-function() {
  library(PBSmodelling)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(gplots)
  library(colorspace)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(matrixStats) 
  library(gridExtra)
  library(grid)
  library(gtools)
  library(TeachingDemos)
  library(snowfall)
  library(parallel)
  library(snow)
  library(foreach)
  library(doSNOW)
  library(spatstat)
  library(alphahull)
  library(beanplot)
  library(png)
  library(PBSadmb)
  library(gtable)
  library(corrplot)
  library(ggforce)
}
load_libraries()

rm(list=(ls()))

###next line sets the wd where you should have the r scripts and the OM/EM folders saved (likely just where this script is saved too)
wd<<-'D:\\NOAA FILES\\Research\\COCOCHA\\Hake Model\\Base Model Development\\_Scenario 1_No Move EM'
setwd(wd)
source('SIM_TIM_editing.r') 

wd<<-'D:\\NOAA FILES\\Research\\COCOCHA\\Hake Model\\Base Model Development\\_Scenario 2_Cnst Move EM'
setwd(wd)
source('SIM_TIM_editing.r') 

wd<<-'D:\\NOAA FILES\\Research\\COCOCHA\\Hake Model\\Base Model Development\\_Scenario 3_Time+Age Move EM'
setwd(wd)
source('SIM_TIM_editing.r') 

wd<<-'D:\\NOAA FILES\\Research\\COCOCHA\\Hake Model\\Base Model Development\\_Scenario 4_Onto Move OM'
setwd(wd)
source('SIM_TIM_editing.r') 

wd<<-'D:\\NOAA FILES\\Research\\COCOCHA\\Hake Model\\Base Model Development\\_Scenario 5_DD Move OM'
setwd(wd)
source('SIM_TIM_editing.r') 

wd<<-'D:\\NOAA FILES\\Research\\COCOCHA\\Hake Model\\Base Model Development\\_Scenario 6_Age-1 Move Only OM'
setwd(wd)
source('SIM_TIM_editing.r') 

#this runs the simulation script, which will perform one iteration of the model, output the plots for the single iteration, then run the loop of simulation runs
#you can adjust 'SIM_TIM_editing.r directly (line 22) to determine how many simulation runs you want to perform (for our papers we use 500)
# everything should be automated so that you don't have to adjust code in any of the other R scripts and all the results will be put into a auto created folder 
# called 'Diagnostics'...this includes folders that have the outputs from both the converged and non-converged runs, as well as, a variety of diagnostic plots
# the pdfs called 'Model_Diagnostics', 'Tag_Residuals', 'Correlation_Matrix', and 'Standard_Error_Table' are the full model outputs from a single run of the model, which
# give an idea of residual fits, likelihood components, uncertainty in parameter estimates, etc... I think these are relatively self-explanatory (Model_Diagnostics is the main file to explore)
# The rest of the created files (mostly csv and the 'Simulation_Summary' pdf) provide the outputs summarizing the simulation results, mainly in terms of % bias (or relative error)
# The pdf gives a fairly complete run down of all the estimated parameters across the time series both in terms of % bias and the distribution of estimates compared to the true value
# For each manuscript we have developed a variety of more sophisticated summary graphics, but they are mostly just higher quality versions of what this code outputs
# Generally I feel like these results give you a good idea of how the model is performing and how your simulation parameters impact relative dynamics
# The .rep file (and .par file) are saved for each run so almost anything that is not graphed here can be grabbed from those and figures created if needed
# if there is a value you want that is not in the .rep file, then you can add a report function to the end of the TIM_EM.tpl file (in the estimation_model folder)
# to output the value and record it in the report file

# to change inputs for the SIM model just open the 'operating_model' folder and then open the TIM_OM.dat...this is the main input file
# I have tried to make this file as detailed as possible and explaining each of the options for all the different model switches
# it is probably not perfect, but I think you can get an idea of what is needed and/or what types of options can be included
# the main tricky part for the .dat file is ensuring the input data matches the indices/subscripts (e.g., if you change the nages or nyears, then each of the inputs will need to change to match these new indices)
# it is probably also worth noting that many of the options are not used (i.e., the switches are set to 0), even though they may have non-zero inputs
# For instance you might turn tagging data simulation off, but you still need to have inputs for all the tag data inputs
# since this is not a disseminated program I don't have the best error checking built in, sorry if it is a bit confusing
#  the Goethel and Berger (2017) and Goethel et al. (2019) have some decent descriptions of the framework and dynamics (the latter especially in the SM for the estimation model dynamics)
# further info can be found at our GitHub pages (https://github.com/dgoethel/Spatially-Integrated-Life-Cycle-SILC-Model) (most recent paper https://github.com/dgoethel/Spatial-Assessment-Simulator), tagging paper (https://github.com/dgoethel/tag-integrated-model), and general SPASAM site (https://github.com/KateBoz/SPASAM) 
#


