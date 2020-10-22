#THIS IS A SCRIPT FOR OPTIMIZED SPIKE INFERENCE FOR A COMPLETE DATASET

#############################################################################################################

#USER PARAMETERS & SETTINGS

#IMPORT DATASET INTO R & GENERATE MATRIX
DATA_CSV = read.csv("C:\\Users\\YUSTE\\Documents\\Github\\caiman_sorter\\Good_Components.csv", header=FALSE)
DATA_SET = as.matrix(DATA_CSV)

#CLEANUP
rm(DATA_CSV)

#SET FRAMERATE & INDICATOR TYPES
#0.7 = f, 1.2 = m, 2 = s
framerate <- 7.725
indicator_type <- 0.7

#SET LAMBDA RANGE
LL = seq(from = 0, to = 0.1, by = 0.0005)
LL <- LL[-1]

##############################################################################################################

#IMPORT LIBRARIES
library(FastLZeroSpikeInference)
library(fractal)
library(STAR)

#FIND NUMBER OF NEURONS
num_neur = nrow(DATA_SET)

#SET NUMBER OF DATAPOINTS ON OPTIMIZATION PLOT
#ROWS = NEURONS, COL = LAMBDA
SPIKES <- matrix(0,nrow=num_neur, ncol=length(LL))
FALSE_SPIKES <- matrix(0, nrow=num_neur, ncol=length(LL))


##############################################################################################################

#MAKE FUNCTIONS
approx_gamma_decay <- function(framerate,indicator_type){
  delta = 1/framerate
  phi = indicator_type
  gamma_decay = 1 - (delta/phi)
  return(gamma_decay)
}

correction <- function(DATA_SET, neuron, order){
  c_dfof <- medianFilter(DATA_SET[neuron,],order)
  return(c_dfof)
}

intrinsic_noise_threshold <- function(c_dfof){
  in_noise_thr <- 2*(stdev(c_dfof))
  return(in_noise_thr)
}

false_spike_finder <- function(in_noise_thr, c_dfof, fit){
  false_spikes_vector <- c_dfof[fit[[1]]]
  adj_fsv <- (false_spikes_vector < in_noise_thr)
  false_spikes <- sum(adj_fsv)
  return(false_spikes)
}

single_loop <- function(DATA_SET, gamma_decay, L_lambda, EST_CAL, neuron, order){
  S = numeric(length=2)
  c_dfof <- correction(DATA_SET, neuron, order)
  in_noise_thr <- intrinsic_noise_threshold(c_dfof)
  fit <- estimate_spikes(dat = c_dfof, gam = gamma_decay, lambda = L_lambda, estimate_calcium = EST_CAL)
  spikes = length(fit[[1]])
  false_spikes = false_spike_finder(in_noise_thr, c_dfof, fit)
  S[1] = spikes
  S[2] = false_spikes
  return(S)
}

final_fit <- function(DATA_SET, gamma_decay, L_lambda, EST_CAL, neuron, order){
  c_dfof <- correction(DATA_SET, neuron, order)
  in_noise_thr <- intrinsic_noise_threshold(c_dfof)
  fit <- estimate_spikes(dat = c_dfof, gam = gamma_decay, lambda = L_lambda, estimate_calcium = EST_CAL)
  fit[["neuron"]]=c(fit[["neuron"]],neuron)
  return(fit)
}

set_order <- function(framerate){
  order <- framerate*3
  return(order)
}

find_spike_times <- function(DATA_SET, gamma_decay, L_lambda, EST_CAL, neuron, order){
  c_dfof <- correction(DATA_SET, neuron, order)
  in_noise_thr <- intrinsic_noise_threshold(c_dfof)
  fit <- estimate_spikes(dat = c_dfof, gam = gamma_decay, lambda = L_lambda, estimate_calcium = EST_CAL)
  STL = fit[[1]]
  return(STL)
}


###############################################################################################################
#FIRST FIND GAMMA DECAY
gamma_decay <- approx_gamma_decay(framerate,indicator_type)

#NEXT SET ORDER
order <- set_order(framerate)

#WE SET EST_CAL TO FALSE TO SAVE COMPUTATION TIME
EST_CAL = F

#WE CREATE VECTORIZED NEURON LIST
neuron_list <- seq(from = 1, to = num_neur, by = 1)

#NOW WE OPTIMIZE NEURON BY NEURON
#WE LOOP FOR EACH LAMBDA ON GRID
for (n in seq_along(neuron_list)){
  for (i in seq_along(LL)){
    S <- single_loop(DATA_SET, gamma_decay, LL[i], EST_CAL, neuron_list[n], order)
    SPIKES[n,i] <- S[1]
    FALSE_SPIKES[n,i] <- S[2]
  }
}

#NOW WE FIND THE OPTIMUM LAMBDA FOR EACH COMPONENT
OPT_LAMBDAS_SPIKES <- apply(FALSE_SPIKES, 1, FUN=min)
OPT_LAMBDAS_INDEX <- apply(FALSE_SPIKES, 1, FUN=which.min)
OPT_LAMBDAS <- LL[OPT_LAMBDAS_INDEX]


#NOW WE FIND SPIKE TIMES FOR EACH COMPONENT
#MAKE LIST OF LISTS FOR SPIKE TIMES
SPIKE_TIMES = list()
for (n in seq_along(neuron_list)){
  STL <- find_spike_times(DATA_SET, gamma_decay, OPT_LAMBDAS[n], EST_CAL, neuron_list[n], order)
  SPIKE_TIMES[[n]]<-STL
}

#NOW WE CONVERT TO TRAINS
SPK_TRAINS <- as.repeatedTrain(SPIKE_TIMES)

#ANALYSIS FINISHED


###############################################################################################################
#HERE WE ARE SOME PLOTTING FUNCTIONS
plot.estimated_spikes <- function(x, xlims = NULL, ...){
  if (sum(is.na(x$estimated_calcium))) {
    stop("Calcium concentration must be estimated before plotting. (Run estimate_calcium(fit).)")
  }
  ind <- 1:length(x$dat)
  rng <- range(c(x$dat, x$estimated_calcium))
  ylims <- rng 
  if (is.null(xlims)){
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "Signal", ylim = ylims, xlab = "Frame", main = paste("Neuron ", x$neuron, " Spike Inference", sep=""))
  } else {
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims, xlab = "Time")
  }
  lines(ind, x$estimated_calcium, col = "blue", lwd = 2)
  
  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "blue", lwd = 1)
  }
}
HIST_LAM <- function(OPT_LAMBDAS, OPT_LAMBDAS_SPIKES){
  #FIRST HISTOGRAM OF OPTOMIZED LAMBDAS
  hist(OPT_LAMBDAS, xlab = "OPTIMIZED LAMBDA", ylab = "FREQUENCY", main = "HISTOGRAM OF OPTIMAL LAMBDAS")
  #NEXT HISTOGRAM OF MIN LAMBDA SPIKES (SANITY CHECK)
  hist(OPT_LAMBDAS_SPIKES, xlab = "NUMBER OF FALSE SPIKES", ylab = "FREQUENCY", main = "HISTOGRAM OF FALSE SPIKES AT OPTIMAL LAMBDAS")
}
PLOT_IND_NEUR <- function(DATA_SET, gamma_decay, neuron, order, OPT_LAMBDAS, SPIKES, FALSE_SPIKES){
  #FIRST WE FIND OPTIMAL LAMBDA
  lam = OPT_LAMBDAS[neuron]
  #NOW WE RUN FIRST PLOT - THE FIT
  #NEXT WE RUN FIT WITH LAMBDA & ESTIMATE CALCIUM CONCENTRATION
  fit <- final_fit(DATA_SET, gamma_decay, lam, EST_CAL=T, neuron, order)
  #NOW WE PRINT THE RESULTS
  print(fit)
  #NOW WE PLOT THE RESULTS
  plot(fit, ylab="Corrected DFoF", xlab="Frame",  main=paste("Neuron ", neuron, " Spike Inference", sep = ""), sub="Blue = Estimated 'True' Ca++ Signal")
  #NOW WE RUN SECOND PLOT - THE OPTIMIZATION PLOT
  plot(LL,FALSE_SPIKES[neuron,], type="o", col="red", pch="o",lty=1, ylab="Spikes Detected", xlab="Lambda", main= paste("Neuron ", neuron, " Spike Inference Optimization", sep=""))
  points(LL,SPIKES[neuron,], col="blue", pch="*")
  lines(LL,SPIKES[neuron,], col="blue", lty=2)
  legend("topright", legend=c("FALSE SPIKES","SPIKES"),col=c("red","blue"),pch=c("o","*"),lty=1:2, cex=0.8, text.font=4, bg='lightblue')
}
MAKE_RASTER_PLOT <- function(SPK_TRAINS){
 plot(SPK_TRAINS, ylab = "Neuron", xlab = "Frames")
 raster(SPK_TRAINS, ylab = "Neuron", xlab = "Frames")
}

#UNCOMMENT IF DESIRED TO PRODUCE ALL PLOTS
#for(i in seq_along(neuron_list)){
 # PLOT_IND_NEUR(DATA_SET, gamma_decay, neuron=neuron_list[i], order, OPT_LAMBDAS, SPIKES, FALSE_SPIKES)
#}


    
    
    
    
    
    
    
    