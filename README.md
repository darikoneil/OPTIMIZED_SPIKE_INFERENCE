# OPTIMIZED_SPIKE_INFERENCE
SCRIPTS FOR OPTIMIZED SPIKE INFERENCE USING JEWELL METHOD OF L0-OPTIMIZATION

## Optimal lambda is defined as that which minimizes the number of ("false") spikes that do not exceed twice the stdev of the intrinsic noise of an adjusted DFoF

## Constrained or different autoregressive orders can be inserted simply by adding correct tag in single_loop & final_fit functions
See following paper for further details:
Fast nonconvex deconvolution of calcium imaging data
Sean W Jewell, Toby Dylan Hocking, Paul Fearnhead, Daniela M Witten
Biostatistics, Volume 21, Issue 4, October 2020, Pages 709–726,

## Gamma decay can be solved for directly by measuring the decay of single AP in a particular genetic line.
Solving the approximate gamma = 1-(Delta/Phi) has shown sufficiency
where Delta = 1/framerate &
Phi is an approximate indicator-specific value 
(i.e., 0.7, 1.2 or 2 for fast, medium, or slow GCaMP6 respectively)

## The size of the median filter can be adjusted in the set order function
It currently is set to approximately 3 seconds to adjust for long-scale shifts in fluorescence


## This script requires the following R packages from GitHub and/or CRAN:
FastLZeroSpikeInference (https://github.com/jewellsean/FastLZeroSpikeInference)
Fractal
ifultools
Sapa
WMTSA
sPlus2R
devtools
scatterplot3D
MASS
methods
