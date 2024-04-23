prodynR :fish:
=====

## Package description

prodynR is an R package that allows you to simulate biomass dynamics of age-structured populations with life-history information as inputs.

## Background

The identification of Essential Fish Habitat (EFH), defined as the substrate needed by fish to complete essential stages of their life cycle, is paramount for habitat conservation strategies, and allows managers to direct monetary resources to conserve high-priority habitats. Monitoring efforts are pivotal to identify EFH. Among several monitoring metrics, NOAA considers production (the amount of biomass produced by a given species at a habitat at a given point in time) to contain the highest level of information. 

This software allows the user to estimate biomass production rates per unit of area for species with available life-history and demographic information. I obtain production rates by stochastic simulation on population dynamics as follows:

The rate of biomass production (GP - denoting Gross Production) can be represented as the product of simulated numbers-at-time and individual biomass-at-time integrated from recruitment (t = 0) to maximum age (t_{max}):

$$ GP = \int_{t = 0}^{t_{max}}\frac{\partial N}{\partial t} \frac{\partial W}{\partial t} $$

## Installation
Download the development version from GitHub:

```R
# install.packages("remotes")
remotes::install_github("matheusbarrosb/prodynR")
```


