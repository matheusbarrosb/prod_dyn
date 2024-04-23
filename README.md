prodynR :fish:
=====

## Package description

prodynR is an R package that allows you to simulate biomass dynamics of age-structured populations with life-history information as inputs.

## Background

The identification of Essential Fish Habitat (EFH), defined as the substrate needed by fish to complete essential stages of their life cycle, is paramount for habitat conservation strategies, and allows managers to direct monetary resources to conserve high-priority habitats. Monitoring efforts are pivotal to identify EFH. Among several monitoring metrics, NOAA considers production (the amount of biomass produced by a given species at a habitat at a given point in time) to contain the highest level of information (Levin et al. 2005). 

This software allows the user to estimate biomass production rates per unit of area for species with available life-history and demographic information. I obtain production rates by stochastic simulation on population dynamics using the following system of differential equations:

The rate of biomass production (GP - denoting Gross Production) can be represented as the product of simulated numbers-at-time and individual biomass-at-time integrated from recruitment (t = 0) to maximum age (t<sub>max</sub>):

$$ GP = \int_{t = 0}^{t_{max}}\frac{\partial N}{\partial t} \frac{\partial W}{\partial t} $$

Numbers-at-time denote the decline in numbers of a cohort through time according to a daily mortality rate:

$$ \frac{\partial N}{\partial t} = N_{t-1}e^{-M_{t}}$$

The mortality rate can be empirical (i.e. field-based, or borrowed from the literature) or estimated from equations from the literature. For instance, we could use an equation that assumes mortality rates are inversely-related to size (Lorenzen 2022):

$$ M_{t} = M_{ref}(\frac{L_{ref}}{L_{t}}) $$

where M<sub>ref</sub> is a reference mortality rate for a given size L<sub>t</sub>.

Individual biomass-at-time is estimated by applying parameters of a length-weight relationship to the predicted size-at-time:

$$ \frac{\partial W}{\partial t} = \begin{cases}
            \alpha(L_{t-1} + GR)^\beta & \text{if L<sub>t</sub> < L<sub>m</sub>} \\
            \alpha[L_{\infty}e^{-K(t - t_{0})}]^\beta & \text{otherwise}
        \end{cases} $$

where GR is a linear absolute growth rate in mm/day, and L<sub>m</sub> is the size at maturity.

## Installation
Download the development version from GitHub:

```R
# install.packages("remotes")
remotes::install_github("matheusbarrosb/prodynR")
```

## Usage

Please refer to function documentations (e.g. ?pred_B, ?stoch_P, and ?sensitivity for examples).

## References

Levin, P. S., & Stunz, G. W. (2005). Habitat triage for exploited fishes: Can we identify essential “Essential Fish Habitat?”. Estuarine, Coastal and Shelf Science, 64(1), 70-78.

Lorenzen, K. (2022). Size-and age-dependent natural mortality in fish populations: Biology, models, implications, and a generalized length-inverse mortality paradigm. Fisheries Research, 255, 106454.
