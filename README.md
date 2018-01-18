# bootylator
Survival Estimation and Bootstrap Program for Comparative Survival Study (CSS)  
*bootylator* is a R package for estimating reach survivals and detection using a closed-form CJS methodology (ie. CJS model with time-variation). The package also contains functions to obtain confident intervals (CIs) by performing non-parametric bootstrapping, weighted or non-weighted. *bootylator* is developed mainly for CSS purposes, so it contains specific functionalities such as caculating SAR, TIR, $C_1$, $C_0$,... ect. However, users have the options to turn off CSS specific functionalities and utilize *bootylator* to only calculate survivals (and bootstrap CIs). *bootylator* package also contains a function to simulate fake capture-recapture data for study purposes.

# Installation

```{r, eval=FALSE}
install.packages("devtools")
```
