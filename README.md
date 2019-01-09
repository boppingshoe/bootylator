# bootylator

Survival Estimation and Bootstrap Program for Comparative Survival Study (CSS)  
*bootylator* is an R package for estimating reach survivals and detection using a closed-form CJS methodology (ie. CJS model with time-variation). The package also contains functions to obtain confident intervals (CIs) by performing non-parametric bootstrapping, weighted or non-weighted. *bootylator* is developed mainly for CSS purposes, so it contains specific functionalities such as calculating SAR, TIR, C1, C0,... ect. However, users have the options to turn off CSS specific functionalities and utilize *bootylator* to only calculate survivals (and bootstrap CIs). *bootylator* package also contains a function to simulate fake capture-recapture data for study purposes.

# Installation

*bootylator* can be installed on your computer directly from GitHub, and Hadley Wickham's *devtools* package can help you do that:
```{r, eval=FALSE}
install.packages("devtools")
```

If you can't install *devtools*, it's probably that you need *Rtools* on your computer first (https://cran.r-project.org/bin/windows/Rtools/). After installing *Rtools* (if you haven't already), you can now install *bootylator* from GitHub:
```{r, eval=FALSE}
devtools::install_github("boppingshoe/bootylator", build_vignettes = TRUE)
```

For further information on survival estimation and bootstrap methods for CSS:  
McCann, J., B. Chockley, E. Cooper, B. Hsu, H. Schaller, S. Haeseker, R. Lessard, C. Petrosky, T. Copeland, E. Tinus, E. Van Dyke, A. Storch and D. Rawding. 2017. Comparative survival study of PIT-tagged spring/summer/fall Chinook, summer steelhead, and sockeye, 2017 Annual Report. Project No. 19960200.
