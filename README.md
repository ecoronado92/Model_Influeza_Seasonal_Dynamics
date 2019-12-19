# Hidden Markov Models (HMM) and Autoregressive Processes (AR)


### Summary

Early detection of the influenza outbreaks is one of the biggest challenges of outbreak surveillance systems. In this paper, a finite, homogeneous two-state Hidden Markov Model (HMM) was developed to determine the epidemic and non-epidemic dynamics influenza-like illnesses (ILI) in a differenced time series. These dynamics were further modeled via a first-order auto-regressive process (AR1) based on the state, and parameters estimated via the Baum-Welch
algorithm. The model was evaluated with US ILI data from 1998-2018.

### Files
- `*.html` final document
- `*.Rmd` generates final document
- `BW_fcns.R` is a script containing the Baum-Welch functions to train the model
- `viterbi.`R is a script containing the decoding algorithm to determin epidemic vs non-epidemic weeks based on training data
- `flu_paper_model.R` is a JAGS model from a similar paper run as a comparison