# Analysis of optical transmitance in FTIR

This code is intended to analyse two different sets of measurements where we study the frustrated total internal reflection (FTIR).

The measurements were taken with a system formed by two coupled triangular prisms separated by a distance of microns:

<p align="center">
  <img src="https://github.com/carolFR13/optical-transmitance-study/blob/main/data/img/prisms_scheme.png" width="450">
</p>

The source used was a laser with a wide spectral range. We took mesurements of the transmitance using a He-Ne filter varying the external angle 
as well as mesurements in all the spectral range for different fixed external angles.

Our final goal is to compare the experimental transmission with the well-known thoretical expression for the FTIR. In order to do that, we need 
to know both the distance between the prisms and the prism's angle (the $\alpha$ angle in the scheme). The code provided here 
makes the study from the experimental data to obain those parameters.
