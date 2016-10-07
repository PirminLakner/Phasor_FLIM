# Phasor_FLIM
Phasor analysis of FLIM data with MATLAB script

## Overview

This repository contains the source codes and raw data that were used in the project:

"A phasor approach analysis of multiphoton FLIM measurements of three-dimensional cell culture models" by Pirmin H. Lakner, Michael G. Monaghan, Yvonne Möller, Monilola A. Olayioye, Katja Schenke-Layland.

This project will be published in Scientific Reports. Additional information is also available in the article supplements. The source code was tested in Matlab R2016a without toolboxes.

### Short Manual

The MATLAB script is the 'phasor.m' file. It must be executed with MATLAB. Copy phasor.m file to the MATLAB folder in e.g. 'Documents', open MATLAB, open phasor.m with MATLAB and press 'run' (F5). 

Meaning of parameters need to be set before:

'thresh' is threshold of pixels considered in calculations based on their intensity compared to the brightest pixel. Inensity of a pixel is the sum of photons over all time bins. Values reach from 0 to 1 standing for 0% to 100%. 0.0 considers all pixels, 1.0 considers only brightest pixel. So 0.3 considers all pixel with 'brightness' above 30% of the brightest pixel.

'harmonic' is the harmonic number of the fourier tranform. Classical values are 0, 1, 2, ... Attention: Code is writen that in all calculations only the first period is considered! So for harmonic number 2 ony the first half of all time bins are considered in lifetime calculations! So the number of contributing time bins is given by number of total time channels divided by harmonic number. Could be very useful for imcomplete decays, e.g. the detector does not record decay information over the whole laser repition period. 

'singleorbatch' is the switch for single or batch calculations. Options are 'single' or 'batch'. When chosed 'single' one file can be selected for analysis. For 'batch' a whole folder gets analysed, file for file. Every file in the folder needs to be a valid .asc-file exported with SPCimage software from Becker&Hickl.

'lifetimefitting' is a switch for type of lifetime calculations. Valid options are 'true' and 'false'. For 'true' all phasor points are fitted with a linear to calculate both contributing lifetimes from intersections with universal circle. 'false' suppress the fitting and calculates the lifetime of the mean point of all phasor data. This method is only valid if the mean point lays on the universal circle. Applicable for lifetime calculations of a mono-exponential decay, e.g. for coumarin(2.5ns). 

'shift_bins' is the parameter for the number of time bins after the maximum that should be ignored for calculations. Could improve results for measurements with bad IRF or flurophores with long lifetimes to suppress influence if IRF. Parameter need to be set in intereger numbers like 0,1,2,... Standard is 0

'freq0' is the laser pulse frequency in Hz, e.g. 8E+7 for 80 MHz.

'delta_t' is the width of one time channel, e.g. 1.5625E-10 for 156.25 ps. The total number of time channels per period is calculated from the frequency and the time channel width.

###Example from paper

The results for the coumarin measurement shown in the paper mentioned above are calculated with the following parameters:

% -- set constants

thresh = 0.3;               % threshold [0,1]

harmonic = 4;               % harmonic number (1,2,3,4,...)

singleorbatch = 'single';   % 'single' or 'batch' processing

lifetimefitting = 'false';   % 'true' for linear fitting of phasor points for lifetime

% 'false' for calculation of lifetime with mean of phasor points


shift_bins = 4;             % number of bins shifted to right from max for minimizing effect of IRF

freq0 = 8E+7;                   % laser frequency, here 80 MHz

delta_t = 1.5625E-10;           % width of one time channel

## Licenses

The License for the source codes is available in the LICENSE file. 

## Contact

Please feel free to contact Pirmin Lakner (pirmin.lakner@googlemail.com).
