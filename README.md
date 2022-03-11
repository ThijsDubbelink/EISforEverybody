# EISforEverybody
This repository contains the following files.
- Import.py
- TR_DRT.py
- calcrq.py
- calculate_EIS.py
- MAINMAIN file

# Import.py
creates import object with:
- list of filenames self.LOF
- list of dataframes self.DFL 

A .mpt biologic file is used for importing.
Furthermore it filters the cyclenumbers and specifies whether a row contains an EIS measurement.

# TR_DRT.py
- Tikhonov regularization is used to obtain the DRT.
- Uses a circuit like -l-R_inf-(DRT)-C to fit the EIS data. 
- l has a fixed value specified in the MAINMAIN file. e.g. 1e-7
- outputs of TR_DRT are:
1. R_inf
2. DRT (gammaRC)
3. C capacitance

Physical interpretation of the following components may be interpreted as:
1. l is the inductance of the wires.
2. R_inf electrolyte + electrical resistance.
3. DRT (Distribution of Relaxation Times) Could contain (SEI,CEI,Charge transfer,diffusive properties).
4. C Capacitive effects since it is a finite length system.

TR_DRT.py was inspired by the DRT libraries of https://github.com/ciuccislab

# calcrq.py
- Use scipy.optimize's curve_fit to Calculate R, time constant τ0, φ from the obtained DRT (gammaRC) of TR_DRT.py according to equation 4 of the following paper:
Distribution (function) of relaxation times, successor to complex nonlinear least squares analysis of electrochemical impedance spectroscopy? by Bernard A Boukamp et al.

# calculate_EIS.py
- Calculate Z_real and Z_imag from gammaRC.
- this function can be used to visually inspect goodness of fit.

# MAINMAIN[std] at ocv.py

This is an example how the above mentioned files can be used to analyze a dataset of EIS data.

# DRT Plot
![Qualitative](https://user-images.githubusercontent.com/70464197/157878490-ea76fc06-cdfc-4f41-a9be-d99c5df5e145.png)
- left figure: The solid dots display the raw data, the black crosses are the calculated Impedance from the DRT.
- right figure: The DRT profile for different batteries.

# rq PLOT
![Cyc=1 0,V=3 4652586](https://user-images.githubusercontent.com/70464197/157879133-5172f222-be39-4ccb-9255-b90682d0b966.png)
- The colored lines are the individual fits of equation 4 to the DRT pattern.
- The dashed line shows the sum of the individual fits.
- The solid dots is the DRT data.

# csv output file fitted parameters
- furthermore the MAINMAIN script will output a csvfile with all the parameters and corresponding errors from the DRT and R_inf.
[EISAnalysis.csv](https://github.com/ThijsDubbelink/EISforEverybody/files/8233162/EISAnalysis.csv)
