# COVID-19.Mexico
Epidemiological Mathematical Models of COVID-19 applied in Mexico
This repository contains the code and input files for the analysis of a mathematical model applied to the outbreak of COVID-19 in Mexico. 
The data was obtained by the daily updates by the National Secretary of Health of Mexico https://coronavirus.gob.mx/datos/#DOView. 

We uploeded 4 files of MathLab where we performed our analysis. The files are the following: 
SEIAR_Covid_Odes_Mainland: In this file we expressed the set of differential equations used to model the spreak of the coranovirus in the country Mexico. 
SEIAR_Covid_run_Mainland: In this file we run our code to solve the differential equations, fit the parameters and obtain the graphs of the results. 
SEIAR_Covid_solver_Mainland: In this file we introduce the value of parameters obtained from our earlier work using the same set of differential equations. 
SEIAR_Covid_sse_Mainland: This file is to reproduce an estimation using the sum of squares error to obtain the parameters of our model. 

To obtain the simulation for the state-level in Mexico just modify the line 13 and 14 and insert the name of the state from the .csv files that are in this repository. 

