# CDWQ-E2
Master's Thesis: Development of the Expanded Comprehensive Disinfection and Water Quality Model

## Introduction
This repository contains my research report submitted to the Faculty of Engineering and the Built Environment at the University of the Witwatersrand, for the degree of Master of Science in Engineering. The MATLAB source code and all input and output files are also provided. <br> <br> The aim of the project was to develop a mathematical model to help researchers and engineers better understand the complex set of biological, chemical and hydraulic factors that affect drinking water quality, after treatment. <br> <br>
A batch version of the model was first developed, which does not account for hydraulic processes, in order to demonstrate the impact of biological and chemical processes on water quality. <br>
Following this, a distribution system version of the model was developed, building on the previous model, which accounts for hydraulic processes. <br>

## Research Report
* The research report is located in `CDWQ-E2-Research-Report.pdf`.

## Batch Version
* The file `CDWQ-E2-Batch.m` is the batch version of model. A user is required to specify inputs within the code and the outputs of the various simulations described in the report are provided: <br> `Batch-Outputs-Varying-CtCO3.mat` <br> `Batch-Outputs-Varying-pH.mat`<br> `Batch-Outputs-Varying-Temperature.mat`

## Distribution System Version
* The file `CDWQ-E2-Distribution.m` is the distribution system version of the model. This code requires the function `Advection.m`.
* A user is prompted to select an input .xlsx or .xls file. The input files for the various simulations described in the report are provided: <br>
`Baseline-Inputs.xlsx` <br>
`Alternative-1-Inputs.xlsx` <br>
`Alternative-2-Inputs.xlsx` <br>
`Alternative-3-Inputs.xlsx` <br>
`Alternative-4-Inputs.xlsx` <br>
`Alternative-5-Inputs.xlsx` <br>
`Alternative-6-Inputs.xlsx` <br>
`Alternative-7-Inputs.xlsx` <br>
`Alternative-8-Inputs.xlsx` <br>
* The outputs of the various simulations described in the report are provided: <br>
`Baseline-Outputs.mat` <br>
`Alternative-1.mat` <br>
`Alternative-2.mat` <br>
`Alternative-3.mat` <br>
`Alternative-4.mat` <br>
`Alternative-5.mat` <br>
`Alternative-6.mat` <br>
`Alternative-7.mat` <br>
`Alternative-8.mat` <br>

