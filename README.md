The project is for Lan Xue, Xinxin Shu, Peibei Shi, Colin O. Wu, Annie Qu. Time-varying feature selection for longitudinal analysis (2019).

1. locVcof_0.1.0.zip is the R package compiled under R 3.6.1, including three major functions

a) The function "splinevcf" is to calculate predicted y, time-varying coefficients and spline coefficients for the proposed model in the paper.

b) The function "psvcf" is to calculate predicted y, time-varying coefficients and spline coefficients for the varying-coefficient model using the spline method WITHOUT penalty term.

c) The function "lambdasel_BIC" is to calucate BIC values and select the best tunning parameters.

You can use ?"function name" after the package installed to get help and details for each function 



2. test_simulation.R is a full example to use the R package and run a simulation. 
