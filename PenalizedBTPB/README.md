# PenalizedBTPB
Penalized Bivariate Tensor Product B-Spline Code
This code is for the paper An Investigation of Actuarial Fair Crop Insurance Rates Using Partial Derivatives of Penalized Bivariate Tensor Product B-splines (BTPB).

There are five R scripts in this project and two simulated data frames to show how the scripts and function work.  

The scripts are: 

FunctionsBTPS.R

GraphSetUp.R

SimulationSetUP.R

Workflow_Script.R

SmallWorkingExample.R

The data frames are in folder ToyData:

BTPBmc- A 1000 rep monte carlo simulation for the penalized bivariate tensor product B-spline.  Code for this is in the Workflow_Script.

KernelMC-A 1000 rep monte carlo simulation for the Kernel Density Estimator.  Code for this is in the Workflow_Script.

########FunctionsBTPS.R###########

The FunctionsBTPS.R r script contains seven functions.  These functions are used for simulating data and and fitting data using the penalized bivariate tensor product B-spline.  There is also a function for estimating the variance of the partial derivatives of the penalized BTPB, the main focas of the paper.

The functions are: fitBTPS, yieldDataSim, coupute_bandwidth, predictDensity, VarianceEstimator, kernelEstimator, and fitKernel.

########fitBTPS

fitBTPS is a function used for fitting a penalized bivariate tensor product B-spline. In addition to fitting the B-spline, which uses the function gam in the mgcv package, it returns estimates, first partial derivative wrt the x value estimates, and second partial derivative wrt the x value estimate for a given set of points.  Knots, degree of the function and degree of the penalty can all be adjusted.  the knots represent the number of interior knots. The residuals of the penalized BTPB are also fitted in this function in order to save a step when getting variance estimates in the function VarianceEstimator.  

########yieldDataSim

This function simulates data according to the simulation study specified in chapter 4.  This function simulates coverage rate (x_i), APH (average historical/actual production history) yield (z_i), current yield (w_i), premium values (mu(x,z)), and observed premium values (y_i). There are three main components, the mean structure (mu(z_i)), the error structure (epsilon_i), and the variance structure (sigma) used to estimate the current yield.  The mean structure is either Linear or Quad, Linear: mu(z_i)=-25+1.3z_i, Quad: mu(z_i)=1.2(z_i-50)+(z_i-150)^2/200. The error structure is Normal or Beta. Normal:10*N(0,1), Beta: 50*(Beta(5,3)-5/8). The variance structure is NonConst or Const. NonConst: abs(z_i)^.2 Const: 10 if error structure in Normal or 50 if error structure is Beta.


########compute_bandwidth

The compute_bandwidth function is a helper function for predictDensity. It is used for calculating bandwidths for bivariate kernal density estimation.   

########predictDensity

The predictDensity function is a helper function for VarianceEstimator.  It is used for calculating the bivariate kernel density for f_hat.  

########VarianceEstimator

This function calculates the variance for any partial derivative up to the fourth using the penalized BTPB.  This function uses the object returned by the fitBTPS function.  The equation used is (22) from the paper.

########kernelDensity

The kernelDensity function is used to estimate the desnity for current year yield data.  It is a helper function for the fitKernel function.  The step by step explination of the steps can be found in the supplimental file section 1. To estimate the kernel density, a Gaussian smoothing kernel is used with bandwidth prescribed by data-based selection used in SHeather and Jones (1991).  This is the default given in r.  

########fitKernel

The fitKernel function estimates the kernel density, estimates the value at new locations, and estimates the variance of the estimated points using a delete one Jackknife method.  

########GraphSetUp.R###########

The GraphSetUp.R script has two functions, one for penalized BTPB fits (myGraphs) and one for kernel density fits (myGraphsKern).  These functions work with monte carlo simulations and gives graphical results for these simulations.  Each function produces two graphs, a line graph showing the 2.5 and 97.5 percentiles and mean of the monte carlo replicates compared to the true value of the function, and a heat map showing the coverage probability of the variance estimators for the penalized BTPB fit and kernel density fit.  myGraphs can be used to produce figures 1 and 2 in the paper while myGraphsKern can be used to produce figures 1 and 2 in the supplimental file.

########myGraphs

This function will create graphs at given quantile values.  It takes output from the simStudyFunc in the SimulationSetUp.R script. It assumes that the first replicate of the estimates, (matrix[ , ,1]) are the true values. It then uses these true values and the rest of the MC replicates to a line graph showing the 2.5 and 97.5 percentiles and mean of the monte carlo replicates compared to the true value of the function.  Then it uses the variance matrix of the penalized BTPB from the simStudyFunc to create a heat map of the coverage rate.  95% is the nominal coverage in this case.  Red is for values above 95%, white is for values near 95% and blue is for values below 95%.  The blue get darker the farther away the coverage rate is from 95%


########myGraphsKern

This function will create graphs at given quantile values.  It takes output from the simStudyKernel in the SimulationSetUp.R script. It assumes that the first replicate of the estimates, (matrix[ , ,1]) are the true values. It then uses these true values and the rest of the MC replicates to a line graph showing the 2.5 and 97.5 percentiles and mean of the monte carlo replicates compared to the true value of the function.  Then it uses the variance matrix of the kernel density from the simStudyFunc to create a heat map of the coverage rate.  The yield value have to be standardized in this case since the only values the w_i values where there is simulated data changes depending on the z_i value.  Standardizing the values allows for us to keep the square set up of the heat map.   95% is the nominal coverage in this case.  Red is for values above 95%, white is for values near 95% and blue is for values below 95%.  The blue get darker the farther away the coverage rate is from 95%.

########SimulationSetUp.R###########

The SimulationSetUp.R script provides a way to perform monte carlo studies for the penalized BTPB and the kernel density estimator.  There are two functions in this script, simStudyFunc and simStudyKernel, for the penalized BTPB and the kernel density estimator respectively.  Both functions allow for their output to be used in their respective functions in the GraphSetUp.R script.  

########simStudyFunc

Creates a MC simulation for any of the sim study combinations discussed in Section 4 as well as a constant variance set up. Fits will be performed using the penalized BTPB. Will create two three dimensional arrays returned in a list for the estimates and the variance estimates.  Dimension 1 is the coverage rate values, dimension 2 is the APH yield values, and dimension three is the replicates.  For the estimates est[ , ,1] will have the true values based on the simulation set up.  Variances are calculated according to (22), so the z^2 and p have to be divided off later.  

########simStudyKernel

Creates a MC simulation for any of the sim study combinations discussed in Section 4 as well as a constant variance set up. Fits will be performed using the kernel density estimation. The function will create two three dimensional arrays returned in a list for the estimates and the variance estimates.  Dimension 1 is the coverage rate values, dimension 2 is the APH yield values, and dimension three is the replicates.  For the estimates est[ , ,1] will have the true values based on the simulation set up.  Variances are calculated using delete one jackknife method, so this function may take a while to run and takes longer the larger the sample size.

########SmallWorkingExample.R###########

This script doesn't contain any functions, but uses fitBTPS and VarianceEstimator from FunctionsBTPS.R to show how the penalized BTPB works on a very simple simulated data, one in a quadratic form.  It is used more as a proof of concept.  

x and z are simulated from 1:100 and y is defined as y=x^2+z^2+z*x^3. The derivatives will be taken wrt x.  So the first derivative is y'=2x+3z*x^2 and the second derivative is y"=2+6zx.  A sample of 500 is fitted by the penalized BTPB.  3-D plots are then used to compare the fitted values to the true values of the actual fit and the first two partial derivatives.  Next the variance is calculated and applied to the second partial derivative wrt x and a 3-D plot with a 95% confidence interval is shown.

########Workflow_Script.R###########

The Workflow_Script.R shows how to use the functions in the other scripts and how to create the graphs shown in the paper.

First data is simulated (using function yieldDataSim from FunctionsBTPS.R).  Then data is fit using the penalized BTPB (using function fitBTPS from FunctionsBTPS.R).  Then the variance is estimated (using function VarianceEstimator from FunctionsBTPS.R). This shows a basic fit of the penalized BTPB.

Next a set up of the monte carlo simulations is shown.  These are used to create the toy data used for creating graphs 1&2 in the paper and supplimental file.  The arguments used are in the comments. Also in the comments are the arguments used when running the script for the paper.  

Using the MC toy data, figures 1&2 from the paper and supplimental file are produced.

Next data is simulated to take the place of the empirical data in order to show how the graphs were created for the paper.

The code for graphs 3 through 6 are not in function and are written out in Workflow_Script.R.  Graph 3 show the comparison of the penalized BTPB fit and the kernel density fit.  Graph 4 shows the whole estimated kernel density curve.  Graph 5 and 6 shows a comparison of the real premium values (or subsidized values) with the fitted BTPB

###################################How to create the Graphs and Table in the Paper################################

Figures 1&2 in the paper can be produced by using functions simStudyFunc in the SimulationSetUp.R script and myGraphs in the GraphSetUp.  Arguments needed for reproducing these figures can be found in the comments of the Workflow_Script.R script.

Figures 1&2 in the supplemental paper can be produced by using functions simStudyKernel in the SimulationSetUp.R script and myGraphsKern in the GraphSetUp.  Arguments needed for reproducing these figures can be found in the comments of the Workflow_Script.R script.

Figure 3 is produced by fitting the data with fitBTPS, VarianceEstimator, and fitKernel functions in FunctionsBTPS.R script.  Code for the graph is provided in Workflow_Script.R script.

Figure 4 is produced by fitting the data with the fitKernel function in FunctionsBTPS.R script.  Code for the graph is provided in Workflow_Script.R script.

Figures 5&6 is produced by fitting the data with the fitBTPS function in FunctionsBTPS.R script.  Code for the graph is provided in Workflow_Script.R script.

Table 2 is produced by fitting the data with the fitBTPS function in FunctionsBTPS.R script.  Code for the graph is provided in Workflow_Script.R script.


###################################Additional R packages needed################################

mgcv_1.7-26

msm_1.5

rgl_0.95.1201

MASS_7.3-29

reshape2_1.2.2

plyr_1.8

ggplot2_0.9.3.1

nlme_3.1-111
