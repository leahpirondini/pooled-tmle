## **pooled-tmle**



Targeted Maximum Likelihood Estimation for Dynamic and Static Longitudinal Marginal Structural Working Models (Peterson et al. 2014 https://www.degruyterbrill.com/document/doi/10.1515/jci-2013-0007/html).



This project contains two subfolders:



**simple\_test\_example**



Contains R code to demonstrate how to use the tmleMSM function within the ltmle package (https://cran.r-project.org/web/packages/ltmle/index.html)

* 02\_dat\_sim.R: simulates a simple dataset with two time-points: O = L0, A0, Y1, L1, A1, Y2
* 03a\_tmleMSM\_package.R: implements pooled TMLE using the tmleMSM function within the ltmle package
* 03b\_tmleMSM\_manual.R: implements pooled TMLE as above using manual R code designed to replicate the steps described in Peterson et al. 2014, section 4.2.2



**tmle\_comparison**



Contains R code to compare pooled TMLE with standard longitudinal TMLE

* 01\_dat\_sim.R: simulates longitudinal data with 20 time-points: O = L0, A0, N0, Y1, ..., L18, A18, N18, Y19, L19, A19, Y20
* 02\_standard\_ltmle.R: implements standard longitudinal TMLE using the lmtp package to target mean counterfactual outcome at each time-point under always and never treat regimes
* 03\_pooled\_tmle.R: implements pooled TMLE using the ltmle package (ltmleMSM) to target the parameters of a chosen marginal structural working model. Fitted MSM coefficients used to estimate mean counterfactual outcome at each time-point under each regime.
* 04\_true\_values.R: generates "true" counterfactual outcomes at each time-point and under each regime using a simulation-based approach by simulating under same DGP as 01\_dat\_sim.R but with A\_k set to values under intervention (for a very large dataset).

