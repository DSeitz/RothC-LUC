# RothC-LUC

Welcome to the official repo for the RothC-LUC model that was tested and described in Seitz et al (2025). The model is a modification of the RothC model described in Coleman & Jenkinson (1996). RothC-LUC is implemented in R and C++ and uses the SoilR package by Sierra et al. (2012). This repo includes scripts needed to run the model, and the artificial data-set described in Seitz et al (2025) with modelling results from all 3 tested versions of RothC.

Scripts:
  - test_RothC_LUC.R : this script shows how to run the model and illustrates the differences of RothC-LUC and RothC-default
  - RothC_LUC_Model.R : R script to call the RothC-LUC model
  - RothC_LUC.cpp : c++ script of the RothC-LUC model
Data:
  - dataset-RothC-LUC.txt : compiled dataset including the main results produced in Seitz et al (2025)


References

  Coleman, K., & Jenkinson, D. (1996). RothC-26.3-A Model for the turnover of carbon in soil. In Evaluation of soil organic matter models,  Springer-Verlag, Berlin, Heidelberg. pp 237-246
  
  Seitz, D., Dechow, R., Emde, D., Schneider, F., Don, A. (2025). Improved broad-scale modelling of soil organic carbon dynamics following land-use changes. European Journal of Soil Science. DOI: 10.1111/ejss.70159
  
  Sierra, C. A., Müller, M., & Trumbore, S. E. (2012). Models of soil organic matter decomposition: the SoilR package, version 1.0. Geosci. Model Dev., 5(4), 1045-1060. https://doi.org/10.5194/gmd-5-1045-2012
