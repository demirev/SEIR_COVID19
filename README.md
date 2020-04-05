# SEIR_COVID19

The code for creating the R Shiny application https://demirev.shinyapps.io/SIRinterventions/ which is based on https://alhill.shinyapps.io/COVID19seir/. 

The R shiny code is in directory COVID19seir. 

This code is an SEIR model for COVID-19 infection, including different age groups and the option to play out multiple interventions. 

The code that produces the interface and functionality of the Shiny App is in files
* **server.R**
* **ui.R**

Files used in the explanatory sections of the app are
* **SEIR.Rmd**
* **About.Rmd**
* **Tutorial.Rmd**

The functions that actually run this variety of the model and process the parameters are in the **code/functions_age.R** file

If you want to run the code to produce the same outputs as Shiny but without dealing with the app structure, you can use the R scripts
* **runSpreadByAge.R**
