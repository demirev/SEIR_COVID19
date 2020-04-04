library(deSolve)
library(plotly)
library(dplyr)
library(reshape)
library(htmltools)
library(dplyr)
library(purrr)
library(tidyr)
source("code/functions.R")
source("code/functions_age.R")

ageBounds = c(
  4, 7, 12, 18, 25, 30, 40, 50, 65, 75, 85, 100
)

ageContactRates <- getMossongData(
  ensureSym = T,
  ageBounds = c(
    4, 7, 12, 18, 25, 30, 40, 50, 65, 75, 85, 100
  )
)

ageContactMatrix <- ageContactRates %>%
  tidyr::spread(age_2, daily_contacts) %>%
  select(-age_1) %>%
  as.matrix()

rownames(ageContactMatrix) <- colnames(ageContactMatrix)

ageBounds <- as.numeric(rownames(ageContactMatrix))

interventions <- bind_rows(
  #setIntervention(30, 60, intervenePercentGenerator, 'b1', 0.2),
  #setIntervention(20, 40, interveneConstantGenerator, 'g2', 0.1),
  #setIntervention(40, 45, interveneAdditionGenerator, 'b3', -0.02),
  setIntervention(
    21, 260, interveneContactPercentGenerator, 
    ageBreaks = ageBounds, targetAges = c(30,85), target = 'We', change = -1
  ),
  setIntervention(
    21, 260, interveneContactPercentGenerator, 
    ageBreaks = ageBounds, targetAges = 30, target = 'We', change = -1
  )
)

resDfBase <- simInterventionsAge(
  interventions = tibble(),
  IncubPeriod= 5, #Duration of incubation period
  DurMildInf= 6, #Duration of mild infections
  FracSevere= c(rep(.05, length(ageBounds) - 1), .15), #% of symptomatic infections that are severe
  FracCritical= c(rep(.01, length(ageBounds) - 1), .5), #% of symptomatic infections that are critical
  ProbDeath= c(rep(.01, length(ageBounds) - 1), .40), #Death rate for critical infections
  DurHosp= 6, #Duration of severe infection/hospitalization
  TimeICUDeath= 8, #Duration critical infection/ICU stay
  
  AllowPresym="Yes",
  AllowAsym="Yes",
  FracAsym=.25, #Fraction of all infections that are asymptomatic
  PresymPeriod=2, #Length of infectious phase of incubation period
  DurAsym=6, #Duration of asympatomatic infection
  
  b1= 0.10, #Transmission rate (mild infections)
  b2= 0.15, #Transmission rate (severe infections)
  b3= 0.17, #Transmission rate (critical infections)
  be = 0.02, #Transmission rate (pre-symptomatic)
  b0 = 0.10, #Transmission rate (asymptomatic infections)
  Tmax= 300, #Maximum time
  InitInf= 1, #Initial # infected
  
  AllowSeason="Yes",
  seas.amp=0.0, #relative amplitude of seasonal fluctuations, in [0,1]
  seas.phase=-90,
  
  NaturalDeathRate = 0,
  InoculationSchedule = 0,
  ReinfectionChance = 0,
  
  ageBounds = as.numeric(rownames(ageContactMatrix)),
  N = rep(1000000, length(ageBounds)), #Total population size
  ExposedContacts = ageContactMatrix,
  AssymContacts = ageContactMatrix,
  MildContacts = 0.2 * ageContactMatrix,
  SevereContacts = 0.05 * ageContactMatrix,
  CritContacts = 0.01 * ageContactMatrix
)

resDfInt <- simInterventionsAge(
  interventions = interventions,
  IncubPeriod= 5, #Duration of incubation period
  DurMildInf= 6, #Duration of mild infections
  FracSevere= c(rep(.05, length(ageBounds) - 1), .15), #% of symptomatic infections that are severe
  FracCritical= c(rep(.01, length(ageBounds) - 1), .5), #% of symptomatic infections that are critical
  ProbDeath= c(rep(.01, length(ageBounds) - 1), .40), #Death rate for critical infections
  DurHosp= 6, #Duration of severe infection/hospitalization
  TimeICUDeath= 8, #Duration critical infection/ICU stay
  
  AllowPresym="Yes",
  AllowAsym="Yes",
  FracAsym=.25, #Fraction of all infections that are asymptomatic
  PresymPeriod=2, #Length of infectious phase of incubation period
  DurAsym=6, #Duration of asympatomatic infection
  
  b1= 0.10, #Transmission rate (mild infections)
  b2= 0.15, #Transmission rate (severe infections)
  b3= 0.17, #Transmission rate (critical infections)
  be = 0.02, #Transmission rate (pre-symptomatic)
  b0 = 0.10, #Transmission rate (asymptomatic infections)
  Tmax= 300, #Maximum time
  InitInf= 1, #Initial # infected
  
  AllowSeason="Yes",
  seas.amp=0.0, #relative amplitude of seasonal fluctuations, in [0,1]
  seas.phase=-90,
  
  NaturalDeathRate = 0,
  InoculationSchedule = 0,
  ReinfectionChance = 0,
  
  ageBounds = as.numeric(rownames(ageContactMatrix)),
  N = rep(1000000, length(ageBounds)), #Total population size
  ExposedContacts = ageContactMatrix,
  AssymContacts = ageContactMatrix,
  MildContacts = 0.2 * ageContactMatrix,
  SevereContacts = 0.05 * ageContactMatrix,
  CritContacts = 0.01 * ageContactMatrix
)

plotSimAggregate(resDfBase, usePlotly = F)
plotSimAggregate(resDfInt, usePlotly = F)


