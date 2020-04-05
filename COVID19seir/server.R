library(deSolve)
library(reshape)
library(plotly)
library(dplyr)
library(htmltools)
library(dplyr)
library(purrr)
library(tidyr)
library(shinydashboard)
source("code/functions.R")
source("code/functions_age.R")

ageBounds = c(
  "< 4" = 3, 
  "4 - 7" = 7, 
  "8 - 12" = 12,
  "13 - 18" = 18, 
  "19 - 25" = 25, 
  "26 - 35" = 35, 
  "36 - 50" = 50, 
  "51 - 65" = 65, 
  "66 - 75" = 75, 
  "76 - 85" = 85, 
  "> 85" = 100
)

ageContactRates <- getMossongData(
  ensureSym = T,
  ageBounds = as.numeric(ageBounds)
) %>%
  mutate(
    age_both = map2_chr(
      age_1, age_2, 
      function(a1, a2) paste(sort(c(a1, a2)), collapse = "_")
    )
  )

ageContactRatesPairs <- ageContactRates %>%
  distinct(age_both, .keep_all = T)


function(input, output, session) {

  ageContacts <- map(seq_len(nrow(ageContactRatesPairs)/3), function(i) {
    fluidRow(
      column(
        4, 
        numericInput(
          paste0("A",ageContactRatesPairs$age_both[(i-1)*3+1]),
          paste0(
            names(ageBounds)[ageBounds == ageContactRatesPairs$age_1[(i-1)*3+1]],
            " and ",
            names(ageBounds)[ageBounds == ageContactRatesPairs$age_2[(i-1)*3+1]]
          ),
          value = round(ageContactRatesPairs$daily_contacts[(i-1)*3+1], 2),
          min = 0, max = 50, step = 0.01
        )
      ),
      column(
        4,
        numericInput(
          paste0("A",ageContactRatesPairs$age_both[(i-1)*3+2]),
          paste0(
            names(ageBounds)[ageBounds == ageContactRatesPairs$age_1[(i-1)*3+2]],
            " and ",
            names(ageBounds)[ageBounds == ageContactRatesPairs$age_2[(i-1)*3+2]]
          ),
          value = round(ageContactRatesPairs$daily_contacts[(i-1)*3+1], 2),
          min = 0, max = 50, step = 0.01
        )
      ),
      column(
        4,
        numericInput(
          paste0("A",ageContactRatesPairs$age_both[(i-1)*3+3]),
          paste0(
            names(ageBounds)[ageBounds == ageContactRatesPairs$age_1[(i-1)*3+3]],
            " and ",
            names(ageBounds)[ageBounds == ageContactRatesPairs$age_2[(i-1)*3+3]]
          ),
          value = round(ageContactRatesPairs$daily_contacts[(i-1)*3+3], 2),
          min = 0, max = 50, step = 0.01
        )
      )
    )
  })
  
  output$ageContacts <- renderUI({ageContacts})
  
  contactMatrix <- reactive({
    ageContactRates %>%
      mutate(daily_contacts = map_dbl(age_both, function(ageKey) {
        res <- input[[paste0("A",ageKey)]]
        if (!length(res)) {
          return(
            ageContactRatesPairs$daily_contacts[ageContactRatesPairs$age_both == ageKey]
          )
        } 
      })) %>%
      select(-age_both) %>%
      tidyr::spread(age_2, daily_contacts) %>%
      select(-age_1) %>%
      as.matrix()
  })
  
  # Plot timecourse of all variables
  resultBase <- reactive({
    simInterventionsAge(
      interventions = tibble(),
      IncubPeriod = input$IncubPeriod,
      DurMildInf = input$DurMildInf, 
      FracSevere = c(
        input$FracSevere1,
        input$FracSevere2,
        input$FracSevere3,
        input$FracSevere4,
        input$FracSevere5,
        input$FracSevere6,
        input$FracSevere7,
        input$FracSevere8,
        input$FracSevere9,
        input$FracSevere10,
        input$FracSevere11
      ),
      FracCritical = c(
        input$FracCritical1,
        input$FracCritical2,
        input$FracCritical3,
        input$FracCritical4,
        input$FracCritical5,
        input$FracCritical6,
        input$FracCritical7,
        input$FracCritical8,
        input$FracCritical9,
        input$FracCritical10,
        input$FracCritical11
      ),
      ProbDeath = c(
        input$ProbDeath1,
        input$ProbDeath2,
        input$ProbDeath3,
        input$ProbDeath4,
        input$ProbDeath5,
        input$ProbDeath6,
        input$ProbDeath7,
        input$ProbDeath8,
        input$ProbDeath9,
        input$ProbDeath10,
        input$ProbDeath11
      ),
      DurHosp = input$DurHosp,
      TimeICUDeath = input$TimeICUDeath,
      AllowPresym = input$AllowPresym,
      AllowAsym = input$AllowAsym,
      FracAsym = input$FracAsym,
      PresymPeriod = input$PresymPeriod, #Length of infectious phase of incubation period
      DurAsym = input$DurAsym,
      
      b1 = input$b1,
      b2 = input$b2,
      b3 = input$b3,
      be = input$be,
      b0 = input$b0,
      Tmax = input$Tmax, 
      InitInf= input$InitInf,
      
      AllowSeason = input$AllowSeason,
      seas.amp = input$seas.amp,
      seas.phase = input$seas.phase,
      
      NaturalDeathRate = input$NaturalDeath,
      InoculationSchedule = 0,
      ReinfectionChance = input$LoseImunity,
      
      ageBounds = ageBounds,
      N = c(
        input$N1,
        input$N2,
        input$N3,
        input$N4,
        input$N5,
        input$N6,
        input$N7,
        input$N8,
        input$N9,
        input$N10,
        input$N11
      ),
      ExposedContacts = contactMatrix(),
      AssymContacts = contactMatrix(),
      MildContacts = input$SocConMild * contactMatrix(),
      SevereContacts = input$SocConSevere * contactMatrix(),
      CritContacts = input$SocConCritical * contactMatrix()
    )
  })
  
  resultInterventions <- reactive({
    
    if (length(interventionInfo$intrv)) {
      interventions = interventionInfo$intrv %>%
        map(~.$intervention) %>%
        reduce(bind_rows)
    } else {
      interventions <- tibble()
    }
    
    simInterventionsAge(
      interventions = interventions,
      IncubPeriod = input$IncubPeriod,
      DurMildInf = input$DurMildInf, 
      FracSevere = c(
        input$FracSevere1,
        input$FracSevere2,
        input$FracSevere3,
        input$FracSevere4,
        input$FracSevere5,
        input$FracSevere6,
        input$FracSevere7,
        input$FracSevere8,
        input$FracSevere9,
        input$FracSevere10,
        input$FracSevere11
      ),
      FracCritical = c(
        input$FracCritical1,
        input$FracCritical2,
        input$FracCritical3,
        input$FracCritical4,
        input$FracCritical5,
        input$FracCritical6,
        input$FracCritical7,
        input$FracCritical8,
        input$FracCritical9,
        input$FracCritical10,
        input$FracCritical11
      ),
      ProbDeath = c(
        input$ProbDeath1,
        input$ProbDeath2,
        input$ProbDeath3,
        input$ProbDeath4,
        input$ProbDeath5,
        input$ProbDeath6,
        input$ProbDeath7,
        input$ProbDeath8,
        input$ProbDeath9,
        input$ProbDeath10,
        input$ProbDeath11
      ),
      DurHosp = input$DurHosp,
      TimeICUDeath = input$TimeICUDeath,
      AllowPresym = input$AllowPresym,
      AllowAsym = input$AllowAsym,
      FracAsym = input$FracAsym,
      PresymPeriod = input$PresymPeriod, #Length of infectious phase of incubation period
      DurAsym = input$DurAsym,
      
      b1 = input$b1,
      b2 = input$b2,
      b3 = input$b3,
      be = input$be,
      b0 = input$b0,
      Tmax = input$Tmax, 
      InitInf= input$InitInf,
      
      AllowSeason = input$AllowSeason,
      seas.amp = input$seas.amp,
      seas.phase = input$seas.phase,
      
      NaturalDeathRate = input$NaturalDeath,
      InoculationSchedule = 0,
      ReinfectionChance = input$LoseImunity,
      
      ageBounds = ageBounds,
      N = c(
        input$N1,
        input$N2,
        input$N3,
        input$N4,
        input$N5,
        input$N6,
        input$N7,
        input$N8,
        input$N9,
        input$N10,
        input$N11
      ),
      ExposedContacts = contactMatrix(),
      AssymContacts = contactMatrix(),
      MildContacts = input$SocConMild * contactMatrix(),
      SevereContacts = input$SocConSevere * contactMatrix(),
      CritContacts = input$SocConCritical * contactMatrix()
    )
  })
  
  plotIandEMax <- reactive({
    getMaxValEandI(resultBase())
  })
  
  output$plotBase <- renderPlotly({
    plotSimAggregate(resultBase())
  })
  output$plotBaseInttab <- renderPlotly({
    plotSimAggregate(resultBase())
  })
  output$plotIntInttab <- renderPlotly({
    plotSimAggregate(resultInterventions())
  })
  
  output$plotEandI <- renderPlotly({
    plotSimEandI(resultBase())
  })
  output$plotEandIBaseInttab <- renderPlotly({
    plotSimEandI(resultBase(), maxY = plotIandEMax() * 1.05)
  })
  output$plotEandIIntInttab <- renderPlotly({
    plotSimEandI(resultInterventions(), maxY = plotIandEMax() * 1.05)
  })
  
  output$plotAge <- renderPlotly({
    plotSimByAge(
      resultBase(), 
      whichGroup = input$agePlotWhat, 
      usePlotly = F, 
      ageBounds = ageBounds
    )
  })
  output$plotAgeBaseInttab <- renderPlotly({
    plotSimByAge(
      resultBase(), 
      whichGroup = input$agePlotWhatInt, 
      usePlotly = F, 
      ageBounds = ageBounds
    )
  })
  output$plotAgeIntInttab <- renderPlotly({
    plotSimByAge(
      resultInterventions(), 
      whichGroup = input$agePlotWhatInt, 
      usePlotly = F, 
      ageBounds = ageBounds
    )
  })
  
  countIntrv   <- reactiveValues(n=0) # will hold a counter of number of blocks
  interventionInfo   <- reactiveValues(intrv = list()) # will hold block info
  selectedIntrv <- reactiveValues(which = 0) # will hold id of selected
  
  #increment counter on button click
  observeEvent(input$addIntrv, {
    if (input$newIntrv != "" & countIntrv$n < 24 &
        !(input$newIntrv %in% sapply(interventionInfo$blks, function(b) b$name))) {
      
      countIntrv$n <- countIntrv$n + 1
      
      entry <- list(
        name = input$newIntrv
      )
      
      interventionInfo$intrv[[length(interventionInfo$intrv)+1]] <- entry
    }
    
    #input$newBlk <- ""
    
  })
  
  #dynamically add UI inputs based on button clicks
  interventionHeaders <- reactive({
    n <- countIntrv$n
    selectedIntrv$which <- n
    if (n > 0){
      lapply(seq_len(n), function(i) {
        fluidRow(
          column(
            2, style = "margin-top: 20px;", #should be relative
            actionButton(
              inputId = paste0("remIntrv",i), 
              label = "", 
              icon = icon("times")
            )     
          ),
          column(
            10, style = "margin-top: 20px;",
            actionButton(
              inputId = paste0("selIntrv",i), 
              label = interventionInfo$intrv[[i]]$name, 
              width = '100%'
            )
          )
        )
      })
    }
  })
  
  output$interventionHeaders <- renderUI({interventionHeaders()})
  
  #dynamically change selected block box
  interventionInputs <- reactive({
    id <- selectedIntrv$which
    hasInfo <- if (id == 0) F else !is.null(interventionInfo$intrv[[id]]$inttype)
    if (id != 0){
      shinydashboard::box(
        title = interventionInfo$intrv[[id]]$name,
        width = 12, collapsible = F, collapsed = F,
        column(
          6,
          selectInput(
            "inttype",
            "Intervention Type",
            choices = c(
              "Restrict All Contacts For All Ages" = "allContacts", 
              "Restrict All Contacts For Specific Age Brackets" = "allForAge",
              "Restrict Contacts Between Two Age Brackets" = "ageForAge",
              "Decrease Probability of infection regardless of contact frequency" = "coef"
            ), 
            selected = if (hasInfo) interventionInfo$intrv[[id]]$inttype else "allContacts"
          ),
          numericInput(
            "intstart", 
            "Intervention starts on day", 
            value = if (hasInfo) interventionInfo$intrv[[id]]$intstart else 0,
            min = 0, max = round(0.9 * input$Tmax)
          ),
          numericInput(
            "intend", 
            "Intervention ends on day", 
            value = if (hasInfo) interventionInfo$intrv[[id]]$intend else 0,
            min = 0, max = round(0.9 * input$Tmax)
          )
        ),
        column(
          6,
          conditionalPanel(
            condition="input.inttype == 'coef'",
            selectInput(
              "intWhichCoef",
              "Select Coefficient to intervene on",
              choices = c("be", "b0", "b1", "b2", "mu"),
              selected = interventionInfo$intrv[[id]]$intWhichCoef
            ),
            sliderInput(
              "intCoefMagnitude",
              "Decrease coefficient by %", 
              value = interventionInfo$intrv[[id]]$intCoefMagnitude,
              min = 0, max = 1, step = 0.01
            )
          ),
          conditionalPanel(
            condition="input.inttype == 'allContacts'",
            selectInput(
              "intAgeAllTarget",
              "Decrease all social contacts for group:",
              choices = c(
                "Exposed" = "We",
                "Asymptomatic Infections" = "W0",
                "Mild Infections" = "W1",
                "Severe Infections" = "W2",
                "Critical Infections" = "W3"
              ),
              selected = interventionInfo$intrv[[id]]$intAgeAllTarget
            ),
            sliderInput(
              "intAgeAllMagnitude",
              "Decrease all social contacts by %",  
              value = interventionInfo$intrv[[id]]$intAgeAllMagnitude,
              min = 0, max = 1, step = 0.01
            )
          ),
          conditionalPanel(
            condition="input.inttype == 'ageForAge'",
            selectInput(
              "intAgeForAgeTarget",
              "Decrease social contacts for group:",
              choices = c(
                "Exposed" = "We",
                "Asymptomatic Infections" = "W0",
                "Mild Infections" = "W1",
                "Severe Infections" = "W2",
                "Critical Infections" = "W3"
              ),
              selected = interventionInfo$intrv[[id]]$intAgeForAgeTarget
            ),
            selectInput(
              "intWhichAge1",
              "Decrease contacts between age group:", 
              choices = ageBounds,
              selected = interventionInfo$intrv[[id]]$intWhichAge1
            ),
            selectInput(
              "intWhichAge2",
              ".. and age group:",
              choices = ageBounds,
              selected = interventionInfo$intrv[[id]]$intWhichAge2
            ),
            sliderInput(
              "intAge12Magnitude",
              "Decrease social contacts between the two age groups by %",  
              value = interventionInfo$intrv[[id]]$intAge12Magnitude,
              min = 0, max = 1, step = 0.01
            )
          ),
          conditionalPanel(
            condition="input.inttype == 'allForAge'",
            selectInput(
              "intAllForAgeTarget",
              "Decrease social contacts for group:",
              choices = c(
                "Exposed" = "We",
                "Asymptomatic Infections" = "W0",
                "Mild Infections" = "W1",
                "Severe Infections" = "W2",
                "Critical Infections" = "W3"
              ),
              selected = interventionInfo$intrv[[id]]$intAllForAgeTarget
            ),
            selectInput(
              "intWhichAge0",
              "Select age group for which social contacts will be decreased:",
              choices = ageBounds,
              selected = interventionInfo$intrv[[id]]$intWhichAge0
            ),
            sliderInput(
              "intAge0Magnitude",
              "Decrease all social contacts for this age group by %", 
              value = interventionInfo$intrv[[id]]$intAge0Magnitude,
              min = 0, max = 1, step = 0.01
            )
          )
        ),
        fluidRow(column(1, actionButton("intrvSave", "Save"), offset = 11))
      )
    }
  })
  
  output$interventionInputs <- renderUI({interventionInputs()})
  
  # commit inputs to block list object
  observeEvent(input$intrvSave, {
    id <- selectedIntrv$which
    
    # set intervetnion
    if (input$inttype == "coef") {
      interventionInfo$intrv[[id]]$intervention <- setIntervention(
        from = input$intstart, 
        to = input$intend, 
        genFunc = intervenePercentGenerator, 
        target = input$intWhichCoef,
        change = -input$intCoefMagnitude
      )
    } else if (input$inttype == "ageForAge") {
      interventionInfo$intrv[[id]]$intervention <- setIntervention(
        from = input$intstart, 
        to = input$intend, 
        genFunc = interveneContactPercentGenerator, 
        ageBreaks = ageBounds, 
        targetAges = c(input$intWhichAge1, input$intWhichAge2), 
        target = input$intAgeForAgeTarget, 
        change = -input$intAge12Magnitude
      )
    } else if (input$inttype == "allForAge") {
      interventionInfo$intrv[[id]]$intervention <- setIntervention(
        from = input$intstart, 
        to = input$intend, 
        genFunc = interveneContactPercentGenerator, 
        ageBreaks = ageBounds, 
        targetAges = c(input$intWhichAge0), 
        target = input$intAllForAgeTarget, 
        change = -input$intAge0Magnitude
      )
    } else if (input$inttype == "allContacts") {
      interventionInfo$intrv[[id]]$intervention <- setIntervention(
        from = input$intstart, 
        to = input$intend, 
        genFunc = interveneContactPercentGenerator, 
        ageBreaks = ageBounds, 
        targetAges = c(), 
        target = input$intAgeAllTarget, 
        change = -input$intAgeAllMagnitude
      )
    } else {
      stop("unknown type")
    }
    
    # save input values
    interventionInfo$intrv[[id]]$inttype <- input$inttype
    interventionInfo$intrv[[id]]$intstart <- input$intstart
    interventionInfo$intrv[[id]]$intend <- input$intend
    interventionInfo$intrv[[id]]$intWhichCoef <- input$intWhichCoef
    interventionInfo$intrv[[id]]$intCoefMagnitude <- input$intCoefMagnitude
    interventionInfo$intrv[[id]]$intWhichAge1 <- input$intWhichAge1
    interventionInfo$intrv[[id]]$intWhichAge2 <- input$intWhichAge2
    interventionInfo$intrv[[id]]$intAge12Magnitude <- input$intAge12Magnitude
    interventionInfo$intrv[[id]]$intWhichAge0 <- input$intWhichAge0
    interventionInfo$intrv[[id]]$intAge0Magnitude <- input$intAge0Magnitude
    interventionInfo$intrv[[id]]$intAgeAllMagnitude <- input$intAgeAllMagnitude
    interventionInfo$intrv[[id]]$intAgeForAgeTarget <- input$intAgeForAgeTarget
    interventionInfo$intrv[[id]]$intAllForAgeTarget <- input$intAllForAgeTarget
    interventionInfo$intrv[[id]]$intAgeAllTarget <- input$intAgeAllTarget
  })
  
  # there's got to be a more elegant way..
  observeEvent(input$remIntrv1, {interventionInfo$intrv[[1]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv2, {interventionInfo$intrv[[2]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv3, {interventionInfo$intrv[[3]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv4, {interventionInfo$intrv[[4]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv5, {interventionInfo$intrv[[5]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv6, {interventionInfo$intrv[[6]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv7, {interventionInfo$intrv[[7]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv8, {interventionInfo$intrv[[8]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv9, {interventionInfo$intrv[[9]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv10, {interventionInfo$intrv[[10]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv11, {interventionInfo$intrv[[11]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv12, {interventionInfo$intrv[[12]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv13, {interventionInfo$intrv[[13]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv14, {interventionInfo$intrv[[14]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv15, {interventionInfo$intrv[[15]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv16, {interventionInfo$intrv[[16]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv17, {interventionInfo$intrv[[17]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv18, {interventionInfo$intrv[[18]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv19, {interventionInfo$intrv[[19]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv20, {interventionInfo$intrv[[20]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv21, {interventionInfo$intrv[[21]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv22, {interventionInfo$intrv[[22]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv23, {interventionInfo$intrv[[23]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  observeEvent(input$remIntrv24, {interventionInfo$intrv[[24]] <- NULL; countIntrv$n <- countIntrv$n - 1})
  
  observeEvent(input$selIntrv1, {selectedIntrv$which <- 1})
  observeEvent(input$selIntrv2, {selectedIntrv$which <- 2})
  observeEvent(input$selIntrv3, {selectedIntrv$which <- 3})
  observeEvent(input$selIntrv4, {selectedIntrv$which <- 4})
  observeEvent(input$selIntrv5, {selectedIntrv$which <- 5})
  observeEvent(input$selIntrv6, {selectedIntrv$which <- 6})
  observeEvent(input$selIntrv7, {selectedIntrv$which <- 7})
  observeEvent(input$selIntrv8, {selectedIntrv$which <- 8})
  observeEvent(input$selIntrv9, {selectedIntrv$which <- 9})
  observeEvent(input$selIntrv10, {selectedIntrv$which <- 10})
  observeEvent(input$selIntrv11, {selectedIntrv$which <- 11})
  observeEvent(input$selIntrv12, {selectedIntrv$which <- 12})
  observeEvent(input$selIntrv13, {selectedIntrv$which <- 13})
  observeEvent(input$selIntrv14, {selectedIntrv$which <- 14})
  observeEvent(input$selIntrv15, {selectedIntrv$which <- 15})
  observeEvent(input$selIntrv16, {selectedIntrv$which <- 16})
  observeEvent(input$selIntrv17, {selectedIntrv$which <- 17})
  observeEvent(input$selIntrv18, {selectedIntrv$which <- 18})
  observeEvent(input$selIntrv19, {selectedIntrv$which <- 19})
  observeEvent(input$selIntrv20, {selectedIntrv$which <- 20})
  observeEvent(input$selIntrv21, {selectedIntrv$which <- 21})
  observeEvent(input$selIntrv22, {selectedIntrv$which <- 22})
  observeEvent(input$selIntrv23, {selectedIntrv$which <- 23})
  observeEvent(input$selIntrv24, {selectedIntrv$which <- 24})
  
}
