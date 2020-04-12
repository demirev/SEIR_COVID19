library(shiny)
library(shinyWidgets)
library(plotly)

fluidPage(
  titlePanel("Age-Specific Pandemic Model"),
  hr(),
  fluidRow(
    column(
      12,
      p(div(HTML("Disclaimer: This dashboard is based on a tool by Alison Hill which can be found at <a href=https://alhill.shinyapps.io/COVID19seir/> this address </a>. As such it is also distributed under a <a href=https://creativecommons.org/licenses/by-sa/4.0/> Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) License </a>. Changes made include the introduction of different age groups and contact matrices as well as the option for multiple interventions. All errors and omissions in the changes introduced in this version are my own."))),
      offset = 0
    )
  ),
  fluidRow(
    column(
      12,
      p(div(HTML("Disclaimer: I am not an epidemiologist and all work on the tool below has been strictly for my own educational purposes. I make no claims that this tool can be used to accurately model any real pandemic. There is considerable uncertainty in all default parameters. This tool is not intended to be taken as policy advice or expert opinion. This dashboard is a work-in-progress. It has not been independently reviewed and some bugs and errors are to be expected."))),
      offset = 0
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      navbarPage(
        "Parameters:",
        tabPanel(
          "General Parameters",
          h4(div(HTML("<em>Set simulation values...</em>"))),
          #sliderInput("LogN", div(HTML("Total population size (log10)")), 1, 9, 3, step=0.1),
          #htmlOutput("N"),
          column(
            width=12,
            numericInput("InitInf","Initial # infected per age group:",value = 3, min = 1, step = 1)
          ),
          sliderInput("Tmax", div(HTML("Maximum time")),100, 3000, 1000, step=10, post=" days"),
          br(),
          h4(div(HTML("<em>Set clinical parameters for all age groups...</em>"))),
          sliderInput("IncubPeriod", "Duration of incubation period", 0, 20, 5, step=0.5, post = " days"),
          sliderInput("DurMildInf", "Duration of mild infections", 0, 20, 6, step=1, post = " days"),
          sliderInput("DurHosp", "Duration of severe infection (hospital stay)", 0, 10, 6, step=1, post = " days"),
          sliderInput("TimeICUDeath", "Duration critical infection (ICU stay)", 0, 30, 8, step=1, post = " days"),
          sliderInput("NaturalDeath", "Daily probability of dying of other causes", 0, 1, 0, step=0.01)
        ),
        tabPanel(
          "Transmission Paramters",
        
          h4(div(HTML("<em>Set probability of transmission...</em>"))),
          sliderInput("b1", div(HTML("Probability of transmission after contact with mild infections")), 0, 1, 0.033, step=0.001, post=""), # source: https://www.medrxiv.org/content/10.1101/2020.03.24.20042606v1
          sliderInput("b2", div(HTML("Probability of transmission after contact with severe infections")),0, 1, 0.062, step=0.001, post=""), # source: https://www.medrxiv.org/content/10.1101/2020.03.24.20042606v1
          sliderInput("b3", div(HTML("Probability of transmission after contact with critical infections")),0, 1, 0.062, step=0.001, post=""), # source: https://www.medrxiv.org/content/10.1101/2020.03.24.20042606v1
          sliderInput("SocConMild", "Social contacts of mild cases relative to healthy individuals", 0, 1, 0.8, step=0.01),
          sliderInput("SocConSevere", "Social contacts of severe cases relative to healthy individuals", 0, 1, 0.3, step=0.01),
          sliderInput("SocConCritical", "Social contacts of critical cases relative to healthy individuals", 0, 1, 0.1, step=0.01),
          sliderInput("LoseImunity", "Daily probability of losing immunity for recovered cases", 0, 1, 0.0, step=0.01),
          radioButtons("AllowAsym", "Allow asymptomatic infections?",
                       choices = list("Yes" = "Yes","No" = "No"),inline=TRUE,selected="Yes"),
          conditionalPanel(
            condition="input.AllowAsym == 'Yes'",
            sliderInput("FracAsym", "% of infections that are asymptomatic", 0, 1, 0.2, step=0.01, pre=),
            sliderInput("DurAsym", "Duration of asymptomatic infections", 1, 20, 6, step=1, post = " days"),
            sliderInput("b0", div(HTML("Asymptomatic transmission probability")), 0, 1, 0.0033, step=0.001, post="")  # source: https://www.medrxiv.org/content/10.1101/2020.03.24.20042606v1
          ),
          radioButtons(
            "AllowPresym", "Allow pre-symptomatic transmission?",
            choices = list("Yes" = "Yes","No" = "No"),inline=TRUE, selected="Yes"
          ),
          conditionalPanel(
            condition="input.AllowPresym == 'Yes'",
            sliderInput("PresymPeriod", "Time before symptom onset at which transmission is possible", 0, 3, 2, step=0.5, post = " days"), #Make reactive
            sliderInput("be", div(HTML("Presymptomatic transmission probability")),0, 1, 0.0033, step=0.001, post="")
          ),
          radioButtons("AllowSeason", "Allow seasonality in transmission?",
                       choices = list("Yes" = "Yes","No" = "No"),inline=TRUE, selected="Yes"),
          conditionalPanel(
            condition="input.AllowSeason == 'Yes'",
            sliderInput("seas.amp", "Amplitude of seasonality", 0, 100, 0, step=10, pre="%"),
            sliderInput("seas.phase", "Day of peak transmission (relative to t=0)", -365, 365, 0, step=1, post = " days")
          )
          
        ),
        tabPanel(
          "Population Per Age Group",
          numericInput("N1", div(HTML("Population < 4:")), value=250105, max=10^10, min=1000, step=1000),
          numericInput("N2", div(HTML("Population 4 - 7:")), value=285689, max=10^10, min=1000, step=1000),
          numericInput("N3", div(HTML("Population 8 - 12:")), value=353589, max=10^10, min=1000, step=1000),
          numericInput("N4", div(HTML("Population 13 - 18:")), value=374042, max=10^10, min=1000, step=1000),
          numericInput("N5", div(HTML("Population 19 - 25:")), value=501729, max=10^10, min=1000, step=1000),
          numericInput("N6", div(HTML("Population 26 - 35:")), value=950932, max=10^10, min=1000, step=1000),
          numericInput("N7", div(HTML("Population 36 - 50:")), value=1556518, max=10^10, min=1000, step=1000),
          numericInput("N8", div(HTML("Population 51 - 65:")), value=1452168, max=10^10, min=1000, step=1000),
          numericInput("N9", div(HTML("Population 66 - 75:")), value=824867, max=10^10, min=1000, step=1000),
          numericInput("N10", div(HTML("Population 76 - 85:")), value=435116, max=10^10, min=1000, step=1000),
          numericInput("N11", div(HTML("Population > 85:")), value=117696, max=10^10, min=1000, step=1000)
        ), # source: https://images.populationpyramid.net/capture/?selector=%23pyramid-share-container&url=https%3A%2F%2Fwww.populationpyramid.net%2Fbulgaria%2F2017%2F%3Fshare%3Dtrue
        tabPanel(
          "Severe Infections Per Age Group",
          sliderInput("FracSevere1", "% of infections that are severe < 4", 0, 1, 0.0205, step=.01),
          sliderInput("FracSevere2", "% of infections that are severe 4 - 7", 0, 1, 0.0205, step=.01),
          sliderInput("FracSevere3", "% of infections that are severe 8 - 12", 0, 1, 0.0205, step=.01),
          sliderInput("FracSevere4", "% of infections that are severe 13 - 18", 0, 1, 0.0205, step=.01),
          sliderInput("FracSevere5", "% of infections that are severe 19 - 25", 0, 1, 0.1534, step=.01),
          sliderInput("FracSevere6", "% of infections that are severe 26 - 35", 0, 1, 0.1755, step=.01),
          sliderInput("FracSevere7", "% of infections that are severe 36 - 50", 0, 1, 0.2043, step=.01),
          sliderInput("FracSevere8", "% of infections that are severe 51 - 65", 0, 1, 0.2587, step=.01),
          sliderInput("FracSevere9", "% of infections that are severe 66 - 75", 0, 1, 0.3691, step=.01),
          sliderInput("FracSevere10", "% of infections that are severe 76 - 85", 0, 1, 0.446, step=.01),
          sliderInput("FracSevere11", "% of infections that are severe > 85", 0, 1, 0.508, step=.01)
        ), # source: https://www.cdc.gov/mmwr/volumes/69/wr/pdfs/mm6912e2-H.pdf
        tabPanel(
          "Critical Infections Per Age Group",
          sliderInput("FracCritical1", "% of infections that are critical < 4", 0, 1, 0.0001, step=.01),
          sliderInput("FracCritical2", "% of infections that are critical 4 - 7", 0, 1, 0.0001, step=.01),
          sliderInput("FracCritical3", "% of infections that are critical 8 - 12", 0, 1, 0.0001, step=.01),
          sliderInput("FracCritical4", "% of infections that are critical 13 - 18", 0, 1, 0.0001, step=.01),
          sliderInput("FracCritical5", "% of infections that are critical 19 - 25", 0, 1, 0.0266, step=.01),
          sliderInput("FracCritical6", "% of infections that are critical 26 - 35", 0, 1, 0.031, step=.01),
          sliderInput("FracCritical7", "% of infections that are critical 36 - 50", 0, 1, 0.0502, step=.01),
          sliderInput("FracCritical8", "% of infections that are critical 51 - 65", 0, 1, 0.0830, step=.01),
          sliderInput("FracCritical9", "% of infections that are critical 66 - 75", 0, 1, 0.1418, step=.01),
          sliderInput("FracCritical10", "% of infections that are critical 76 - 85", 0, 1, 0.2075, step=.01),
          sliderInput("FracCritical11", "% of infections that are critical > 85", 0, 1, 0.2075, step=.01)
        ), # source: https://www.cdc.gov/mmwr/volumes/69/wr/pdfs/mm6912e2-H.pdf
        tabPanel(
          "Death Rate of Critical Infections Per Age Group",
          h4(div(HTML("<em>Change probability of lethal outcome of critical casess...</em>"))),
          sliderInput("ProbDeath1", "Death probability for critical infections < 4", 0, 1, 0.0001, step=.01),
          sliderInput("ProbDeath2", "Death probability for critical infections 4 - 7", 0, 1, 0.0001, step=.01),
          sliderInput("ProbDeath3", "Death probability for critical infections 8 - 12", 0, 1, 0.0001, step=.01),
          sliderInput("ProbDeath4", "Death probability for critical infections 13 - 18", 0, 1, 0.0001, step=.01),
          sliderInput("ProbDeath5", "Death probability for critical infections 19 - 25", 0, 1, 0.0415, step=.01),
          sliderInput("ProbDeath6", "Death probability for critical infections 26 - 35", 0, 1, 0.0484, step=.01),
          sliderInput("ProbDeath7", "Death probability for critical infections 36 - 50", 0, 1, 0.0619, step=.01),
          sliderInput("ProbDeath8", "Death probability for critical infections 51 - 65", 0, 1, 0.2085, step=.01),
          sliderInput("ProbDeath9", "Death probability for critical infections 66 - 75", 0, 1, 0.2899, step=.01),
          sliderInput("ProbDeath10", "Death probability for critical infections 76 - 85", 0, 1, 0.3566, step=.01),
          sliderInput("ProbDeath11", "Death probability for critical infections > 85", 0, 1, 0.9084, step=.01)
        ), # source: https://www.cdc.gov/mmwr/volumes/69/wr/pdfs/mm6912e2-H.pdf
        tabPanel(
          "Contact Rate Per Age Group",
          h4(div(HTML("<em>Change default values of daily contacts between age groups...</em>"))),
          uiOutput("ageContacts")
        )
      ),
      width = 5
    ),

    
    mainPanel(
      
      navbarPage(
        "Output:",
     
        tabPanel("Spread",
          fluidPage(
            fluidRow(
              h3("Predicted cases by clinical outcome"),
              p(HTML("Simulate the natural course of an epidemic in a single population without any interventions.")),
              br(),
              column(
                width=12,
                radioButtons(
                  "plotType", "Plot Type:",
                  choices = list(
                    "All Groups" = "all",
                    "Infected and Exposed" = "IandE",
                    "By Age" = "byage"
                  ),
                  inline=TRUE,
                  selected = "IandE"
                )
              ),
              br(),
              br(),
              br(),
              conditionalPanel(
                condition="input.plotType == 'all'",
                plotlyOutput("plotBase")
              ),
              conditionalPanel(
                condition="input.plotType == 'IandE'",
                plotlyOutput("plotEandI")
              ),
              conditionalPanel(
                condition="input.plotType == 'byage'",
                selectInput(
                  "agePlotWhat", 
                  "Group to plot", 
                  choices = c("S", "E0", "E1", "I0", "I1", "I2", "I3", "R", "D")
                ),
                plotlyOutput("plotAge")
              ),
              br(),
              p(HTML("<b>User instructions:</b> The graph shows the expected numbers of individuals over time who are infected, recovered, susceptible, or dead over time. Infected individuals first pass through an exposed/incubation phase where they are asymptomatic and not infectious, and then move into a symptomatic and infections stage classified by the clinical status of infection (mild, severe, or critical). A more detailed description of the model is provided in the Model Description tab. Use the radio buttons to switch between viewing the entire population, just the exosed and infected or a breakdown of any group by age. The population size, initial condition, and parameter values used to simulate the spread of infection can be specified through the sliders located in the left-hand panel. The plot is interactive: Hover over it to get values, double-click a curve in the legend to isolate it, or single-click to remove it. Dragging over a range allows zooming."))
            )
          )
        ),
       
        tabPanel("Intervention",
          fluidPage(
            fluidRow(
              h3("Reduction in predicted infections after interventions"),
              p(HTML("Simulate the change in the time course of cases after applying a set of interventions. See the Tutorial tab for instructions.")),
              br(),
              h4("Define Interventions"),
              fluidRow(
                column(
                  3,
                  shinydashboard::box(
                    title = "",
                    width = 12, collapsible = F, collapsed = F,
                    uiOutput('interventionHeaders'),
                    fluidRow(
                      column(
                        2, style = "margin-top: 20px;", #should be relative
                        actionButton("addIntrv", "", icon = icon("plus"))     
                      ),
                      column(
                        10, 
                        textInput("newIntrv", "", placeholder = "add new intervention...")
                      )
                    )
                    #background = "olive"
                  )
                ),
                column(
                  9,
                  uiOutput('interventionInputs')       
                )
              ),
              br(),
              br(),
              column(
                width=12,
                radioButtons(
                  "plotTypeInt", "Plot Type:",
                  choices = list(
                    "All Groups" = "all",
                    "Infected and Exposed" = "IandE",
                    "By Age" = "byage"
                  ),
                  inline=TRUE,
                  selected = "IandE"
                )
              ),
              br(),
              br(),
              conditionalPanel(
                condition="input.plotTypeInt == 'all'",
                h4("no interventions"),
                plotlyOutput("plotBaseInttab"),
                h4("with interventions"),
                plotlyOutput("plotIntInttab")
              ),
              conditionalPanel(
                condition="input.plotTypeInt == 'IandE'",
                h4("no intervention"),
                plotlyOutput("plotEandIBaseInttab"),
                h4("with interventions"),
                plotlyOutput("plotEandIIntInttab")
              ),
              conditionalPanel(
                condition="input.plotTypeInt == 'byage'",
                selectInput(
                  "agePlotWhatInt", 
                  "Group to plot", 
                  choices = c("S", "E0", "E1", "I0", "I1", "I2", "I3", "R", "D")
                ),
                h4("no intervention"),
                plotlyOutput("plotAgeBaseInttab"),
                h4("with interventions"),
                plotlyOutput("plotAgeIntInttab")
              )
            ),
            fluidRow(
              br(),
              p(HTML("<b>User instructions:</b> This tab allows the user to define a series of time limited interventions in order to simulate policy impact on the spread of the disease. To create a new intervention, first enter a name and then click the plus button. Fill in parameters and click the 'save' button. Interventions can be in the form of reduction of social contacts (e.g. social distancing measures) or in the form of reduction of overall probability of infection (e.g. mask wearing measures). For a social distancing measure select the age groups affected and the population targeted. The graph above shows the ourcome without the interventions, and the graph below shows it after the interventions."))
            )
          )
        ),
         
        tabPanel(
          "Model", 
          br(),
          fluidRow(
            column(
              12,
              withMathJax(),
              h2("Model Description"),
              includeMarkdown("SEIR.Rmd"),
              br()
            )
          )
        ),
         
        tabPanel(
          "Sources",
          fluidPage(
            br(),
            includeMarkdown("Sources.Rmd")
          )
        ),
         
        tabPanel(
          "About",
          fluidPage(
          br(),
          includeMarkdown("About.Rmd")
          )
        ),
        
        tabPanel(
          "Tutorial",
          fluidPage(
            br(),
            includeMarkdown("Tutorial.Rmd")
          )
        )
         
      ),
      width=7
    )
    
  )
  
)

