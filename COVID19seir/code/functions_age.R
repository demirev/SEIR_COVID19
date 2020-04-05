getMossongData <- function(
  ensureSym = T,
  ageBounds = c(
    4, 7, 12, 18, 25, 30, 40, 50, 65, 75, 85, 100
  )
) {
  #' Derive default per-age contact rates using the Mossong data set
  
  findBound <- function(x, bounds) {
    if (x > max(bounds)) return(max(bounds))
    bounds[which(x <= bounds)][1]
  }
  
  result <- socialmixr::polymod$participants %>% 
    as_tibble() %>%
    select(part_id, part_age) %>%
    left_join(socialmixr::polymod$contacts, by = "part_id") %>%
    mutate(cnt_age_exact = ifelse(is.na(cnt_age_exact), round((cnt_age_est_min + cnt_age_est_max) / 2), cnt_age_exact)) %>%
    group_by(part_id, part_age, cnt_age_exact) %>%
    count() %>%
    group_by(part_age, cnt_age_exact) %>%
    summarise(n = mean(n)) %>%
    filter(!is.na(part_age) & !is.na(cnt_age_exact))
  
  result <- expand.grid(
    unique(result$part_age), unique(result$cnt_age_exact)
  ) %>%
    as_tibble() %>%
    dplyr::rename(part_age = Var1, cnt_age_exact = Var2) %>%
    left_join(result) %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  if (ensureSym) {
    result <- result %>%
      dplyr::rename(
        cnt_age_exact = part_age, part_age = cnt_age_exact, n_rev = n
      ) %>%
      inner_join(result) %>%
      mutate(n = (n + n_rev) / 2)
  }
  
  result <- result  %>%
    select(
      age_1 = cnt_age_exact, age_2 = part_age, daily_contacts = n
    )
  
  if (!is.null(ageBounds)) {
    result <- result %>%
      ungroup() %>%
      mutate(
        age_1 = map_dbl(age_1, findBound, ageBounds), 
        age_2 = map_dbl(age_2, findBound, ageBounds)
      ) %>%
      group_by(age_1, age_2) %>%
      summarise(daily_contacts = mean(daily_contacts))
  }
  
  ungroup(result)
}

SetODEs_SEIRage = function(t, y, p) {
  #' main function for setting the ODEs of the model
  nBrackets = length(p$a)
  
  # set state
  S = y[(0*nBrackets + 1) : (1*nBrackets)]
  E0 = y[(1*nBrackets + 1) : (2*nBrackets)]
  E1 = y[(2*nBrackets + 1) : (3*nBrackets)]
  I0 = y[(3*nBrackets + 1) : (4*nBrackets)]
  I1 = y[(4*nBrackets + 1) : (5*nBrackets)]
  I2 = y[(5*nBrackets + 1) : (6*nBrackets)]
  I3 = y[(6*nBrackets + 1) : (7*nBrackets)]
  R = y[(7*nBrackets + 1) : (8*nBrackets)]
  D = y[(8*nBrackets + 1) : (9*nBrackets)]
  
  # pull out parameters
  #N = p$N
  seasAmp = p$seasAmp
  seasPhase = p$seasPhase
  We = p$We
  W0 = p$W0
  W1 = p$W1
  W2 = p$W2
  W3 = p$W3
  be = p$be
  b0 = p$b0
  b1 = p$b1
  b2 = p$b2
  b3 = p$b3
  a = p$a
  eta = p$eta
  pIm = p$pIm
  phi = p$phi
  eps0 = p$eps0
  eps1 = p$eps1
  f = p$f
  theta1 = p$theta1
  theta2 = p$theta2
  gamma0 = p$gamma0
  gamma1 = p$gamma1
  gamma2 = p$gamma2
  gamma3 = p$gamma3
  mu = p$mu
  
  # define time derivatives
  seas=(1 + seasAmp*cos(2*pi*(t-seasPhase)/365))
  N = S + E0 + E1 + I0 + I1 + I2 + I3 + R
  
  dSdt = - be * We %*% E1 * S / N * seas  -
    b0 * W0 %*% I0 * S / N * seas  -
    b1 * W1 %*% I1 * S / N * seas  -
    b2 * W2 %*% I2 * S / N * seas  -
    b3 * W3 %*% I3 * S / N * seas  +
    c(eta * sum(N), rep(0, nBrackets - 1))  -
    eta * S +
    c(0, a[1 : (nBrackets - 1)]) * c(0, S[1 : (nBrackets - 1)]) * (1 - pIm) - 
    a * S +
    phi * R
  
  dE0dt = be * We %*% E1 * S / N * seas  +
    b0 * W0 %*% I0 * S / N * seas  +
    b1 * W1 %*% I1 * S / N * seas  +
    b2 * W2 %*% I2 * S / N * seas  +
    b3 * W3 %*% I3 * S / N * seas -
    eps0 * E0 -
    eta * E0 +
    c(0, a[1 : (nBrackets - 1)]) * c(0, E0[1 : (nBrackets - 1)]) - 
    a * E0
  
  dE1dt = eps0 * E0 -
    eps1 * E1 -
    eta * E1 +
    c(0, a[1 : (nBrackets - 1)]) * c(0, E1[1 : (nBrackets - 1)]) - 
    a * E1
  
  dI0dt = f * eps1 * E1 - 
    gamma0 * I0 -
    eta * I0 +
    c(0, a[1 : (nBrackets - 1)]) * c(0, I0[1 : (nBrackets - 1)]) - 
    a * I0
  
  dI1dt = (1 - f) * eps1 * E1 -
    gamma1 * I1 -
    theta1 * I1 -
    eta * I1 +
    c(0, a[1 : (nBrackets - 1)]) * c(0, I1[1 : (nBrackets - 1)]) - 
    a * I1
  
  dI2dt = theta1 * I1 -
    gamma2 * I2 -
    theta2 * I2 -
    eta * I2 +
    c(0, a[1 : (nBrackets - 1)]) * c(0, I2[1 : (nBrackets - 1)]) - 
    a * I2
  
  dI3dt = theta2 * I2 -
    gamma3 * I3 -
    eta * I3 +
    c(0, a[1 : (nBrackets - 1)]) * c(0, I3[1 : (nBrackets - 1)]) - 
    a * I3 -
    mu * I3
  
  dRdt = gamma0 * I0 +
    gamma1 * I1 +
    gamma2 * I2 +
    gamma3 * I3 -
    phi * R +
    c(0, a[1 : (nBrackets - 1)]) * c(0, S[1 : (nBrackets - 1)]) * pIm -
    eta * R +
    c(0, a[1 : (nBrackets - 1)]) * c(0, R[1 : (nBrackets - 1)]) - 
    a * R 
  
  dDdt = mu * I3
  
  return(
    list(c(dSdt, dE0dt, dE1dt, dI0dt, dI1dt, dI2dt, dI3dt, dRdt, dDdt))
  )
}


GetModelParamsAge = function(input) {
  #' converty shiny inputs to parameters
  
  IncubPeriod=input$IncubPeriod  #Incubation period, days
  DurMildInf=input$DurMildInf #Duration of mild infections, days
  FracSevere=input$FracSevere #Fraction of infections that are severe
  FracCritical=input$FracCritical #Fraction of infections that are critical
  FracMild=1-FracSevere-FracCritical  #Fraction of infections that are mild
  ProbDeath=input$ProbDeath  #Probability of dying given critical infection
  CFR=ProbDeath*FracCritical #Case fatality rate (fraction of infections resulting in death)
  TimeICUDeath=input$TimeICUDeath #Time from ICU admission to death, days
  DurHosp=input$DurHosp #Duration of hospitalization, days
  ageBounds = input$ageBounds
  ExposedContacts = input$ExposedContacts
  AssymContacts = input$AssymContacts
  MildContacts = input$MildContacts
  SevereContacts = input$SevereContacts
  CritContacts = input$CritContacts
  NaturalDeathRate = input$NaturalDeathRate 
  InoculationSchedule = input$InoculationSchedule
  ReinfectionChance = input$ReinfectionChance
  
  N=input$N
  
  # If seasonality is allowed. If there is seasonality, the input beta values correspond to the current values. Must be adjusted to find the true (average) beta values
  
  if(input$AllowSeason=="Yes"){
    seas.amp=input$seas.amp/100 #relative amplitude of seasonal fluctuations, in [0,1]
    seas.phase=input$seas.phase #phase of seasonal fluctuations, measuered in days relative to time zero when peak will occur (0=peak occurs at time zero, 30 = peak occurs one month after time zero). Can be negative
  }else{
    seas.amp=0.0 
    seas.phase=0
  }
  seas0=(1 + seas.amp*cos(2*pi*seas.phase/365)) #value of seasonality coefficient at time zero
  
  # The transmission rates are changed from values per time to values per capita per time
  b1=input$b1/(seas0)
  b2=input$b2/(seas0)
  b3=input$b3/(seas0)
  
  #If asymptomatic infection is allowed
  if(input$AllowAsym=="Yes"){
    FracAsym=input$FracAsym #Fraction of all infections that are asymptomatic
    DurAsym=input$DurAsym #Duration of asympatomatic infection
    b0=input$b0/(seas0)
  }else{
    FracAsym=0 #Fraction of all infections that are asymptomatic
    DurAsym=7 #Duration of asympatomatic infection
    b0 = 0 #Transmission rate (asymptomatic infections)
  }
  
  # If presymptomatic transmission is allowed
  if(input$AllowPresym=="Yes"){
    PresymPeriod=input$PresymPeriod #Length of infections phase of incubation period
    be=input$be/(seas0)
  }else{
    PresymPeriod=0 #Length of infectious phase of incubation period
    be = 0 #Transmission rate (pre-symptomatic)
  }
  
  # Turn these clinical parameters into the rate constants of the model
  a1=min(10^6,1/PresymPeriod) #presymptomatic period of transmission
  a0=min(10^6,(IncubPeriod-PresymPeriod)^(-1)) # true latent period, avoid infinity when no presymptomatic phase
  
  f=FracAsym
  
  g0=1/DurAsym
  
  g1=(1/DurMildInf)*FracMild
  p1=(1/DurMildInf)-g1
  
  p2=(1/DurHosp)*(FracCritical/(FracSevere+FracCritical))
  g2=(1/DurHosp)-p2
  
  u = ifelse(FracCritical == 0, 0, (1/TimeICUDeath)*(CFR/FracCritical))
  
  g3=(1/TimeICUDeath)-u
  
  # convert to daily
  ageBounds <- ageBounds * 365
  dailyAgeRate = c(1/diff(ageBounds), 0)
  
  pModelAge = list(
    N = N,
    seasAmp = seas.amp,
    seasPhase = seas.phase,
    We = ExposedContacts,
    W0 = AssymContacts,
    W1 = MildContacts,
    W2 = SevereContacts,
    W3 = CritContacts,
    be = be,
    b0 = b0,
    b1 = b1,
    b2 = b2,
    b3 = b3,
    a = dailyAgeRate,
    eta = NaturalDeathRate,
    pIm = InoculationSchedule,
    phi = ReinfectionChance,
    eps0 = a0,
    eps1 = a1,
    f = f,
    theta1 = p1,
    theta2 = p2,
    gamma0 = g0,
    gamma1 = g1,
    gamma2 = g2,
    gamma3 = g3,
    mu = u
  )
}

simOnceAge <- function(
  pModelAge,
  InitInf,
  ageBounds,
  Tmax
) {
  #' solve the system odf ODEs oncew
  N=pModelAge$N
  
  # Set initial conditions and time interval
  E00=InitInf
  S0 = N-E00
  y0 = c(
    S=S0, 
    E0 = if(length(E00) == 1) rep(E00, length(ageBounds)) else E00,  
    E1=rep(0, length(ageBounds)), 
    I0=rep(0, length(ageBounds)), 
    I1=rep(0, length(ageBounds)), 
    I2=rep(0, length(ageBounds)), 
    I3=rep(0, length(ageBounds)), 
    R=rep(0, length(ageBounds)),
    D=rep(0, length(ageBounds))
  )
  
  #run ODEs
  outDf=GetSpread_SEIR(pModelAge,Tmax,y0, setterFunc = SetODEs_SEIRage)
}

plotSimAggregate <- function(
  odeRes, usePlotly = T
) {
  #' plot aggregate values of each SEIRD group over the timeframe of simulation
  odeRes <- odeRes  %>%
    gather("key", "value", -time) %>%
    mutate(N = sum(odeRes[1,-1]), value_perc = value / N * 100) %>%
    mutate(
      group = NA,
      group = ifelse(stringr::str_detect(key, "S"), "S", group),
      group = ifelse(stringr::str_detect(key, "E0"), "E0", group),
      group = ifelse(stringr::str_detect(key, "E1"), "E1", group),
      group = ifelse(stringr::str_detect(key, "I0"), "I0", group),
      group = ifelse(stringr::str_detect(key, "I1"), "I1", group),
      group = ifelse(stringr::str_detect(key, "I2"), "I2", group),
      group = ifelse(stringr::str_detect(key, "I3"), "I3", group),
      group = ifelse(stringr::str_detect(key, "R"), "R", group),
      group = ifelse(stringr::str_detect(key, "D"), "D", group),
      age = stringr::str_replace(key, group, "")
    ) %>%
    mutate(
      group = ifelse(group %in% c("E0", "E1"), "E", group),
      group = ifelse(group %in% c("I0", "I1", "I2", "I3"), "I", group)
    ) %>%
    group_by(time, group) %>%
    summarize(value = sum(value), value_perc = sum(value_perc)) %>%
    ungroup() %>%
    mutate(
      ord_ = NA,
      ord_ = ifelse(group == "S", 1, ord_),
      ord_ = ifelse(group == "E", 2, ord_),
      ord_ = ifelse(group == "I", 3, ord_),
      ord_ = ifelse(group == "R", 4, ord_),
      ord_ = ifelse(group == "D", 5, ord_),
      group = ifelse(group == "S", "Susceptible", group),
      group = ifelse(group == "E", "Exposed", group),
      group = ifelse(group == "I", "Infected", group),
      group = ifelse(group == "R", "Recovered", group),
      group = ifelse(group == "D", "Deceased", group)
    ) %>%
    arrange(time, ord_) %>%
    mutate(
      group = factor(
        group, levels = c(
          "Susceptible", 
          "Exposed", 
          "Infected",
          "Recovered",
          "Deceased"
        )
      )
    )
  
  if (usePlotly) {
    odeRes %>%
      mutate(
        caption = paste0(
          round(value), " people (",
          round(value_perc, 2),
          "%) are ", group, " at day ", time
        )
      ) %>%
      plot_ly(
        x=~time, 
        y=~value, 
        color=~group, 
        type='scatter', 
        mode='lines',
        text = ~caption,
        hoverinfo = 'text'
      ) %>%
      layout(
        xaxis=list(title="Time since introduction (days)"),
        yaxis=list(title="Number of people")
      )
  } else {
    odeRes %>% 
      ggplot(aes(x = time, y = value, color = group)) + 
      geom_line() +
      theme_bw() +
      xlab("Time since introduction (days)") +
      ylab("Number of people")
  }
  
}

getMaxValEandI <- function(
  odeRes
) {
  #' get the maximum value of the y-axis for the EandI plot
  #' TODO: there is code duplication with this function and the plot function
  odeRes %>%
    gather("key", "value", -time) %>%
    mutate(N = sum(odeRes[1,-1]), value_perc = value / N * 100) %>%
    mutate(
      group = NA,
      group = ifelse(stringr::str_detect(key, "S"), "S", group),
      group = ifelse(stringr::str_detect(key, "E0"), "Exposed but cannot spread", group),
      group = ifelse(stringr::str_detect(key, "E1"), "Exposed and can spread", group),
      group = ifelse(stringr::str_detect(key, "I0"), "Asymptomatic Infections", group),
      group = ifelse(stringr::str_detect(key, "I1"), "Mild Infections", group),
      group = ifelse(stringr::str_detect(key, "I2"), "Severe Infections (Hospital stay)", group),
      group = ifelse(stringr::str_detect(key, "I3"), "Critical Infections (ICU)", group),
      group = ifelse(stringr::str_detect(key, "R"), "R", group),
      group = ifelse(stringr::str_detect(key, "D"), "D", group),
      age = stringr::str_replace(key, group, "")
    ) %>%
    filter(
      group %in% c("Exposed but cannot spread", "Exposed and can spread") | 
        group %in% c(
          "Asymptomatic Infections", 
          "Mild Infections",
          "Severe Infections (Hospital stay)",
          "Critical Infections (ICU)"
        )
    ) %>%
    group_by(time, group) %>%
    summarize(value = sum(value), value_perc = sum(value_perc)) %>%
    ungroup() %>%
    mutate(
      ord_ = NA,
      ord_ = ifelse(group == "Exposed but cannot spread", 1, ord_),
      ord_ = ifelse(group == "Exposed and can spread", 2, ord_),
      ord_ = ifelse(group == "Asymptomatic Infections", 3, ord_),
      ord_ = ifelse(group == "Mild Infections", 4, ord_),
      ord_ = ifelse(group == "Severe Infections (Hospital stay)", 5, ord_),
      ord_ = ifelse(group == "Critical Infections (ICU)", 5, ord_)
    ) %>%
    arrange(time, ord_) %>%
    mutate(
      group = factor(
        group, levels = c(
          "Exposed but cannot spread", 
          "Exposed and can spread", 
          "Asymptomatic Infections",
          "Mild Infections",
          "Severe Infections (Hospital stay)",
          "Critical Infections (ICU)"
        )
      )
    ) %>%
    pull(value) %>%
    max()
}

plotSimEandI <- function(
  odeRes, usePlotly = T, maxY = NULL
) {
  #' detailed plot of the I and E groups
  odeRes <- odeRes %>%
    gather("key", "value", -time) %>%
    mutate(N = sum(odeRes[1,-1]), value_perc = value / N * 100) %>%
    mutate(
      group = NA,
      group = ifelse(stringr::str_detect(key, "S"), "S", group),
      group = ifelse(stringr::str_detect(key, "E0"), "Exposed but cannot spread", group),
      group = ifelse(stringr::str_detect(key, "E1"), "Exposed and can spread", group),
      group = ifelse(stringr::str_detect(key, "I0"), "Asymptomatic Infections", group),
      group = ifelse(stringr::str_detect(key, "I1"), "Mild Infections", group),
      group = ifelse(stringr::str_detect(key, "I2"), "Severe Infections (Hospital stay)", group),
      group = ifelse(stringr::str_detect(key, "I3"), "Critical Infections (ICU)", group),
      group = ifelse(stringr::str_detect(key, "R"), "R", group),
      group = ifelse(stringr::str_detect(key, "D"), "D", group),
      age = stringr::str_replace(key, group, "")
    ) %>%
    filter(
      group %in% c("Exposed but cannot spread", "Exposed and can spread") | 
        group %in% c(
          "Asymptomatic Infections", 
          "Mild Infections",
          "Severe Infections (Hospital stay)",
          "Critical Infections (ICU)"
        )
    ) %>%
    group_by(time, group) %>%
    summarize(value = sum(value), value_perc = sum(value_perc)) %>%
    ungroup() %>%
    mutate(
      ord_ = NA,
      ord_ = ifelse(group == "Exposed but cannot spread", 1, ord_),
      ord_ = ifelse(group == "Exposed and can spread", 2, ord_),
      ord_ = ifelse(group == "Asymptomatic Infections", 3, ord_),
      ord_ = ifelse(group == "Mild Infections", 4, ord_),
      ord_ = ifelse(group == "Severe Infections (Hospital stay)", 5, ord_),
      ord_ = ifelse(group == "Critical Infections (ICU)", 5, ord_)
    ) %>%
    arrange(time, ord_) %>%
    mutate(
      group = factor(
        group, levels = c(
          "Exposed but cannot spread", 
          "Exposed and can spread", 
          "Asymptomatic Infections",
          "Mild Infections",
          "Severe Infections (Hospital stay)",
          "Critical Infections (ICU)"
        )
      )
    )
  
  if (is.null(maxY)) maxY <- 1.1 * max(odeRes$value)
  
  if (usePlotly) {
    odeRes %>%
      mutate(
        caption = paste0(
          round(value), " people (",
          round(value_perc, 2),
          "%) are ", group, " at day ", time
        )
      ) %>%
      plot_ly(
        x=~time, 
        y=~value, 
        color=~group, 
        type='scatter',
        mode='lines',
        text = ~caption,
        hoverinfo = 'text'
      ) %>%
      layout(
        xaxis=list(title="Time since introduction (days)"),
        yaxis=list(range = c(0, maxY), title="Number of people")
      )
  } else {
    odeRes %>% 
      ggplot(aes(x = time, y = value, color = group)) + 
      geom_line() +
      theme_bw() +
      xlab("Time since introduction (days)") +
      ylab("Number of people")
  }
  
}

plotSimByAge <- function(
  odeRes, usePlotly = T, whichGroup = "S",
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
) {
  #' by age plot
  odeRes <- odeRes %>%
    gather("key", "value", -time) %>%
    mutate(
      group = NA,
      group = ifelse(stringr::str_detect(key, "S"), "S", group),
      group = ifelse(stringr::str_detect(key, "E0"), "E0", group),
      group = ifelse(stringr::str_detect(key, "E1"), "E1", group),
      group = ifelse(stringr::str_detect(key, "I0"), "I0", group),
      group = ifelse(stringr::str_detect(key, "I1"), "I1", group),
      group = ifelse(stringr::str_detect(key, "I2"), "I2", group),
      group = ifelse(stringr::str_detect(key, "I3"), "I3", group),
      group = ifelse(stringr::str_detect(key, "R"), "R", group),
      group = ifelse(stringr::str_detect(key, "D"), "D", group),
      age = stringr::str_replace(key, group, "")
    ) %>%
    filter(
      group == whichGroup
    ) %>%
    mutate(
      ord_ = NA,
      ord_ = ifelse(group == "E0", 1, ord_),
      ord_ = ifelse(group == "E1", 2, ord_),
      ord_ = ifelse(group == "I0", 3, ord_),
      ord_ = ifelse(group == "I1", 4, ord_),
      ord_ = ifelse(group == "I2", 5, ord_),
      ord_ = ifelse(group == "I3", 5, ord_)
    ) %>%
    arrange(time, ord_) %>%
    mutate(
      group = factor(
        group, levels = c(
          "E0", 
          "E1", 
          "I0",
          "I1",
          "I2",
          "I3"
        )
      )
    ) %>%
    mutate(
      age = map_chr(age, function(ag) names(ageBounds)[as.numeric(ag)])
    )
  
  if (usePlotly) {
    odeRes %>%
      plot_ly(
        x=~time, y=~value, fill=~age, type='scatter', mode='lines'
      ) %>%
      layout(
        xaxis=list(title="Time since introduction (days)"),
        yaxis=list(title="Number of people")
      )
  } else {
    odeRes %>% 
      mutate(age = factor(age, levels = names(ageBounds))) %>%
      ggplot(aes(x = time, y = value, fill = age, color = age)) + 
      #scale_fill_brewer(palette = "Set1") +
      scale_fill_hue(c=45, l=80) +
      geom_area() +
      theme_bw() +
      xlab("Time since introduction (days)") +
      ylab("Number of people")
  }
  
}

simInterventionsAge <- function(
  interventions = tibble(),
  IncubPeriod= 5, #Duration of incubation period
  DurMildInf= 6, #Duration of mild infections
  FracSevere= 15, #% of symptomatic infections that are severe
  FracCritical= 5, #% of symptomatic infections that are critical
  ProbDeath= 40, #Death rate for critical infections
  DurHosp= 6, #Duration of severe infection/hospitalization
  TimeICUDeath= 8, #Duration critical infection/ICU stay
  
  AllowPresym="No",
  AllowAsym="No",
  FracAsym=25, #Fraction of all infections that are asymptomatic
  PresymPeriod=2, #Length of infectious phase of incubation period
  DurAsym=6, #Duration of asympatomatic infection
  
  b1= 0.5, #Transmission rate (mild infections)
  b2= 0.01, #Transmission rate (severe infections)
  b3= 0.01, #Transmission rate (critical infections)
  be = 0.5, #Transmission rate (pre-symptomatic)
  b0 = 0.5, #Transmission rate (asymptomatic infections)
  Tmax= 300, #Maximum time
  InitInf= 1, #Initial # infected
  
  AllowSeason="No",
  seas.amp=0.0, #relative amplitude of seasonal fluctuations, in [0,1]
  seas.phase=-90,
  
  NaturalDeathRate = 0,
  InoculationSchedule = 0,
  ReinfectionChance = 0,
  
  ageBounds = c(
    4, 7, 12, 15, 18, 22, 25, 30, 35, 40, 45, 
    50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100
  ),
  N = rep(1000, length(ageBounds)), #Total population size
  ExposedContacts = matrix(
    0.03, nrow = length(ageBounds), ncol = length(ageBounds)
  ),
  AssymContacts = matrix(
    0.03, nrow = length(ageBounds), ncol = length(ageBounds)
  ),
  MildContacts = matrix(
    0.03, nrow = length(ageBounds), ncol = length(ageBounds)
  ),
  SevereContacts = matrix(
    0.03, nrow = length(ageBounds), ncol = length(ageBounds)
  ),
  CritContacts = matrix(
    0.03, nrow = length(ageBounds), ncol = length(ageBounds)
  )
) {
  #' solve ODEs given inputs and a set of interventions described in tibble format
  #' @param interventions a tibble of three columns: start and end that define
  #' when the intervention starts and ends and func. The column func is a list of functions,
  #' where each entry is a function that takes in the list of parameters, does some alteration
  #' and returns the same format list
  pModelAge=GetModelParamsAge(
    list(
      "IncubPeriod" = IncubPeriod,
      "DurMildInf" = DurMildInf,
      "FracSevere" = FracSevere,
      "FracCritical" = FracCritical,
      "ProbDeath" = ProbDeath,
      "DurHosp" = DurHosp,
      "TimeICUDeath"=TimeICUDeath,
      "FracAsym"=FracAsym, 
      "PresymPeriod"=PresymPeriod, 
      "DurAsym"=DurAsym, 
      "be"=be,
      "b0"=b0, 
      "b1"=b1,
      "b2"=b2,
      "b3"=b3,
      "seas.amp"=seas.amp, 
      "seas.phase"=seas.phase,
      "N"=N,
      "Tmax"=Tmax,
      "InitInf"=InitInf,
      #"yscale"=yscale,
      "AllowPresym"=AllowPresym,
      "AllowAsym"=AllowAsym,
      "AllowSeason"=AllowSeason,
      #"PlotCombine"=PlotCombine
      "ageBounds"=ageBounds,
      "ExposedContacts" = ExposedContacts,
      "AssymContacts" = AssymContacts,
      "MildContacts" = MildContacts,
      "SevereContacts" = SevereContacts,
      "CritContacts" = CritContacts,
      "NaturalDeathRate" = NaturalDeathRate,
      "InoculationSchedule" = InoculationSchedule,
      "ReinfectionChance" = ReinfectionChance
    )
  ) # calculate base params
  
  #pModel=ParamStruct$pModel
  N=pModelAge$N
  
  # Set initial conditions and time interval
  E00=InitInf
  S0 = N-E00
  y0 = c(
    S=S0, 
    E0 = if(length(E00) == 1) rep(E00, length(ageBounds)) else E00,  
    E1=rep(0, length(ageBounds)), 
    I0=rep(0, length(ageBounds)), 
    I1=rep(0, length(ageBounds)), 
    I2=rep(0, length(ageBounds)), 
    I3=rep(0, length(ageBounds)), 
    R=rep(0, length(ageBounds)),
    D=rep(0, length(ageBounds))
  )
  
  if (!length(interventions)) {
    interventions <- tibble(start = 0, end = 0, func = list(function(x) x))[0,]
  }
  
  if (nrow(interventions)) {
    interventions <- filter(interventions, start <= Tmax) %>%
      mutate(end = ifelse(end > Tmax, Tmax - 5, end))
  }
  
  tBreaks <- sort(unique(c(0, interventions$start, interventions$end)))
  
  yt = y0
  outDf = tibble()
  
  for (i in seq_along(tBreaks)) {
    pModelT <- pModelAge    
    
    tfrom <- tBreaks[[i]]
    tto <- if (i == length(tBreaks)) Tmax else tBreaks[[i+1]]
    
    activeInterventions <- interventions %>%
      filter(start <= tfrom & end >= tto)
    
    for (interventionFunc in activeInterventions$func) {
      pModelT <- interventionFunc(pModelT)
    }
    
    #print(list(from = tfrom, to = tto, p = pModelT))
    
    resT <- GetSpread_SEIR(
      pModelT, 
      (tto - tfrom), 
      yt, 
      setterFunc = SetODEs_SEIRage
    )
    
    outDf=bind_rows(
      outDf,
      resT[-1,]
    )
    
    yt = as.numeric(outDf[nrow(outDf), -1])
    names(yt) = colnames(outDf)[-1]
  }
  
  mutate(outDf, time = 1:nrow(outDf))
}

intervenePercentGenerator <- function(target = 'b0', change = -0.1) {
  #' generator function for creating funcs in intervention tibble
  func <- function(inputs) {
    inputs[[target]] <- inputs[[target]] + change * inputs[[target]]
    inputs
  }
}

interveneAdditionGenerator <- function(target = 'b0', change = -0.1) {
  #' generator function for creating funcs in intervention tibble
  func <- function(inputs) {
    inputs[[target]] <- inputs[[target]] + change
    inputs
  }
}

interveneConstantGenerator <- function(target = 'b0', change = -0.1) {
  #' generator function for creating funcs in intervention tibble
  func <- function(inputs) {
    inputs[[target]] <- change
    inputs
  }
}

interveneContactPercentGenerator <- function(
  ageBreaks, targetAges, target = 'Wo', change = -0.1
) {
  #' generator function for creating funcs in intervention tibble
  if (length(targetAges) == 2) {
    func <- function(inputs) {
      i1 <- ageBreaks == targetAges[1]
      i2 <- ageBreaks == targetAges[2]
      inputs[[target]][i1, i2] <- inputs[[target]][i1, i2] + change * inputs[[target]][i1, i2]
      inputs[[target]][i2, i1] <- inputs[[target]][i2, i1] + change * inputs[[target]][i2, i1]
      inputs
    }
  } else if (length(targetAges) == 1) {
    func <- function(inputs) {
      i1 <- ageBreaks == targetAges[1]
      inputs[[target]][i1, !i1] <- inputs[[target]][i1, !i1] + change * inputs[[target]][i1, !i1]
      inputs[[target]][!i1, i1] <- inputs[[target]][!i1, i1] + change * inputs[[target]][!i1, i1] 
      inputs[[target]][i1, i1] <- inputs[[target]][i1, i1] + change * inputs[[target]][i1, i1]
      inputs
    }
  } else {
    func <- function(inputs) {
      inputs[[target]] <- inputs[[target]] + change * inputs[[target]]
      inputs
    }
  }
  
  func
}

interveneContactAdditionGenerator <- function(
  ageBreaks, targetAges, target = 'Wo', change = -0.1
) {
  #' generator function for creating funcs in intervention tibble
  if (length(targetAges) == 2) {
    func <- function(inputs) {
      i1 <- ageBreaks == targetAges[1]
      i2 <- ageBreaks == targetAges[2]
      inputs[[target]][i1, i2] <- inputs[[target]][i1, i2] + change
      inputs[[target]][i2, i1] <- inputs[[target]][i2, i1] + change
      inputs
    }
  } else if (length(targetAges) == 1) {
    func <- function(inputs) {
      i1 <- ageBreaks == targetAges[1]
      inputs[[target]][i1, !i1] <- inputs[[target]][i1, !i1] + change
      inputs[[target]][!i1, i1] <- inputs[[target]][!i1, i1] + change
      inputs[[target]][i1, i1] <- inputs[[target]][i1, i1] + change 
      inputs
    }
  } else {
    func <- function(inputs) {
      inputs[[target]] <- inputs[[target]] + change 
      inputs
    }
  }
  
  func
}

interveneContactConstantGenerator <- function(
  ageBreaks, targetAges, target = 'Wo', change = .1
) {
  #' generator function for creating funcs ii intervention tibble
  if (length(targetAges) == 2) {
    func <- function(inputs) {
      i1 <- ageBreaks == targetAges[1]
      i2 <- ageBreaks == targetAges[2]
      inputs[[target]][i1, i2] <- change
      inputs[[target]][i2, i1] <- change 
      inputs
    }
  } else if (length(targetAges) == 1) {
    func <- function(inputs) {
      i1 <- ageBreaks == targetAges[1]
      inputs[[target]][i1, !i1] <- change
      inputs[[target]][!i1, i1] <- change 
      inputs[[target]][i1, i1] <- change 
      inputs
    }
  } else {
    func <- function(inputs) {
      inputs[[target]] <- change
      inputs
    }
  }
  
  func
}

setIntervention <- function(from, to, genFunc, ...) {
  #' given start and end dates and a generator function 
  #' (and the arguments to the generator function) creaate a single entry
  #' into the interventions tibble
  tibble(
    start = from,
    end = to,
    func = list(genFunc(...))
  )
}

