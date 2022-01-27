library(SHELF)
library(shiny)
library(survival)
library(shinydashboard)
library(readxl)
library(rsconnect)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(shinyjs)

ui <- fluidPage(
  
  # Application title
  titlePanel("Delayed Treatment Effects - Weibull parameterisation"),
  
  # sidebarLayout(
  mainPanel(tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  tabsetPanel(
    # Control UI ---------------------------------
    
    tabPanel("Control", 
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 fileInput("uploadSample", "Upload your control sample"),
                 numericInput("lambda2", "lambda2", value=0.05),
                 numericInput("gamma2", "gamma2", value=1)
               ), 
               mainPanel = mainPanel(
                 plotOutput("plotControl"),
                 htmlOutput("recommendedParams")
               )
             ),
    ),
    
    # T UI ---------------------------------
    tabPanel("Eliciting T",
             fluidRow(
               column(4, 
                      textInput("limits1", label = h5("T limits"), 
                                value = "0, 6")
               ),
               column(4,
                      textInput("values1", label = h5("T values"), 
                                value = "2, 3, 4")
               ),
               column(4,
                      textInput("probs1", label = h5("Cumulative probabilities"), 
                                value = "0.25, 0.5, 0.75")
               )
             ),
             fluidRow(
               column(4, 
                      selectInput("dist1", label = h5("Distribution"), 
                                  choices =  list(Histogram = "hist",
                                                  Normal = "normal", 
                                                  'Student-t' = "t",
                                                  Gamma = "gamma",
                                                  'Log normal' = "lognormal",
                                                  'Log Student-t' = "logt",
                                                  Beta = "beta",
                                                  'Mirror gamma' = "mirrorgamma",
                                                  'Mirror log normal' = "mirrorlognormal",
                                                  'Mirror log Student-t' = "mirrorlogt",
                                                  'Best fitting' = "best"),
                                  #choiceValues = 1:8,
                                  selected = 1
                      )),
               column(4,conditionalPanel(
                 condition = "input.dist1 == 't' || input.dist1 == 'logt' || input.dist1 == 'mirrorlogt'",
                 numericInput("tdf1", label = h5("Student-t degrees of freedom"),
                              value = 3)
               )
               )
               
             ),
             plotOutput("distPlot1")
    ),
    # HR UI ---------------------------------
    tabPanel("Eliciting HR",
             fluidRow(
               column(4, 
                      textInput("limits2", label = h5("HR limits"), 
                                value = "0, 1")
               ),
               column(4,
                      textInput("values2", label = h5("HR values"), 
                                value = "0.5, 0.6, 0.7")
               ),
               column(4,
                      textInput("probs2", label = h5("Cumulative probabilities"), 
                                value = "0.25, 0.5, 0.75")
               )
             ),
             fluidRow(
               column(4, 
                      selectInput("dist2", label = h5("Distribution"), 
                                  choices =  list(Histogram = "hist",
                                                  Normal = "normal", 
                                                  'Student-t' = "t",
                                                  Gamma = "gamma",
                                                  'Log normal' = "lognormal",
                                                  'Log Student-t' = "logt",
                                                  Beta = "beta",
                                                  'Mirror gamma' = "mirrorgamma",
                                                  'Mirror log normal' = "mirrorlognormal",
                                                  'Mirror log Student-t' = "mirrorlogt",
                                                  'Best fitting' = "best"),
                                  #choiceValues = 1:8,
                                  selected = 1
                      )),
               column(4,
                      conditionalPanel(
                        condition = "input.dist2 == 't' || input.dist2 == 'logt' || input.dist1 == 'mirrorlogt'",
                        numericInput("tdf2", label = h5("degrees of freedom"),
                                     value = 3)
                        
                        
                      ))
               
             ),
             
             plotOutput("distPlot2"),
             htmlOutput("HRProportion")
    ),
    
    # Feedback UI ---------------------------------
    tabPanel("Feedback", 
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 checkboxGroupInput("showfeedback", "Add to plot", choices = c("Median survival line", "95% CI for T", "CI for Survival Curves (0.1 and 0.9)")),
                 #numericInput("triallength", "How long do you expect to run the trial for? (months)", value=18),
                 numericInput("clinicaldiff", "What is the minimum clinical difference you need to see here?", value = 10)
               ),
               mainPanel = mainPanel(
                 plotOutput("plotFeedback"),
                 htmlOutput("TrialFeedback")
               )
             ),
    ),
    
    # Assurance UI ---------------------------------
    tabPanel("Assurance",
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 useShinyjs(),
                 numericInput("minlength", "Minimum trial length?", value=24),
                 numericInput("maxlength", "Maximum trial length?", value=48),
                 box(width = 10, title = "Ratio of patients in each group?",
                     splitLayout(
                       numericInput("n1", "Control", value=1, min=1),
                       numericInput("n2", "Treatment", value=1, min=1)
                     )
                 ),
                 numericInput("chosenassurance", "What assurance do we want?", value=0.5),
                 actionButton("drawSSvMonths", "Produce plot"),
                 numericInput("chosenlength", "How long will you run the trial for?", value=45)
               ), 
               mainPanel = mainPanel(
                 plotOutput("samplevmonths"),
                 htmlOutput("samplesizerequired"),
                 plotOutput("eventsPlot")
               )
             ),
             
  ),
  wellPanel(
    fluidRow(
      column(3, selectInput("outFormat", label = "Report format",
                            choices = list('html' = "html_document",
                                           'pdf' = "pdf_document",
                                           'Word' = "word_document"))
      ),
      column(3, offset = 1,
             numericInput("fs", label = "Font size", value = 12)
      )),
    fluidRow(
      column(3, downloadButton("report", "Download report")
      ),
      column(3, downloadButton("downloadData", "Download sample")
      ),
      column(3, actionButton("exit", "Quit")
      )
    )
    
  )
  )
  )
)
  


server = function(input, output, session) {
  
# Functions for the control tab ---------------------------------
  
  inputData <- reactive({
    chosenFile <- input$uploadSample
    if (is.null(chosenFile)){
      return(NULL)
    } else {
      controlSample <- read_excel(chosenFile$datapath, sheet=1)
      weibfit <- survreg(Surv(time, cens)~1, data = controlSample, dist = "weibull")
      updateTextInput(session, "lambda2", value = round(as.numeric(1/(exp(weibfit$icoef[1]))), 3))
      updateTextInput(session, "gamma2", value = round(as.numeric(exp(-weibfit$icoef[2])), 3))
      return(list(gamma2 = as.numeric(exp(-weibfit$icoef[2])), lambda2 = as.numeric(1/(exp(weibfit$icoef[1]))), controltime = controlSample$time, controlcens = controlSample$cens))
    }
    
  })
  
  output$recommendedParams <- renderUI({
    if (is.null(inputData())){
      
    } else {
      str1 <- paste0("For your uploaded sample, the best fitting parameters are:")
      str2 <- paste0("Lambda2 = ", round(inputData()$lambda2, 3))
      str3 <- paste0("Gamma2 = ", round(inputData()$gamma2, 3))
      HTML(paste(str1, str2, str3, sep = '<br/>'))
    }
    
  })
  
  output$plotControl <- renderPlot({
    
    controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    controlsurv <- exp(-(input$lambda2*controltime)^input$gamma2)
    controldf <- data.frame(controltime = controltime,
                            controlsurv = controlsurv)
    theme_set(theme_grey(base_size = input$fs))
    p1 <- ggplot(data=controldf, aes(x=controltime, y=controlsurv)) +
      geom_line(colour="blue") + xlab("Time") + ylab("Survival") + ylim(0,1)
    
    print(p1)
    
    if (is.null(inputData())){
      
    } else {
      controlSample <- data.frame(time = inputData()$controltime, cens = inputData()$controlcens)
      km <- survival::survfit(Surv(time, cens)~1, data = controlSample)
      autoplot(km, conf.int = F, surv.colour = "red", xlab = "Time", ylab="Survival")  + 
        geom_line(data = controldf, aes(x = controltime, y = controlsurv), colour = "blue")
      
    }
    
  })
  
# Functions for the eliciting distributions tabs ---------------------------------

  
  limits1 <- reactive({
    eval(parse(text = paste("c(", input$limits1, ")")))
  })
  
  limits2 <- reactive({
    eval(parse(text = paste("c(", input$limits2, ")")))
  })
  
  p1 <- reactive({
    eval(parse(text = paste("c(", input$probs1, ")")))
  })
  
  p2 <- reactive({
    eval(parse(text = paste("c(", input$probs2, ")")))
  })
  
  v1 <- reactive({
    eval(parse(text = paste("c(", input$values1, ")")))
  })
  
  v2 <- reactive({
    eval(parse(text = paste("c(", input$values2, ")")))
  })
  
  m1 <- reactive({
    approx(p1(), v1(), 0.5)$y
  })
  
  m2 <- reactive({
    approx(p2(), v2(), 0.5)$y
  })
  
  myfit1 <- reactive({
    fitdist(vals = v1(), probs = p1(), lower = limits1()[1],
            upper = limits1()[2], 
            tdf = input$tdf1)
  })
  
  myfit2 <- reactive({
    fitdist(vals = v2(), probs = p2(), lower = limits2()[1],
            upper = limits2()[2], 
            tdf = input$tdf2)
  })

 
# Functions for the T tab ---------------------------------
  
  output$distPlot1 <- renderPlot({
    
    
    #d = dist[as.numeric(input$radio1)]
    # dist<-c("hist","normal", "t", "gamma", "lognormal", "logt","beta", "best")
    suppressWarnings(plotfit(myfit1(), d = input$dist1,
                             ql = 0.05, qu = 0.95,
                             xl = limits1()[1], xu = limits1()[2], 
                             fs = input$fs))
    
  })
  
# Functions for the HR tab ---------------------------------
  
  # output$TrialFeedback <- renderUI({
  #   
  #   gamma1 <- input$gamma2
  #   controlcurve <- exp(-(input$lambda2*input$triallength)^input$gamma2)
  #   bigTMedian <- feedback(myfit1(), quantiles = 0.5)$fitted.quantiles[input$dist1][, 1]
  #   HRMedian <- feedback(myfit2(), quantiles = 0.5)$fitted.quantiles[input$dist2][, 1]
  #   lambda1 <- exp((log(HRMedian)/input$gamma2)+log(input$lambda2))
  #   treatmentsurv2 <- exp(-(input$lambda2*bigTMedian)^input$gamma2 - lambda1^gamma1*(input$triallength^gamma1-bigTMedian^gamma1))
  #   
  #   str1 <- paste0("If the trial runs for ", input$triallength, " months, we expect:")
  #   str2 <- paste0(round(controlcurve*100, 1), "% of patients to be alive in the control group")
  #   str3 <- paste0(round(treatmentsurv2*100, 1), "% of patients to be alive in the treatment group")
  #   str4 <- paste0("This means there should be an absolute treatment effect of ", round((treatmentsurv2-controlcurve)*100, 1), "% after ", input$triallength, " months" )
  #   HTML(paste(str1, str2, str3, str4, sep = '<br/>'))
  #   
  # })
  
  output$distPlot2 <- renderPlot({
    
    
    #  dist<-c("hist","normal", "t", "gamma", "lognormal", "logt","beta", "best")
    suppressWarnings(plotfit(myfit2(), d = input$dist2,
                             ql = 0.05, qu = 0.95,
                             xl = limits2()[1], xu = limits2()[2], 
                             fs = input$fs))
    
    
  })

  # HRProportionCalc <- reactive({
  #   
  #   gamma1 <- input$gamma2
  #   controlcurve <- exp(-(input$lambda2*input$triallength)^input$gamma2)
  #   bigTMedian <- feedback(myfit1(), quantiles = 0.5)$fitted.quantiles[input$dist1][, 1]
  #   HR <- seq(0.1, 1, by=0.01)
  #   diff <- rep(NA, length(HR))
  #   for (i in 1:length(HR)){
  #     lambda1 <- exp((log(HR[i])/input$gamma2)+log(input$lambda2))
  #     treatmentsurv2 <- exp(-(input$lambda2*bigTMedian)^input$gamma2 - lambda1^gamma1*(input$triallength^gamma1-bigTMedian^gamma1))
  #     diff[i] <- treatmentsurv2 - controlcurve
  #   }
  #   return(HR[sum(diff>(input$clinicaldiff)/100)])
  #   
  # })
  
  # output$HRProportion <- renderUI({
  #   str1 <- paste0("The probability that HR is less than 1 is: ", feedback(myfit2(), values = 1)$fitted.probabilities[input$dist2][, 1])
  #   str2 <- paste0("For clinical difference, HR needs to be no bigger than: ", HRProportionCalc())
  #   str3 <- paste0("Therefore, probabiliy than HR is lower than target treatment effect: ", feedback(myfit2(), values = HRProportionCalc())$fitted.probabilities[input$dist2][, 1])
  #   HTML(paste(str1, str2, str3, sep = '<br/>'))
  # })
  
# Functions for the Feedback tab ---------------------------------
  
  drawsimlines <- reactive({
    
    gamma1 <- input$gamma2
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    mySample <- data.frame(copulaSample(myfit1(), myfit2(), cp = conc.probs, n = 500, d = c(input$dist1, input$dist2)))
    
    time <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    SimMatrix <- matrix(NA, nrow = 500, ncol=length(time))
    
    for (i in 1:500){
      bigT <- mySample[i,1]
      HR <- mySample[i,2]
      lambda1 <- exp((log(HR)/input$gamma2)+log(input$lambda2))
      
      controltime <- seq(0, bigT, by=0.01)
      controlsurv <- exp(-(input$lambda2*controltime)^input$gamma2)
      
      treatmenttime <- seq(bigT, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
      treatmentsurv <- exp(-(input$lambda2*bigT)^input$gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
      
      timecombined <- c(controltime, treatmenttime)[1:length(time)]
      survcombined <- c(controlsurv, treatmentsurv)[1:length(time)]
      
      
      SimMatrix[i,] <- survcombined
      
    }
    
    lowerbound <- rep(NA, length(time))
    upperbound <- rep(NA, length(time))
    for (j in 1:length(time)){
      lowerbound[j] <- quantile(SimMatrix[,j], 0.1)
      upperbound[j] <- quantile(SimMatrix[,j], 0.9)
    }
    
    return(list(lowerbound=lowerbound, upperbound=upperbound, time=time))
    
  })
  
  output$plotFeedback <- renderPlot({
    
    gamma1 <- input$gamma2
    controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    controlcurve <- exp(-(input$lambda2*controltime)^input$gamma2)
    controldf <- data.frame(controltime = controltime, controlcurve = controlcurve)
    theme_set(theme_grey(base_size = input$fs))
    p1 <- ggplot(data=controldf, aes(x=controltime, y=controlcurve)) +
      geom_line(colour="blue") + xlab("Time") + ylab("Survival") + ylim(0,1)
    
    #legend("topright", legend = c("Same fit before changepoint", "Control", "Treatment"),
    #col=c("green", "blue", "red"), lty=c(1), cex=0.75)
    
    
    bigTMedian <- feedback(myfit1(), quantiles = 0.5)$fitted.quantiles[input$dist1][, 1]
    HRMedian <- feedback(myfit2(), quantiles = 0.5)$fitted.quantiles[input$dist2][, 1]
    lambda1 <- exp((log(HRMedian)/input$gamma2)+log(input$lambda2))
    
    treatmenttime1 <- seq(0, bigTMedian, by=0.01)
    treatmentsurv1 <- exp(-(input$lambda2*treatmenttime1)^input$gamma2)
    treatmenttime1df <- data.frame(treatmenttime1 = treatmenttime1, treatmentsurv1 = treatmentsurv1)
    p1 <-  p1 + geom_line(data = treatmenttime1df, aes(x = treatmenttime1, y = treatmentsurv1), colour = "green") 
    
    
    treatmenttime2 <- seq(bigTMedian, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    treatmentsurv2 <- exp(-(input$lambda2*bigTMedian)^input$gamma2 - lambda1^gamma1*(treatmenttime2^gamma1-bigTMedian^gamma1))
    treatmenttime2df <- data.frame(treatmenttime2 = treatmenttime2, treatmentsurv2 = treatmentsurv2)
    p1 <-  p1 + geom_line(data = treatmenttime2df, aes(x = treatmenttime2, y = treatmentsurv2), colour = "red")
    #scale_color_manual(name='James', breaks=c('Same fit before changepoint', 'Control', 'Treatment'),
    # values=c('Same fit before changepoint'='green', 'Control'='blue', 'Treatment'='red'))
    
    print(p1)
    
    
    addfeedback <- input$showfeedback 
    
    if (!is.null(addfeedback)){
      for (i in 1:length(addfeedback)){
        if (addfeedback[i]=="Median survival line"){
          if (exp(-(input$lambda2*bigTMedian)^input$gamma2)<0.5){
            mediandf <- data.frame(x = seq(0, controltime[sum(controlcurve>0.5)], length=2), y = rep(0.5, 2))
            mediandf1 <- data.frame(x = rep(controltime[sum(controlcurve>0.5)], 2), y = seq(0, 0.5, length=2))
            p1 <- p1 + geom_line(data = mediandf, aes(x = x, y=y), linetype = "dashed") + geom_line(data = mediandf1, aes(x = x, y=y), linetype="dashed") 
          } else {
            mediandf <- data.frame(x = seq(0, treatmenttime2[sum(treatmentsurv2>0.5)], length=2), y = rep(0.5, 2))
            mediandf1 <<- data.frame(x = rep(treatmenttime2[sum(treatmentsurv2>0.5)], 2), y = seq(0, 0.5, length=2))
            mediandf2 <- data.frame(x = rep(controltime[sum(controlcurve>0.5)], 2), y = seq(0, 0.5, length=2))
            p1 <- p1 + geom_line(data = mediandf, aes(x = x, y=y), linetype = "dashed") + geom_line(data = mediandf1, aes(x = x, y=y), linetype="dashed") +
              geom_line(data = mediandf2, aes(x = x, y=y), linetype="dashed")
          }
        } else if (addfeedback[i]=="95% CI for T"){
          p1 <- p1 + geom_point(aes(x = feedback(myfit1(), quantiles = 0.025)$fitted.quantiles[input$dist1][, 1], y = controlcurve[sum(controltime<feedback(myfit1(), quantiles = 0.025)$fitted.quantiles[input$dist1][, 1])]), colour="orange", size = 4) +
            geom_point(aes(x = feedback(myfit1(), quantiles = 0.975)$fitted.quantiles[input$dist1][, 1], y = controlcurve[sum(controltime<feedback(myfit1(), quantiles = 0.975)$fitted.quantiles[input$dist1][, 1])]), colour="orange", size = 4)
        } else if (addfeedback[i]=="CI for Survival Curves (0.1 and 0.9)"){
          #Function to draw simulated lines
          simlineslower <- data.frame(x = drawsimlines()$time, y = drawsimlines()$lowerbound)
          simlinesupper <- data.frame(x = drawsimlines()$time, y = drawsimlines()$upperbound)
          p1 <- p1 + geom_line(data = simlineslower, aes(x=x, y=y), linetype="dashed")+
            geom_line(data = simlinesupper, aes(x=x, y=y), linetype="dashed")
        }
      }
    }
    print(p1)
  })
  
# Functions for the Assurance tab ---------------------------------
  
  
  observe({
    hide("chosenlength")
  })
  
  observeEvent(input$drawSSvMonths, {
   show("chosenlength")
  })
  
  calculateSSvMonths <- eventReactive(input$drawSSvMonths, {
    
    gamma1 <- input$gamma2
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    assnum <- 100
    
    assuranceneeded <- 0.5
    triallengthvec <- floor(seq(input$minlength, input$maxlength, length=10))
    ssneededvec <- rep(NA, length(triallengthvec))
    eventsneededvec <- rep(NA, length(triallengthvec))
    
    
    AssFunc <- function(n1, n2, months){
      
      assvec <- rep(NA, assnum)
      eventvec <- rep(NA, assnum)
      mySample <- data.frame(copulaSample(myfit1(), myfit2(), cp = conc.probs, n = assnum, d = c(input$dist1, input$dist2)))
      
      for (i in 1:assnum){
        
        bigT <- mySample[i,1]
        HR <- mySample[i,2]
        lambda1 <- exp((log(HR)/input$gamma2)+log(input$lambda2))
        
        controldata <- data.frame(time = rweibull(n1, input$gamma2, 1/input$lambda2))
        
        CP <- exp(-(input$lambda2*bigT)^input$gamma2)[[1]]
        u <- runif(n2)
        suppressWarnings(z <- ifelse(u>CP, (1/input$lambda2)*exp(1/input$gamma2*log(-log(u))), exp((1/gamma1)*log(1/(lambda1^gamma1)*(-log(u)-(input$lambda2*bigT)^input$gamma2+lambda1^gamma1*bigT*gamma1)))))
        
        treatmentdata <- data.frame(time = z)
        DataCombined <- data.frame(time = c(controldata$time, treatmentdata$time),
                                   group = c(rep("Control", n1), rep("Treatment", n2)), cens = rep(1, n1+n2))
        
        DataCombined <- DataCombined %>%
          filter(time<months)
        
        test <- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)
        eventvec[i] <- nrow(DataCombined)
      } 
      
      return(list(assurance = sum(assvec)/assnum, events = sum(eventvec)/assnum))
      
    } 
    
    for (i in 1:length(triallengthvec)){
      calcass <- 0
      ratiosum <- input$n1+input$n2
      total <- 100
      n1 <- round((total/ratiosum)*input$n1)
      n2 <- round((total/ratiosum)*input$n2)
      while (calcass<input$chosenassurance){
        calcass <- AssFunc(n1, n2, triallengthvec[i])$assurance
        total <- total + 1000
        n1 <- round((total/ratiosum)*input$n1)
        n2 <- round((total/ratiosum)*input$n2)
      }
      ssneededvec[i] <- total
      eventsneededvec[i] <- AssFunc(n1, n2, triallengthvec[i])$events
    }
    
    
    asssmooth <- loess(ssneededvec~triallengthvec)
    eventssmooth <- loess(eventsneededvec~triallengthvec)
    
    return(list(triallengthvec=triallengthvec, asssmooth=asssmooth, eventssmooth=eventssmooth))
    
  })    


  output$samplesizerequired <- renderUI({
    
    str1 <- paste0("If you run the trial for ", input$chosenlength, " months, you will require a total sample size of ", round(predict(calculateSSvMonths()$asssmooth, newdata = input$chosenlength)))
    str2 <- paste0("This is ", floor(input$n1*(round(predict(calculateSSvMonths()$asssmooth, newdata = input$chosenlength))/(input$n1+input$n2))), " patients in the control group and ", ceiling(input$n2*(round(predict(calculateSSvMonths()$asssmooth, newdata = input$chosenlength))/(input$n1+input$n2))), " patients in the treatment group")
    HTML(paste(str1, str2, sep = '<br/>'))
  })
  

  output$eventsPlot <- renderPlot({
    
    theme_set(theme_grey(base_size = input$fs))
    assurancedf <- data.frame(x = calculateSSvMonths()$triallengthvec, y = predict(calculateSSvMonths()$eventssmooth))
    p1 <- ggplot(data = assurancedf) + geom_line(aes(x = x, y = y), linetype="dashed") + xlab("Trial length") +
      ylab("Sample size needed")
    print(p1) 
    
    
  })

  
  output$samplevmonths <- renderPlot({
    
    theme_set(theme_grey(base_size = input$fs))
    assurancedf <- data.frame(x = calculateSSvMonths()$triallengthvec, y = predict(calculateSSvMonths()$asssmooth))
    p1 <- ggplot(data = assurancedf) + geom_line(aes(x = x, y = y), linetype="dashed") + xlab("Trial length") +
     ylab("Sample size needed")
    print(p1)    
    
    
  })
  
  
  # observeEvent(input$exit, {
  #   stopApp(list(parameter1 = myfit1(), parameter2 = myfit2(), 
  #                cp = input$concProb))
  # }) 
  
  # output$downloadData <- downloadHandler(
  #   filename = "joint-sample.csv",
  #   content = function(file) {
  #     utils::write.csv(df1(), file, row.names = FALSE)
  #   }
  # )
  
  # output$report <- downloadHandler(
  #   filename = function(){switch(input$outFormat,
  #                                html_document = "distributions-report.html",
  #                                pdf_document = "distributions-report.pdf",
  #                                word_document = "distributions-report.docx")},
  #   content = function(file) {
  #     # Copy the report file to a temporary directory before processing it, in
  #     # case we don't have write permissions to the current working dir (which
  #     # can happen when deployed).
  #     tempReport <- file.path(tempdir(), "elicitationShinySummaryBivariate.Rmd")
  #     file.copy(system.file("shinyAppFiles", "elicitationShinySummaryBivariate.Rmd",
  #                           package="SHELF"),
  #               tempReport, overwrite = TRUE)
  #     
  #     # Set up parameters to pass to Rmd document
  #     params <- list(fit1 = myfit1(), fit2 = myfit2(), cp = input$concProb,
  #                    d = c(input$dist1, input$dist2), m1 = m1(), m2 = m2())
  #     
  #     # Knit the document, passing in the `params` list, and eval it in a
  #     # child of the global environment (this isolates the code in the document
  #     # from the code in this app).
  #     rmarkdown::render(tempReport, output_file = file,
  #                       params = params,
  #                       output_format = input$outFormat,
  #                       envir = new.env(parent = globalenv())
  #     )
  #   }
  # )
  
}

shinyApp(ui, server)

###Breaks when T is allowed to be less than 0 - need to fix

