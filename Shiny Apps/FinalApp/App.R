library(SHELF)
library(shiny)
library(survival)
library(shinydashboard)
library(readxl)
library(rsconnect)


ui <- fluidPage(
  
  # Application title
  titlePanel("Delayed Treatment Effects - Weibull parameterisation"),
  
  # sidebarLayout(
  mainPanel(tags$style(type="text/css",
                       ".shiny-output-error { visibility: hidden; }",
                       ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  tabsetPanel(
    tabPanel("Control", 
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 fileInput("uploadSample", "Upload your control sample"),
                 numericInput("lambda2", "lambda2", value=0.06),
                 numericInput("gamma2", "gamma2", value=0.8)
               ), 
               mainPanel = mainPanel(
                 plotOutput("plotControl"),
                 htmlOutput("recommendedParams")
               )
             ),
    ),
    tabPanel("Eliciting T",
             fluidRow(
               column(4, 
                      textInput("limits1", label = h5("T limits"), 
                                value = "0, 50")
               ),
               column(4,
                      textInput("values1", label = h5("T values"), 
                                value = "5, 6, 7")
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
    tabPanel("Eliciting HR",
             fluidRow(
               column(4, 
                      textInput("limits2", label = h5("HR limits"), 
                                value = "0, 1")
               ),
               column(4,
                      textInput("values2", label = h5("HR values"), 
                                value = "0.5, 0.6, 0.75")
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
             
             plotOutput("distPlot2")
    ),
    
    tabPanel("Feedback", 
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 checkboxGroupInput("showfeedback", "Add to plot", choices = c("Median survival line", "Hazard Ratio & 95% CI's", "95% CI for T", "CI for Survival Curves (0.1 and 0.9)"))
               ), 
               mainPanel = mainPanel(
                 plotOutput("plotFeedback")
               )
             ),
    ),
    
    
    tabPanel("Assurance",
             sidebarLayout(
               sidebarPanel = sidebarPanel(
                 numericInput("maxss", "Maximum sample size?", value=1000),
                 box(width = 10, title = "Ratio of patients in each group?",
                     splitLayout(
                       numericInput("n1", "Control", value=1, min=1),
                       numericInput("n2", "Treatment", value=1, min=1)
                     )
                 ),
                 actionButton("drawassurance", "Plot assurance line"),
                 numericInput("samplesize", "Assurance at sample size:", value=100),
               ), 
               mainPanel = mainPanel(
                 plotOutput("plotAssurance"),
                 htmlOutput("assuranceSS")
               )
             ),
             
    )
    
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

server = function(input, output, session) {
  
  # Hack to avoid CRAN check NOTE
  
  X1 <- X2 <- xpos <- ypos <- hjustvar <- vjustvar <- annotateText <- NULL
  
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
      
    } else{
      str1 <- paste0("For your uploaded sample, the best fitting parameters are:")
      str2 <- paste0("Lambda2 = ", round(inputData()$lambda2, 3))
      str3 <- paste0("Gamma2 = ", round(inputData()$gamma2, 3))
      HTML(paste(str1, str2, str3, sep = '<br/>'))
    }
    
  })
  
  output$plotControl <- renderPlot({
    
    controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    controlsurv <- exp(-(input$lambda2*controltime)^input$gamma2)
    plot(controltime, controlsurv, type="l", col="blue", xlab="Time", ylab="Survival", ylim=c(0,1))
    
    if (is.null(inputData())){
      
    } else{
      controlSample <- data.frame(time = inputData()$controltime, cens = inputData()$controlcens)
      km <- survfit(Surv(time, cens)~1, data = controlSample)
      lines(km, conf.int = F)
    }
    
    
    legend("topright", legend=c("Kaplan-Meier", "Weibull"), col=c("black", "blue"), lty=1)
    
    
  })
  
  
  output$distPlot1 <- renderPlot({
    
    
    #d = dist[as.numeric(input$radio1)]
    # dist<-c("hist","normal", "t", "gamma", "lognormal", "logt","beta", "best")
    suppressWarnings(plotfit(myfit1(), d = input$dist1,
                             ql = 0.05, qu = 0.95,
                             xl = limits1()[1], xu = limits1()[2], 
                             fs = input$fs))
    
  })
  
  
  output$distPlot2 <- renderPlot({
    
    
    #  dist<-c("hist","normal", "t", "gamma", "lognormal", "logt","beta", "best")
    suppressWarnings(plotfit(myfit2(), d = input$dist2,
                             ql = 0.05, qu = 0.95,
                             xl = limits2()[1], xu = limits2()[2], 
                             fs = input$fs))
    
    
  })
  
  
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
  
  
  calculateAssurance <- eventReactive(input$drawassurance, {
    
    gamma1 <- input$gamma2
    conc.probs <- matrix(0, 2, 2)
    conc.probs[1, 2] <- 0.5
    assnum <- 100
    
    AssFunc <- function(n1, n2){
      
      assvec <- rep(NA, assnum)
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
        test <- survdiff(Surv(time, cens)~group, data = DataCombined)
        assvec[i] <- test$chisq > qchisq(0.95, 1)
      }
      
      return(sum(assvec)/assnum)
      
    }
    
    ratiosum <- input$maxss/(input$n1+input$n2)
    
    n1vec <- floor(seq(10, ratiosum*input$n1, length=50))
    n2vec <- floor(seq(10, ratiosum*input$n2, length=50))
    assvec <- rep(NA, 50)
    
    for (i in 1:50){
      assvec[i] <- AssFunc(n1vec[i], n2vec[i])
    }
    
    sumvec <- n1vec+n2vec
    asssmooth <- loess(assvec~sumvec)
    
    return(list(sumvec=sumvec, asssmooth=asssmooth))
    
    
  })
  
  
  output$plotFeedback <- renderPlot({
    
    gamma1 <- input$gamma2
    controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    controlcurve <- exp(-(input$lambda2*controltime)^input$gamma2)
    plot(controltime, controlcurve, type="l", col="blue", xlab="Time", ylab="Survival", main="Elicitation of treatment curve")
    legend("topright", legend = c("Same fit before changepoint", "Control", "Treatment"),
           col=c("green", "blue", "red"), lty=c(1), cex=0.75)
    
    
    bigTMedian <- feedback(myfit1(), quantiles = 0.5)$fitted.quantiles[input$dist1][, 1]
    HRMedian <- feedback(myfit2(), quantiles = 0.5)$fitted.quantiles[input$dist2][, 1]
    lambda1 <- exp((log(HRMedian)/input$gamma2)+log(input$lambda2))
    
    treatmenttime1 <- seq(0, bigTMedian, by=0.01)
    treatmentsurv1 <- exp(-(input$lambda2*treatmenttime1)^input$gamma2)
    lines(treatmenttime1, treatmentsurv1, col="green")
    
    treatmenttime2 <- seq(bigTMedian, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
    treatmentsurv2 <- exp(-(input$lambda2*bigTMedian)^input$gamma2 - lambda1^gamma1*(treatmenttime2^gamma1-bigTMedian^gamma1))
    lines(treatmenttime2, treatmentsurv2, col="red")
    
    addfeedback <- input$showfeedback 
    
    if (!is.null(addfeedback)){
      for (i in 1:length(addfeedback)){
        if (addfeedback[i]=="Median survival line"){
          lines(seq(-1, treatmenttime2[sum(treatmentsurv2>0.5)], length=2), rep(0.5, 2), lty=3)
          lines(rep(treatmenttime2[sum(treatmentsurv2>0.5)], 2), seq(-1, 0.5, length=2), lty=3)
          lines(rep(controltime[sum(controlcurve>0.5)], 2), seq(-1, 0.5, length=2), lty=3)
        } else if (addfeedback[i]=="Hazard Ratio & 95% CI's"){
          #Top line
          lines(seq(0, bigTMedian, length=2), rep(1, 2))
          #Vertical line
          lines(rep(bigTMedian, 2), seq(HRMedian, 1, length=2))
          #Bottom line
          lines(seq(bigTMedian, exp((1.527/input$gamma2)-log(input$lambda2))*1.1,length=2), rep(HRMedian, 2))
          #Conf intervals
          lines(seq(bigTMedian, exp((1.527/input$gamma2)-log(input$lambda2))*1.1,length=2), rep(feedback(myfit2(), quantiles = 0.975)$fitted.quantiles[input$dist2][, 1], 2), lty=2)
          lines(seq(bigTMedian, exp((1.527/input$gamma2)-log(input$lambda2))*1.1,length=2), rep(feedback(myfit2(), quantiles = 0.025)$fitted.quantiles[input$dist2][, 1], 2), lty=2)
        } else if (addfeedback[i]=="95% CI for T"){
          points(feedback(myfit1(), quantiles = 0.025)$fitted.quantiles[input$dist1][, 1], controlcurve[sum(controltime<feedback(myfit1(), quantiles = 0.025)$fitted.quantiles[input$dist1][, 1])], cex=1.5, col="orange", pch=19)
          points(feedback(myfit1(), quantiles = 0.975)$fitted.quantiles[input$dist1][, 1], controlcurve[sum(controltime<feedback(myfit1(), quantiles = 0.975)$fitted.quantiles[input$dist1][, 1])], cex=1.5, col="orange", pch=19)
        } else if (addfeedback[i]=="CI for Survival Curves (0.1 and 0.9)"){
          #Function to draw simulated lines
          lines(drawsimlines()$time, drawsimlines()$lowerbound, lty=2)
          lines(drawsimlines()$time, drawsimlines()$upperbound, lty=2)
          
        }
      }
    }
  })
  
  output$plotAssurance <- renderPlot({
    
    plot(calculateAssurance()$sumvec, predict(calculateAssurance()$asssmooth), type="l", lty=2, ylim=c(0,1),
         xlab="Total sample size", ylab="Assurance")
    
  })
  
  output$assuranceSS <- renderUI({
    
    paste0("With a sample size of ", input$samplesize, " assurance is: ", round(predict(calculateAssurance()$asssmooth, newdata = input$samplesize), 2))
    
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


