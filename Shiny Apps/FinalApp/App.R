library(SHELF)
library(shiny)
library(survival)
#' Elicit a bivariate distribution using a Gaussian copula
#' 
#' Opens up a web browser (using the shiny package), from which you can specify
#' judgements, fit distributions, plot the fitted density functions, and plot samples 
#' from the joint distributions. A joint distribution is constructed using a Gaussian
#' copula, whereby the correlation parameter is determined via the elicitation of a 
#' concordance probability (a probability that the two uncertain quantities are either
#' both greater than their medians, or both less than their medians.)
#' 
#' Click on the "Help" tab for instructions. Click the "Quit" button to exit the app and return
#' the results from the \code{fitdist} command. Click "Download report" to generate a report
#' of all the fitted distributions for each uncertain quantity, and "Download sample" to
#' generate a csv file with a sample from the joint distribution.
#'   
#' @return A list, with two objects of class \code{elicitation}, and the 
#' elicited concordance probability. See \code{\link{fitdist}} for details.
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' 
#' \dontrun{
#' 
#' elicit()
#' 
#' }
#' @import shiny
#' @export
elicitBivariate<- function(){
  runApp(list(
    ui = shinyUI(fluidPage(
      
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
                     numericInput("lambda2", "lambda2", value=0.06),
                     numericInput("gamma2", "gamma2", value=0.8)
                   ), 
                   mainPanel = mainPanel(
                     plotOutput("plotControl")
                   )
                 ),
        ),
        tabPanel("Eliciting T",
                 fluidRow(
                   column(4, 
                          textInput("limits1", label = h5("Parameter 1 limits"), 
                                    value = "0, 50")
                   ),
                   column(4,
                          textInput("values1", label = h5("Parameter 1 values"), 
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
                          textInput("limits2", label = h5("Parameter limits"), 
                                    value = "0, 1")
                   ),
                   column(4,
                          textInput("values2", label = h5("Parameter values"), 
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
                     checkboxGroupInput("showfeedback", "Add to plot", choices = c("Median survival line", "Hazard Ratio & 95% CI's", "95% CI for T", "Simulation curves"))
                   ), 
                   mainPanel = mainPanel(
                     plotOutput("plotFeedback")
                   )
                 ),
        ),
                 
          
        tabPanel("Assurance",
                 plotOutput("plotAssurance")
                 
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
    ),
    
    server = function(input, output) {
      
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
      
      output$plotControl <- renderPlot({
        gamma1 <- input$gamma2
        controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
        controlcurve <- exp(-(input$lambda2*controltime)^input$gamma2)
        plot(controltime, controlcurve, type="l", col="blue", xlab="Time", ylab="Survival", main="Control curve", ylim=c(0,1))
        
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
      
      # df1 <- reactive({
      #   conc.probs <- matrix(0, 2, 2)
      #   conc.probs[1, 2] <- 0.5
      #   data.frame(copulaSample(myfit1(), myfit2(), cp = conc.probs,
      #                           n = input$sampleSize,
      #                           d = c(input$dist1, input$dist2)))
      # })
      
      drawsimlines <- reactive({
    
        gamma1 <- input$gamma2
        linelist <- list()
        for (i in 1:10){
          bigT <- rnorm(1, mean = as.numeric(myfit1()$Normal[1]), sd = as.numeric(myfit1()$Normal[2]))
          HR <- rbeta(1, as.numeric(myfit2()$Beta[1]), as.numeric(myfit2()$Beta[2]))
          lambda1 <- exp((log(HR)/input$gamma2)+log(input$lambda2))
          treatmenttime <- seq(bigT, 120, by=0.01)
          treatmentsurv <- exp(-(input$lambda2*bigT)^input$gamma2 - lambda1^gamma1*(treatmenttime^gamma1-bigT^gamma1))
          linelist[[(i*2)-1]] <- treatmenttime
          linelist[[i*2]] <- treatmentsurv
        }
        list(linelist = linelist)
      })
      
      output$plotFeedback <- renderPlot({
  
        gamma1 <- input$gamma2
        controltime <- seq(0, exp((1.527/input$gamma2)-log(input$lambda2))*1.1, by=0.01)
        controlcurve <- exp(-(input$lambda2*controltime)^input$gamma2)
        plot(controltime, controlcurve, type="l", col="blue", xlab="Time", ylab="Survival", main="Elicitation of treatment curve")
        legend("topright", legend = c("Same fit before changepoint", "Control", "Treatment"),
               col=c("green", "blue", "red"), lty=c(1), cex=0.75)
        
    
        
        # theta <<- feedback(myfit1(), quantiles = 0.5)
        # print(theta$fitted.quantiles[input$dist1][, 1])
        
        bigTMedian <- feedback(myfit1(), quantiles = 0.5)$fitted.quantiles[input$dist1][, 1]
        HRMedian <- feedback(myfit2(), quantiles = 0.5)$fitted.quantiles[input$dist2][, 1]
        lambda1 <- as.numeric(exp((log(HRMedian)/input$gamma2)+log(input$lambda2)))
        
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
                lines(seq(bigTMedian, 120,length=2), rep(HRMedian, 2))
                #Conf intervals
                #lines(seq(bigTMedian, 120,length=2), rep(qbeta(0.025, as.numeric(HRFit[1]), as.numeric(HRFit[2])), 2), lty=2)
                #lines(seq(bigTMedian, 120,length=2), rep(qbeta(0.975, as.numeric(HRFit[1]), as.numeric(HRFit[2])), 2), lty=2)
              } else if (addfeedback[i]=="95% CI for T"){
                points(qnorm(0.025, mean = bigTMedian, sd = as.numeric(myfit1()$Normal[2])), controlcurve[sum(controltime<qnorm(0.025, mean = bigTMedian, sd = as.numeric(myfit1()$Normal[2])))], cex=1.5, col="orange", pch=19)
                points(qnorm(0.975, mean = bigTMedian, sd = as.numeric(myfit1()$Normal[2])), controlcurve[sum(controltime<qnorm(0.975, mean = bigTMedian, sd = as.numeric(myfit1()$Normal[2])))], cex=1.5, col="orange", pch=19)
              } else if (addfeedback[i]=="Simulation curves"){
                Curves <- drawsimlines()$linelist
                for (i in 1:10){
                  lines(Curves[[(i*2)-1]], Curves[[i*2]], col="purple", lwd=0.25, lty=2)
                }
              }
            }
          }
      })
      
      output$plotAssurance <- renderPlot({
        
        AssFunc <- function(n1, n2){
          assnum <- 100
          assvec <- rep(NA, assnum)
        
        assvec <- rep(NA, 100)
        
        for (i in 1:100){
        
          gamma1 <- input$gamma2
          bigT <- rnorm(1, mean = as.numeric(myfit1()$Normal[1]), sd = as.numeric(myfit1()$Normal[2]))
          HR <- rbeta(1, as.numeric(myfit2()$Beta[1]), as.numeric(myfit2()$Beta[2]))
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
        
        return(sum(assvec)/100)
        
        }
        
        n1vec <- seq(10, 500, length=50)
        n2vec <- seq(10, 500, length=50)
        assvec <- rep(NA, 50)
        
        for (i in 1:50){
          assvec[i] <- AssFunc(n1vec[i], n2vec[i])
        }
        
        sumvec <- n1vec+n2vec
        asssmooth <- loess(assvec~sumvec)
        plot(sumvec, predict(asssmooth), type="l", lty=2, ylim=c(0,1),
             xlab="Total sample size", ylab="Assurance")
       
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
  ))
}

elicitBivariate()

