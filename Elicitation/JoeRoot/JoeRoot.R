library(tidyverse)
#Original data set can be found: https://stats.espncricinfo.com/ci/engine/player/303669.html?class=1;template=results;type=allround;view=match

#Import dataset
JoeRootData <- readRDS("Elicitation/JoeRoot/Joe_Root_Test_Batting_Stats.rds")

mean(JoeRootData$Runs)
#Mean of 46.39

median(JoeRootData$Runs)
#Median is lower at 29.5

#Plotting the number of runs as a bar chart - blue = out, yellow = not out
#Also rolling career average is plotted as a red line
p<-ggplot() +
  geom_bar(data = JoeRootData, aes(x=`Test no`, y=Runs),stat="identity", fill =JoeRootData$OutColour)  + 
  geom_line(data = JoeRootData,aes(x=`Test no`, y =RollingCareerAverage), col="red")
p

#Extracting only Ashes tests
Ashes <- JoeRootData %>%
  filter(`Ashes?`=="Y")

mean(Ashes$Runs)
#36.82609
#A little lower mean than other tests - to be expected

median(Ashes$Runs)
#Again, lower at 19.5

#Fitting survival curves to data
library(survival)

#Fitting a survival curve for all tests
f1 <- survfit(Surv(Runs, Cens) ~ 1, data = JoeRootData)
plot(f1)

#Just Ashes tests
f2 <- survfit(Surv(Runs, Cens) ~ 1, data = Ashes)
plot(f2)

#Home or away doesn't really make a difference
f3 <- survfit(Surv(Runs, Cens) ~ Home, data = JoeRootData)
plot(f3, col=c("red", "blue"))

#Scores lower in ashes tests than other tests
f4 <- survfit(Surv(Runs, Cens) ~ Ashes, data = JoeRootData)
plot(f4, col=c("red", "blue"))

#We are interested in innings 1 and 2, seems to score higher in 2
f5 <- survfit(Surv(Runs, Cens) ~ Inns, data = JoeRootData)
plot(f5, col=c("red", "blue", "yellow", "black"))

#Again, home or away doesn't really make a difference in Ashes tests
f6 <- survfit(Surv(Runs, Cens) ~ Home, data = Ashes)
plot(f6, col=c("red", "blue"))


#So main findings:
#Home/Away does not matter, scores less in Ashes tests, difference in innings 1 and 2

f7 <- survfit(Surv(Runs, Cens) ~ Inns, data = Ashes)
plot(f7, col=c("red", "blue", "yellow", "black"))

#However, difference is less clear in Ashes tests

#So we will just consider data from all the Ashes tests, aka f2

plot(f2)

Brisbane <- Ashes %>%
  filter(Ground=="Brisbane")

mean(Brisbane$Runs)

#The first test is in Brisbane, his first innings scores there was previously 15 and 2

#Using survival curve f2, we can see that the median value is ~20
#Similarly, UQ = 9, LQ = 63
#Lowest possible value is 0
#Highest value 400? Very improbable, but guess that's the point

#Can plug these values into SHELF

library(SHELF)

SHELF::elicit()


#We need to choose a distribution that only takes non-negative real values
#Possibly Gamma, log-Normal
#Can look at quantiles to determine whether suitable

#For gamma(0.657, 0.0153) the 0.1 quantile is 1.7 which seems low
#The 0.9 quantile is 109 which seems plausible

#For log-Normal(3.09, 1.46) the 0.1 quantile is 3.41 which again is low but a bit more plausible
#(he scored 2 in 2013 in Brisbane)
#The 0.9 quantile is 143 which seems plausible (just)

#I would choose the log-Normal(3.09, 1.46) distribution here


x <- seq(0,400, by=0.01)
y <- dlnorm(x,3.09,1.46)
plot(x,y, type="l")

#Perhaps a truncated Normal/student-t distribution could be a better choice?

















