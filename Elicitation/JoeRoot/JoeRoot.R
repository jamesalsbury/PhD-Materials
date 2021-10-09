library(readxl)
library(xlsx)
library(tidyverse)
library(ggplot2)
library(zoo)
library(tidyverse)

#Import dataset
JoeRootData <- readRDS("Elicitation/Joe_Root_Test_Batting_Stats.rds")

mean(JoeRootData$Runs)
#Mean of 46.39

median(JoeRootData$Runs)
#Median is lower at 29.5

#Plotting the number of runs as a bar chart - blue = out, yellow = not out
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

#PDF and CDF plots are found in the folder
#Think 400 may be a little high, gives 0.9 quantile of being 265, which is very high for 1 in 10 chance
#But needs to be high as is still plausible (just)

