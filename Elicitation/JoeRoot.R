library(readxl)
library(xlsx)
library(tidyverse)
library(ggplot2)
library(zoo)
install.packages("zoo")
JoeRootData <- read_excel("PhD Materials/Elicitation/Joe_Root_Test_Batting_Stats.xlsx")

JoeRootTibble <- JoeRootData %>%
  as_tibble()

RollingCareerAverage <- cumsum(JoeRootTibble$Runs) / seq_along(JoeRootTibble$Runs) 

JoeRootTibble <- JoeRootTibble %>%
  add_column(RollingCareerAverage)


p<-ggplot() +
  geom_bar(data = JoeRootTibble, aes(x=`Test no`, y=Runs),stat="identity", fill =JoeRootTibble$OutColour)  + 
  geom_line(data = JoeRootTibble,aes(x=`Test no`, y =RollingCareerAverage), col="red")
p

JoeRootTibble


write.xlsx(JoeRootTibble, file = "PhD Materials/Elicitation/Joe_Root_Test_Batting_Stats.xlsx", append = FALSE)
library("xlsx") 
install.packages("xlsx")
saveRDS(JoeRootTibble, file = "PhD Materials/Elicitation/Joe_Root_Test_Batting_Stats.rds")

