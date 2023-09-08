library(tidyverse)
library(pracma)

thresholds <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
data_list <- list()

for (threshold in thresholds) {
  file_path <- paste("~/R Projects/PhD-Materials/Papers/InterimDTEPaper/Data/HRData/HR", threshold, ".rds", sep = "")
  data_list[[as.character(threshold)]] <- readRDS(file_path)
  
  # Save the HR dataset as a separate object (e.g., HR0.5, HR0.6, etc.)
  assign(paste("HR", threshold, sep = ""), data_list[[as.character(threshold)]])
}

# Function to calculate TP, FP, FN, TN, FPR, and TPR
calculate_metrics <- function(data, rule, HRIA1, HRIA2) {
  
  Stopped <- data[[paste0(rule, "OHR1")]] > HRIA1 | data[[paste0(rule, "OHR2")]] > HRIA2
  
  data <- cbind(data, Stopped)
  
    SucHR <- data %>%
      filter(FinalOutcome == "Successful trial")
    
    NotSuc <- data %>%
      filter(FinalOutcome != "Successful trial")
    
    TP <- sum(SucHR$Stopped==F)
    FP <- sum(SucHR$Stopped==T)
    FN <- sum(NotSuc$Stopped==F)
    TN <-  sum(NotSuc$Stopped==T)
  
  return(list(TP = TP, FP = FP, TN = TN, FN = FN))
}

ROCFunction <- function(Rule, HRIA1, HRIA2){
  TP <- 0
  FP <- 0
  TN <- 0
  FN <- 0
  
  for (threshold in thresholds){
    filtered_data <- get(paste("HR", threshold, sep = ""))
    outcome <- calculate_metrics(filtered_data, Rule, HRIA1, HRIA2)
    TP <- TP + outcome$TP
    FP <- FP + outcome$FP
    TN <- TN + outcome$TN
    FN <- FN + outcome$FN
  }
  
  FPR <- FP/(FP+TN)
  TPR <- TP/(TP+FN)
  
  return(list(FPR = FPR, TPR = TPR))
  
  
}

#Changing HRIA2
HRIA1 <- seq(0.5, 1.5, by=0.1)

AUCFunc <- function(HRIA2){
  WieandRule <- sapply(X = HRIA1, FUN = ROCFunction, Rule = "Prop", 0.9)
  WieandRuleFPR <- c(1, unlist(WieandRule[1,]), 0)
  WieandRuleTPR <- c(1, unlist(WieandRule[2,]), 0)
  abs(trapz(WieandRuleFPR, WieandRuleTPR))
}

optimize(AUCFunc, lower = 0.3, upper = 2, maximum = T)

#Changing HRIA1
HRIA2 <- seq(0.5, 1.5, by=0.1)

AUCFunc <- function(HRIA1){
  WieandRule <- sapply(X = HRIA2, FUN = ROCFunction, Rule = "Wieand", HRIA1)
  WieandRuleFPR <- c(1, unlist(WieandRule[1,]), 0)
  WieandRuleTPR <- c(1, unlist(WieandRule[2,]), 0)
  abs(trapz(WieandRuleFPR, WieandRuleTPR))
}

optimize(AUCFunc, lower = 0.3, upper = 2, maximum = T)


 plot(WieandRuleFPR, WieandRuleTPR, xlim = c(0,1), ylim=c(0,1), type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate")
# lines(OBFRuleFPR, OBFRuleTPR, xlim = c(0,1), ylim=c(0,1), col = "red")
# lines(PropRuleFPR, PropRuleTPR, xlim = c(0,1), ylim=c(0,1), col = "blue")
# legend("bottomright", legend = c("Wieand", "OBF", "Proposed"), col = c("black", "red", "blue"), lty = 1)
# curve(x^1, 0, 1, add = T, lty = 2)


OBFRule <- sapply(X = HRIA1, FUN = ROCFunction, Rule = "OBF", 1)
PropRule <- sapply(X = HRIA1, FUN = ROCFunction, Rule = "Prop", 1)
OBFRuleFPR <- c(1, unlist(OBFRule[1,]), 0)
OBFRuleTPR <- c(1, unlist(OBFRule[2,]), 0)
PropRuleFPR <- c(1, unlist(PropRule[1,]), 0)
PropRuleTPR <- c(1, unlist(PropRule[2,]), 0)
abs(trapz(OBFRuleFPR, OBFRuleTPR))
abs(trapz(PropRuleFPR, PropRuleTPR))

#################
calculate_metrics_BPP <- function(data, BPP1, BPP2) {
  
  Stopped <- data$WBPP1 < BPP1 | data$WBPP2 < BPP2
  
  data <- cbind(data, Stopped)
  
  SucHR <- data %>%
    filter(FinalOutcome == "Successful trial")
  
  NotSuc <- data %>%
    filter(FinalOutcome != "Successful trial")
  
  TP <- sum(SucHR$Stopped==F)
  FP <- sum(SucHR$Stopped==T)
  FN <- sum(NotSuc$Stopped==F)
  TN <-  sum(NotSuc$Stopped==T)
  
  return(list(TP = TP, FP = FP, TN = TN, FN = FN))
}

ROCFunction_BPP <- function(BPP1, BPP2){
  TP <- 0
  FP <- 0
  TN <- 0
  FN <- 0
  
  for (threshold in thresholds){
    filtered_data <- get(paste("HR", threshold, sep = ""))
    outcome <- calculate_metrics_BPP(filtered_data, BPP1, BPP2)
    TP <- TP + outcome$TP
    FP <- FP + outcome$FP
    TN <- TN + outcome$TN
    FN <- FN + outcome$FN
  }
  
  FPR <- FP/(FP+TN)
  TPR <- TP/(TP+FN)
  
  return(list(FPR = FPR, TPR = TPR))
  
}

#Changing BPP2
BPP1 <- seq(0, 1, by=0.05)

AUCFunc_BPP <- function(BPP2){
  BPPRule <- sapply(X = BPP1, FUN = ROCFunction_BPP, BPP2)
  BPPRuleFPR <- c(0, unlist(BPPRule[1,]), 1)
  BPPRuleTPR <- c(0, unlist(BPPRule[2,]), 1)
  abs(trapz(BPPRuleFPR, BPPRuleTPR))
}

optimize(AUCFunc_BPP, lower = 0.2, upper = 1, maximum = T)

BPPRule <- sapply(X = BPP1, FUN = ROCFunction_BPP, 0.3889344)
BPPRuleFPR <- c(0, unlist(BPPRule[1,]), 1)
BPPRuleTPR <- c(0, unlist(BPPRule[2,]), 1)

lines(BPPRuleFPR, BPPRuleTPR, col = "blue")

#Changing BPP1
BPP2 <- seq(0.2, 1, by=0.05)

AUCFunc_BPP <- function(BPP1){
  BPPRule <- sapply(X = BPP2, FUN = ROCFunction_BPP, BPP1)
  BPPRuleFPR <- c(0, unlist(BPPRule[1,]), 1)
  BPPRuleTPR <- c(0, unlist(BPPRule[2,]), 1)
  abs(trapz(BPPRuleFPR, BPPRuleTPR))
}

optimize(AUCFunc_BPP, lower = 0.2, upper = 1, maximum = T)





