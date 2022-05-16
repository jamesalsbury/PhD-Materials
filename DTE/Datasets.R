#Data sets that Satrajit sent me for DTE

#First data set
DTEDataSet1 <- read.csv(file = "DTE/CM017PFS.csv")

TreatmentData <- data.frame(time = DTEDataSet1$NIVO_Time[1:135], cens = DTEDataSet1$NIVO_Event_Flag[1:135])
 
ControlData <- data.frame(time = DTEDataSet1$DOX_Time, cens = DTEDataSet1$DOX_Event_Flag)

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")




#Second data set
DTEDataSet2 <- read.csv(file = "DTE/CM141OS.csv")

DTEDataSet2$SOC_Event

TreatmentData <- data.frame(time = DTEDataSet2$Nivolumab_Time, cens = DTEDataSet2$Nivolumab_Event)

ControlData <- data.frame(time = DTEDataSet2$SOC_Time[1:121], cens = DTEDataSet2$SOC_Event[1:121])

TreatmentFit <- survfit(Surv(time, cens)~1, data = TreatmentData)

ControlFit <- survfit(Surv(time, cens)~1, data = ControlData)

plot(TreatmentFit, conf.int = F, col="red")

lines(ControlFit, conf.int = F, col="blue")
