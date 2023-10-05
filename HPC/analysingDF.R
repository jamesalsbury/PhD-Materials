# 
par(mfrow = c(3,4))

for (i in 1:length(ScenarioList)){
  plot(outcomeDF$`No IA Power`[i], outcomeDF$`No IA Duration`[i], pch = 19, col = "red", ylim= c(0, 60), xlab = "Power", ylab = "Duration", xlim = c(0, 1),
       main = paste0("HR ", ScenarioList[[i]]$HR1, " for ", ScenarioList[[i]]$T1, " months, then ", ScenarioList[[i]]$HR2, ", R: ", ScenarioList[[i]]$recTime))
  points(outcomeDF$`Wieand Power`[i], outcomeDF$`Wieand Duration`[i], pch = 19, col = "blue")
  points(outcomeDF$`OBF Power`[i], outcomeDF$`OBF Duration`[i], pch = 19, col = "yellow")
  points(outcomeDF$`Prop Power`[i], outcomeDF$`Prop Duration`[i], pch = 19, col = "green")
 # legend("topleft", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("red", "blue", "yellow", "green"), pch = 19)
}

par(mfrow = c(3,4))

for (i in 1:length(ScenarioList)){
  plot(outcomeDF$`No IA Power`[i], outcomeDF$`No IA SS`[i], pch = 19, col = "red", ylim= c(300, 700), xlab = "Power", ylab = "Sample size", xlim = c(0, 1), 
       main = paste0("HR ", ScenarioList[[i]]$HR1, " for ", ScenarioList[[i]]$T1, " months, then ", ScenarioList[[i]]$HR2, ", R: ", ScenarioList[[i]]$recTime))
  points(outcomeDF$`Wieand Power`[i], outcomeDF$`Wieand SS`[i], pch = 19, col = "blue")
  points(outcomeDF$`OBF Power`[i], outcomeDF$`OBF SS`[i], pch = 19, col = "yellow")

  points(outcomeDF$`Prop Power`[i], outcomeDF$`Prop SS`[i], pch = 19, col = "green")
 

  #legend("bottomright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("red", "blue", "yellow", "green"), pch = 19)
}
# 
# 
# par(mfrow = c(3,4))
# 
# for (i in 1:length(ScenarioList)){
#   plot(outcomeDF$`No IA Duration`[i], outcomeDF$`No IA SS`[i], pch = 19, col = "red", ylim= c(300, 700), 
#        xlab = "Duration", ylab = "Sample size", xlim = c(10, 60), main = paste0("HR ", ScenarioList[[i]]$HR1, " for ", ScenarioList[[i]]$T1, " months, then ", ScenarioList[[i]]$HR2))
#   points(outcomeDF$`OBF Duration`[i], outcomeDF$`OBF SS`[i], pch = 19, col = "yellow")
#   
#   points(outcomeDF$`Prop Duration`[i], outcomeDF$`Prop SS`[i], pch = 19, col = "green")
#   points(outcomeDF$`Wieand Duration`[i], outcomeDF$`Wieand SS`[i], pch = 19, col = "blue")
#   
#   legend("bottomright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("red", "blue", "yellow", "green"), pch = 19)
# }





# hist(outcomeDF$NoIAPowerRank)
# hist(outcomeDF$WieandPowerRank)
# hist(outcomeDF$OBFPowerRank)
# hist(outcomeDF$PropPowerRank)
# 
# hist(outcomeDF$NoIADurationRank)
# hist(outcomeDF$WieandDurationRank)
# hist(outcomeDF$OBFDurationRank)
# hist(outcomeDF$PropDurationRank)
# 
# hist(outcomeDF$NoIASSRank)
# hist(outcomeDF$WieandSSRank)
# hist(outcomeDF$OBFSSRank)
# hist(outcomeDF$PropSSRank)


#saveRDS(outcomeDF, file = "HPC/outcomeDF.rds")

#Power vs sample size 
plot(outcomeDF$`No IA Power`, outcomeDF$`No IA SS`, xlab = "Power", ylab = "Sample size", 
     pch = 19, col = "black", ylim = c(200,700), xlim = c(0,1))
points(outcomeDF$`Wieand Power`, outcomeDF$`Wieand SS`, pch = 19, col = "green")
points(outcomeDF$`OBF Power`, outcomeDF$`OBF SS`, pch = 19, col = "red")
points(outcomeDF$`Prop Power`, outcomeDF$`Prop SS`, pch = 19, col = "blue")
legend("bottomright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("black", "green", "red", "blue"), pch = 19)


#Power vs duration 
plot(outcomeDF$`No IA Power`, outcomeDF$`No IA Duration`, xlab = "Power", ylab = "Duration", 
     pch = 19, col = "black", xlim = c(0,1), ylim = c(0,55))
points(outcomeDF$`Wieand Power`, outcomeDF$`Wieand Duration`, pch = 19, col = "green")
points(outcomeDF$`OBF Power`, outcomeDF$`OBF Duration`, pch = 19, col = "red")
points(outcomeDF$`Prop Power`, outcomeDF$`Prop Duration`, pch = 19, col = "blue")
legend("bottomright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("black", "green", "red", "blue"), pch = 19)

#Power vs sample size (rankings)
plot(outcomeDF$NoIAPowerRank, outcomeDF$NoIASSRank, xlab = "Power", ylab = "Sample size", 
     pch = 19, col = "black", ylim = c(0,6), xlim = c(0,5))
points(outcomeDF$WieandPowerRank, outcomeDF$WieandSSRank, pch = 19, col = "green")
points(outcomeDF$OBFPowerRank, outcomeDF$OBFSSRank, pch = 19, col = "red")
points(outcomeDF$PropPowerRank, outcomeDF$PropSSRank, pch = 19, col = "blue")
legend("topright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("black", "green", "red", "blue"), pch = 19)

#Power vs duration (rankings)
plot(outcomeDF$NoIAPowerRank, outcomeDF$NoIADurationRank, xlab = "Power", ylab = "Duration", 
     pch = 19, col = "black", ylim = c(0,6), xlim = c(0,5))
points(outcomeDF$WieandPowerRank, outcomeDF$WieandDurationRank, pch = 19, col = "green")
points(outcomeDF$OBFPowerRank, outcomeDF$OBFDurationRank, pch = 19, col = "red")
points(outcomeDF$PropPowerRank, outcomeDF$PropDurationRank, pch = 19, col = "blue")
legend("topright", legend = c("No IA", "Wieand", "OBF", "Proposed"), col = c("black", "green", "red", "blue"), pch = 19)

#Power vs duration
df1_summary = ddply(
  outcomeDF,
  .(NoIAPowerRank, NoIADurationRank),
  summarize,
  count=length(NoIADurationRank)
)


df2_summary = ddply(
  outcomeDF,
  .(WieandPowerRank, WieandDurationRank),
  summarize,
  count=length(WieandDurationRank)
)

df3_summary = ddply(
  outcomeDF,
  .(OBFPowerRank, OBFDurationRank),
  summarize,
  count=length(OBFDurationRank)
)

df4_summary = ddply(
  outcomeDF,
  .(PropPowerRank, PropDurationRank),
  summarize,
  count=length(PropDurationRank)
)

ggplot(df1_summary, aes(NoIAPowerRank, NoIADurationRank, size=count)) +
  geom_point(col = "red") + 
  geom_point(data = df2_summary, aes(WieandPowerRank, WieandDurationRank, size=count), col = "blue") +
  geom_point(data = df3_summary, aes(OBFPowerRank, OBFDurationRank, size=count), col = "green") +
  geom_point(data = df4_summary, aes(PropPowerRank, PropDurationRank, size=count), col = "yellow")

#Power vs sample size
df1_summary = ddply(
  outcomeDF,
  .(NoIAPowerRank, NoIASSRank),
  summarize,
  count=length(NoIASSRank)
)


df2_summary = ddply(
  outcomeDF,
  .(WieandPowerRank, WieandSSRank),
  summarize,
  count=length(WieandSSRank)
)

df3_summary = ddply(
  outcomeDF,
  .(OBFPowerRank, OBFSSRank),
  summarize,
  count=length(OBFSSRank)
)

df4_summary = ddply(
  outcomeDF,
  .(PropPowerRank, PropSSRank),
  summarize,
  count=length(PropSSRank)
)

ggplot(df1_summary, aes(NoIAPowerRank, NoIASSRank, size=count)) +
  geom_point(col = "red") + 
  geom_point(data = df2_summary, aes(WieandPowerRank, WieandSSRank, size=count), col = "blue") +
  geom_point(data = df3_summary, aes(OBFPowerRank, OBFSSRank, size=count), col = "green") +
  geom_point(data = df4_summary, aes(PropPowerRank, PropSSRank, size=count), col = "yellow")

#Power vs duration
df1_summary = ddply(
  outcomeDF,
  .(`No IA Power`, `No IA Duration`),
  summarize,
  count=length(`No IA Duration`)
)


df2_summary = ddply(
  outcomeDF,
  .(`Wieand Power`, `Wieand Duration`),
  summarize,
  count=length(`Wieand Duration`)
)

df3_summary = ddply(
  outcomeDF,
  .(`OBF Power`, `OBF Duration`),
  summarize,
  count=length(`OBF Duration`)
)

df4_summary = ddply(
  outcomeDF,
  .(`Prop Power`, `Prop Duration`),
  summarize,
  count=length(`Prop Duration`)
)

ggplot(df1_summary, aes(`No IA Power`, `No IA Duration`, size=count)) +
  geom_point(col = "red") + 
  geom_point(data = df2_summary, aes(`Wieand Power`, `Wieand Duration`, size=count), col = "blue") +
  geom_point(data = df3_summary, aes(`OBF Power`, `OBF Duration`, size=count), col = "green") +
  geom_point(data = df4_summary, aes(`Prop Power`, `Prop Duration`, size=count), col = "yellow")   

