library(gridExtra)
library(ggplot2)
DelayTimeSeq <- seq(10, 130, by = 5)
delayProseq <- seq(11, 131, by = 5)
powerseq <- seq(12, 132, by = 5)
durationseq <- seq(13, 133, by = 5)
SSSeq <- seq(14, 134, by = 5)


#Power vs Duration
par(mfrow  = c(4,3))
for (k in 1:12){
  plot(outcomeDF[k,powerseq[1]], outcomeDF[k,durationseq[1]], xlab = "Power", 
       xlim = c(0,1), ylab = "Duration", ylim = c(10, 50), pch = 19, main = paste0("HR", outcomeDF$HR1Vec[k], 
                                  " for ", outcomeDF$T1Vec[k], " months, then HR", outcomeDF$HR2Vec[k], ": R = ", outcomeDF$recTimeVec[k]))
  for (j in 2:25){
    points(outcomeDF[k,powerseq[j]], outcomeDF[k,durationseq[j]], pch = 19)
  }
}

#Power vs SS
par(mfrow  = c(4,3))
for (k in 1:12){
  plot(outcomeDF[k,powerseq[1]], outcomeDF[k,SSSeq[1]], xlab = "Power", 
       xlim = c(0,1), ylab = "Sample Size", ylim = c(450, 700), pch = 19, main = paste0("HR", outcomeDF$HR1Vec[k], 
                                   " for ", outcomeDF$T1Vec[k], " months, then HR", outcomeDF$HR2Vec[k], ": R = ", outcomeDF$recTimeVec[k]))
  for (j in 2:25){
    points(outcomeDF[k,powerseq[j]], outcomeDF[k,SSSeq[j]], pch = 19)
  }
}

#Duration vs SS
par(mfrow  = c(4,3))
for (k in 1:12){
  plot(outcomeDF[k,durationseq[1]], outcomeDF[k,SSSeq[1]], xlab = "Duration", 
       xlim = c(10,50), ylab = "Sample Size", ylim = c(450, 700), pch = 19, main = paste0("HR", outcomeDF$HR1Vec[k], 
                                 " for ", outcomeDF$T1Vec[k], " months, then HR", outcomeDF$HR2Vec[k], ": R = ", outcomeDF$recTimeVec[k]))
  for (j in 2:25){
    points(outcomeDF[k,durationseq[j]], outcomeDF[k,SSSeq[j]], pch = 19)
  }
}

# Plotting using ggplot2
plot_data <- data.frame(
  Power = as.numeric(outcomeDF[k, powerseq]),
  Duration = as.numeric(outcomeDF[k, durationseq]),
  ColorColumn = as.factor(as.numeric(outcomeDF[k, DelayTimeSeq])),
  ShapeColumn = as.factor(as.numeric(outcomeDF[k, delayPropseq])) 
)


# Define custom colors and shapes
custom_colors <- c("red", "blue", "green", "purple", "orange")
custom_shapes <- 0:4

p <- ggplot(plot_data, aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Duration"
  ) +
  theme_minimal() +
  xlim(0, 1) +
  ylim(10, 50) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Power = as.numeric(outcomeDF[k, powerseq[j]]),
    Duration = as.numeric(outcomeDF[k, durationseq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p <- p + geom_point(
    data = plot_data_j,
    aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)
  )
}

# Display the plot
print(p)

