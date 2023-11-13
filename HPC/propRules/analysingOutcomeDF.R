library(gridExtra)
library(ggplot2)
DelayTimeSeq <- seq(10, 130, by = 5)
delayPropseq <- seq(11, 131, by = 5)
powerseq <- seq(12, 132, by = 5)
durationseq <- seq(13, 133, by = 5)
SSSeq <- seq(14, 134, by = 5)


#png(file="PowerVsDuration.png",
  #  width=200, height=280, units = "mm", res = 300)
#Power vs Duration
par(mfrow  = c(4,3))
for (k in 1:12){
  plot(outcomeDF[k,powerseq[1]], outcomeDF[k,durationseq[1]], xlab = "Power", 
       xlim = c(0,1), ylab = "Duration", cex.main = 0.75,  ylim = c(10, 50), pch = 19, main = paste0("HR", outcomeDF$HR1Vec[k], 
                                  " for ", outcomeDF$T1Vec[k], " months, then HR", outcomeDF$HR2Vec[k], ": R = ", outcomeDF$recTimeVec[k]))
  for (j in 2:25){
    points(outcomeDF[k,powerseq[j]], outcomeDF[k,durationseq[j]], pch = 19)
  }
}
#dev.off()

#png(file="PowerVsSS.png",
   # width=200, height=280, units = "mm", res = 300)
#Power vs SS
par(mfrow  = c(4,3))
for (k in 1:12){
  plot(outcomeDF[k,powerseq[1]], outcomeDF[k,SSSeq[1]], xlab = "Power", 
       xlim = c(0,1), ylab = "Sample Size", cex.main = 0.75, ylim = c(450, 700), pch = 19, main = paste0("HR", outcomeDF$HR1Vec[k], 
                                   " for ", outcomeDF$T1Vec[k], " months, then HR", outcomeDF$HR2Vec[k], ": R = ", outcomeDF$recTimeVec[k]))
  for (j in 2:25){
    points(outcomeDF[k,powerseq[j]], outcomeDF[k,SSSeq[j]], pch = 19)
  }
}
#dev.off()

#png(file="DurationVsSS.png",
   # width=200, height=280, units = "mm", res = 300)
#Duration vs SS
par(mfrow  = c(4,3))
for (k in 1:12){
  plot(outcomeDF[k,durationseq[1]], outcomeDF[k,SSSeq[1]], xlab = "Duration", 
       xlim = c(10,50), ylab = "Sample Size", cex.main=0.75, ylim = c(450, 700), pch = 19, main = paste0("HR", outcomeDF$HR1Vec[k], 
                                 " for ", outcomeDF$T1Vec[k], " months, then HR", outcomeDF$HR2Vec[k], ": R = ", outcomeDF$recTimeVec[k]))
  for (j in 2:25){
    points(outcomeDF[k,durationseq[j]], outcomeDF[k,SSSeq[j]], pch = 19)
  }
}
#dev.off()

# Row 8 -------------------------------------------------------------------

k <- 8

#Power vs Duration

# Plotting using ggplot2
plot_data <- data.frame(
  Power = as.numeric(outcomeDF[k, powerseq]),
  Duration = as.numeric(outcomeDF[k, durationseq]),
  SS = as.numeric(outcomeDF[k, durationseq]),
  ColorColumn = as.factor(as.numeric(outcomeDF[k, DelayTimeSeq])),
  ShapeColumn = as.factor(as.numeric(outcomeDF[k, delayPropseq])) 
)

# Define custom colors and shapes
custom_colors <- c("red", "blue", "green", "purple", "orange")
custom_shapes <- 0:4

p1 <- ggplot(plot_data, aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Duration"
  ) +
  theme_minimal() +
  xlim(0, 0.05) +
  ylim(15, 30) +
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
  p1 <- p1 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 1000 months, then HR1: R = 12")
}

# Display the plot
#print(p1)


#Power vs SS

p2 <- ggplot(plot_data, aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(0, 0.05) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Power = as.numeric(outcomeDF[k, powerseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p2 <- p2 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 1000 months, then HR1: R = 12")
}

# Display the plot
#print(p2)


#Duration vs SS

p3 <- ggplot(plot_data, aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Duration",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(10,40) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Duration = as.numeric(outcomeDF[k, durationseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p3 <- p3 + geom_point(
    data = plot_data_j,
    aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 1000 months, then HR1: R = 12")
}

# Display the plot
#print(p3)

# Consolidate legends
p1 <- p1 + guides(color = "none", shape = "none")
p2 <- p2 + guides(color = "none", shape = "none")
p3 <- p3 + guides(color = "none", shape = "none")

# Arrange the plots
png(file="Row8.png", width=280, height=200, units = "mm", res = 300)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()



# Row 9 -------------------------------------------------------------------

k <- 9

#Power vs Duration

# Plotting using ggplot2
plot_data <- data.frame(
  Power = as.numeric(outcomeDF[k, powerseq]),
  Duration = as.numeric(outcomeDF[k, durationseq]),
  SS = as.numeric(outcomeDF[k, durationseq]),
  ColorColumn = as.factor(as.numeric(outcomeDF[k, DelayTimeSeq])),
  ShapeColumn = as.factor(as.numeric(outcomeDF[k, delayPropseq])) 
)

# Define custom colors and shapes
custom_colors <- c("red", "blue", "green", "purple", "orange")
custom_shapes <- 0:4

p1 <- ggplot(plot_data, aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Duration"
  ) +
  theme_minimal() +
  xlim(0, 0.05) +
  ylim(10, 30) +
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
  p1 <- p1 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1.3 for 1000 months, then HR1.3: R = 12")
}

# Display the plot
#print(p1)


#Power vs SS

p2 <- ggplot(plot_data, aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(0, 0.05) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Power = as.numeric(outcomeDF[k, powerseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p2 <- p2 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1.3 for 1000 months, then HR1.3: R = 12")
}

# Display the plot
#print(p2)


#Duration vs SS

p3 <- ggplot(plot_data, aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Duration",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(10,30) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Duration = as.numeric(outcomeDF[k, durationseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p3 <- p3 + geom_point(
    data = plot_data_j,
    aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1.3 for 1000 months, then HR1.3: R = 12")
}

# Display the plot
#print(p3)

# Consolidate legends
p1 <- p1 + guides(color = "none", shape = "none")
p2 <- p2 + guides(color = "none", shape = "none")
p3 <- p3 + guides(color = "none", shape = "none")

# Arrange the plots
png(file="Row9.png", width=280, height=200, units = "mm", res = 300)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()


# Row 11 -------------------------------------------------------------------

k <- 11

#Power vs Duration

# Plotting using ggplot2
plot_data <- data.frame(
  Power = as.numeric(outcomeDF[k, powerseq]),
  Duration = as.numeric(outcomeDF[k, durationseq]),
  SS = as.numeric(outcomeDF[k, durationseq]),
  ColorColumn = as.factor(as.numeric(outcomeDF[k, DelayTimeSeq])),
  ShapeColumn = as.factor(as.numeric(outcomeDF[k, delayPropseq])) 
)

# Define custom colors and shapes
custom_colors <- c("red", "blue", "green", "purple", "orange")
custom_shapes <- 0:4

p1 <- ggplot(plot_data, aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Duration"
  ) +
  theme_minimal() +
  xlim(0.6, 1) +
  ylim(20, 45) +
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
  p1 <- p1 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 6 months, then HR0.62: R = 12")
}

# Display the plot
#print(p1)


#Power vs SS

p2 <- ggplot(plot_data, aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(0.6, 1) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Power = as.numeric(outcomeDF[k, powerseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p2 <- p2 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 6 months, then HR0.62: R = 12")
}

# Display the plot
#print(p2)


#Duration vs SS

p3 <- ggplot(plot_data, aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Duration",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(20,45) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Duration = as.numeric(outcomeDF[k, durationseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p3 <- p3 + geom_point(
    data = plot_data_j,
    aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 6 months, then HR0.62: R = 12")
}

# Display the plot
#print(p3)

# Consolidate legends
p1 <- p1 + guides(color = "none", shape = "none")
p2 <- p2 + guides(color = "none", shape = "none")
p3 <- p3 + guides(color = "none", shape = "none")

# Arrange the plots
png(file="Row11.png", width=280, height=200, units = "mm", res = 300)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()

# Row 12 -------------------------------------------------------------------

k <- 12

#Power vs Duration

# Plotting using ggplot2
plot_data <- data.frame(
  Power = as.numeric(outcomeDF[k, powerseq]),
  Duration = as.numeric(outcomeDF[k, durationseq]),
  SS = as.numeric(outcomeDF[k, durationseq]),
  ColorColumn = as.factor(as.numeric(outcomeDF[k, DelayTimeSeq])),
  ShapeColumn = as.factor(as.numeric(outcomeDF[k, delayPropseq])) 
)

# Define custom colors and shapes
custom_colors <- c("red", "blue", "green", "purple", "orange")
custom_shapes <- 0:4

p1 <- ggplot(plot_data, aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Duration"
  ) +
  theme_minimal() +
  xlim(0.5, 1) +
  ylim(20, 45) +
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
  p1 <- p1 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = Duration, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 3 months, then HR0.628: R = 12")
}

# Display the plot
#print(p1)


#Power vs SS

p2 <- ggplot(plot_data, aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Power",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(0.5, 1) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Power = as.numeric(outcomeDF[k, powerseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p2 <- p2 + geom_point(
    data = plot_data_j,
    aes(x = Power, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 3 months, then HR0.628: R = 12")
}

# Display the plot
#print(p2)


#Duration vs SS

p3 <- ggplot(plot_data, aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)) +
  geom_point() +
  labs(
    x = "Duration",
    y = "Sample size"
  ) +
  theme_minimal() +
  xlim(20,45) +
  ylim(650,700) +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors)

# Adding points for j from 2 to 25
for (j in 2:25) {
  plot_data_j <- data.frame(
    Duration = as.numeric(outcomeDF[k, durationseq[j]]),
    SS = as.numeric(outcomeDF[k, SSSeq[j]]),
    ColorColumn = as.factor(outcomeDF[k, DelayTimeSeq[j]]),
    ShapeColumn = as.factor(outcomeDF[k, delayPropseq[j]])
  )
  p3 <- p3 + geom_point(
    data = plot_data_j,
    aes(x = Duration, y = SS, color = ColorColumn, shape = ShapeColumn)
  ) + ggtitle("HR1 for 3 months, then HR0.628: R = 12")
}

# Display the plot
#print(p3)

# Consolidate legends
p1 <- p1 + guides(color = "none", shape = "none")
p2 <- p2 + guides(color = "none", shape = "none")
p3 <- p3 + guides(color = "none", shape = "none")

# Arrange the plots
png(file="Row12.png", width=280, height=200, units = "mm", res = 300)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()







