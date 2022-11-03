



visualisingMCMCOutput <- function(data_vec){
  len_dat <- length(data_vec)/2
  len_dat_1 <- length(outputvec[,1]$TPPval)
  h <- hist(data_vec[,1]$TPPval[2:len_dat_1], xlim=c(0,1), 
            main=paste0("Histogram, t = ", data_vec[,1]$TPPval[1]), xlab = "Proportion of posterior < TPP", freq=F)
  densityMat <- matrix(NA, nrow = len_dat, ncol = length(h$density))
  if (len_dat%%2==0){
    par(mfrow=c(2,len_dat/2))
  } else {
    par(mfrow=c(2,(len_dat+1)/2))
  }
  densityVec <- rep(NA, length = length(len_dat))
  for (i in 1:len_dat){
     h <- hist(data_vec[,i]$TPPval[2:len_dat_1], xlim=c(0,1), 
                                main=paste0("Histogram, t = ", data_vec[,i]$TPPval[1]), xlab = "Proportion of posterior < TPP", freq=F)
     densityVec[i] <- max(h$density)
  }
  for (i in 1:len_dat){
    h <- hist(data_vec[,i]$TPPval[2:len_dat_1], xlim=c(0,1), ylim=c(0, max(densityVec)),
         main=paste0("Histogram, t = ", data_vec[,i]$TPPval[1]), xlab = "Proportion of posterior < TPP", freq=F)
    #print(mean(data_vec[,i]$TPPval[2:len_dat_1]))
  }
  
  # if (len_dat%%2==0){
  #   par(mfrow=c(2,len_dat/2))
  # } else {
  #   par(mfrow=c(2,(len_dat+1)/2))
  # }
  # for (j in 1:len_dat){
  #   plot(data_vec[,j]$HRval[2:len_dat_1], data_vec[,j]$TPPval[2:len_dat_1], ylab = "Proportion of posterior < TPP",
  #        xlab = "Sampled HR", xlim = c(0,1), main = paste0("Plot, t = ", data_vec[,j]$TPPval[1]))
  # }
}



visualisingMCMCOutput(`20october1`)

visualisingMCMCOutput(`18october`)

visualisingMCMCOutput(`19october`)

visualisingMCMCOutput(`20october`)


visualisingMCMCOutput(outputvec)


