



visualisingMCMCOutput <- function(data_vec){
  len_dat <- length(data_vec)/2
  if (len_dat%%2==0){
    par(mfrow=c(2,len_dat/2))
  } else {
    par(mfrow=c(2,(len_dat+1)/2))
  }
  for (i in 1:len_dat){
    hist(data_vec[,i]$TPPval[2:31], ylim=c(0, 20), xlim=c(0,1), 
         main=paste0("Histogram, t = ", data_vec[,i]$TPPval[1]), xlab = "Proportion of posterior < TPP")
    print(mean(data_vec[,i]$TPPval[2:31]))
  }
  
  if (len_dat%%2==0){
    par(mfrow=c(2,len_dat/2))
  } else {
    par(mfrow=c(2,(len_dat+1)/2))
  }
  for (j in 1:len_dat){
    plot(data_vec[,j]$HRval[2:31], data_vec[,j]$TPPval[2:31], ylim=c(0,1), ylab = "Proportion of posterior < TPP",
         xlab = "Sampled HR", xlim = c(0,1), main = paste0("Plot, t = ", data_vec[,j]$TPPval[1]))
  }
}



visualisingMCMCOutput(`20october1`)

visualisingMCMCOutput(`18october`)

visualisingMCMCOutput(`19october`)

visualisingMCMCOutput(`20october`)



