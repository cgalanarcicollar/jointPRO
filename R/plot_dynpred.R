plot_dynpred <- function(dynpred_obj = list(), id_data = list(), m , future_time) {
 
    plot_list <- list()
  
  
    n.meas <- nrow(id_data$longitudinal)
      
    # Loop through each measurement to generate the plots
    for (i in seq_len(n.meas)) {
      plot.new()
      par(mar = c(5, 4, 4, 4) + 0.1)
      
      # First subplot: Longitudinal data
      plot.window(xlim = c(0, id_data$longitudinal$time[n.meas] + future_time), ylim = c(0, m))
      points(id_data$longitudinal$time[1:i], id_data$longitudinal$y[1:i], pch = 8)
      if(i!=1){
      points(id_data$longitudinal$time[1:i], dynpred_obj$Long$y_mean[i, 1:i], type = "l", col = "black")
      }
      axis(1)
      axis(2, col.axis = "black")
      box()
      
      # Second subplot: Survival probabilities
      plot.window(xlim = c(0, id_data$longitudinal$time[n.meas] + 3), ylim = c(0, 1))
      lines(seq(id_data$longitudinal$time[i], id_data$longitudinal$time[i] + 3, length.out = 10), dynpred_obj$Surv$Surv_cond_mean[, i])
      lines(seq(id_data$longitudinal$time[i], id_data$longitudinal$time[i] + 3, length.out = 10), dynpred_obj$Surv$Surv_cond_upper[, i], col = "red", lty = 2)
      lines(seq(id_data$longitudinal$time[i], id_data$longitudinal$time[i] + 3, length.out = 10), dynpred_obj$Surv$Surv_cond_lower[, i], col = "red", lty = 2)
      abline(v = id_data$longitudinal$time[i], lty = 2)
      axis(4)
      
      # Titles and labels
      mtext("Survival probability", side = 4, las = 3, line = 2, col = "black")
      mtext("y", side = 2, las = 3, line = 2.5, col = "black")
      mtext("time", side = 1, line = 2, col = "black")
      
      # Record the current plot
      plot_list[[i]] <- recordPlot()
      
    }
    
    # Return the list of recorded plots
    return(plot_list)
  }
  