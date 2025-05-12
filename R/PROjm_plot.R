plot_longitudinal_survival <- function(data_long, outcome_col = "y", time_col = "time", data_surv, Time_col = "Time", 
                                       status_col = "status", id_col = "id", y_max, ids = NULL, max_time = 6.5) {
  
  # Subset data for selected IDs (if specified)
  if (!is.null(ids)) {
    data_long <- data_long[data_long[[id_col]] %in% ids, ]
    data_surv <- data_surv[data_surv[[id_col]] %in% ids, ]
  }
  
  library(dplyr)
  
  # Join the datasets
  data_joined <- data_long %>%
    left_join(data_surv, by = "id")
  
  # Update the status column
  data_adjusted <- data_joined %>%
    group_by(id) %>%                               # Group by subject ID
    mutate(
      status = ifelse(status == 1 & row_number() != n(), 0, status)
    ) %>%
    ungroup()
  
  # Create a working copy of the data with relevant columns
  new_data <- data_adjusted[, c(id_col, outcome_col, time_col, Time_col, status_col)]
  colnames(new_data) <- c("id", "y", "time", "Time", "status")
  
  # Add "Failure" rows for deceased patients
  for (i in seq_len(nrow(new_data))) {
    id <- new_data[["id"]][i]
    Time <- new_data[["Time"]][i]
    y <- new_data[["y"]][i]
    status <- 2  # Mark as "Failure"
    new_row <- data.frame(id = id, y = y, time = Time, Time = Time, status = status, 
                          stringsAsFactors = FALSE)
    if (new_data[i, "status"] == 1) {
      new_data <- rbind(new_data, new_row)
    }
  }
  
  # Sort by ID for continuity
  new_data <- new_data[order(new_data[["id"]]), ]
  
  # Add measurement type
  new_data$meas <- as.factor(ifelse(new_data$status == 2, "Failure", 
                                    ifelse(round(new_data[["time"]], 2) == round(new_data[["Time"]], 2), 
                                           "Censored", "Measurement")))
  
  # Prepare dashed lines for "Failure" measurements
  fails <- which(new_data$meas == "Failure")
  dash <- new_data[c(fails - 1, fails), ]
  dash <- dash[order(dash[["id"]]), ]
  
  # Create a mapping of possible shapes and labels
  shape_mapping <- list(
    "Failure" = list(label = "Deceased", value = 4),
    "Measurement" = list(label = "Measurement", value = 16),
    "Censored" = list(label = "Censored", value = 1)
  )
  
  # Determine which shapes are actually present in the data
  present_meas <- levels(droplevels(new_data$meas))
  present_shapes <- shape_mapping[present_meas]
  
  # Extract labels and values for present shapes
  shape_labels <- sapply(present_shapes, function(x) x$label)
  shape_values <- sapply(present_shapes, function(x) x$value)
  
  # Generate the plot
  plot <- ggplot() +
    # Dashed lines for "Failure" data
    geom_point(data = dash, aes(x = time, y = y, group = as.factor(id), 
                                colour = as.factor(id), shape = meas), size = 2) +
    geom_line(data = dash, aes(x = time, y = y, group = as.factor(id), 
                               colour = as.factor(id)), linetype = "dashed") +
    # Solid lines for other data
    geom_point(data = new_data %>% filter(meas != "Failure"), aes(x = time, y = y, 
                                                                  group = as.factor(id), 
                                                                  colour = as.factor(id), shape = meas), size = 2) +
    geom_line(data = new_data %>% filter(meas != "Failure"), aes(x = time, y = y, 
                                                                 group = as.factor(id), 
                                                                 colour = as.factor(id)), linetype = "solid") +
    # Dynamic shape scale based on what's present
    scale_shape_manual(labels = shape_labels, values = shape_values) +
    scale_y_continuous(limits = c(0, y_max), breaks = seq(0, y_max, 1)) +
    xlim(0, max_time) +
    xlab("Time") +
    ylab(outcome_col) +
    labs(color = "Patient ID", shape = "Data Type") +
    theme_classic()
  
  return(plot)
}
