# Requires the deSolve package. Install it once with install.packages("deSolve").
library(deSolve)

eutrophication_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Supporting equations
    sigma_t = (1 - epsilon * cos(2 * pi * t / 365)) / (1 + epsilon)
    sigma_tz = (1 - epsilon * cos(2 * pi * t / 365 - 0.5)) / (1 + epsilon)
    PO4_hypo = 11 + 3 * sin(2 * pi * (t / 365 + 0.3))
    PO4_exog = 0.047 + 0.02 * sin(2 * pi * (t / 365 + 0.12))
    DET_exog = 0.2 + 0.12 * sin(2 * pi * (t / 365 + 0.16))
    outflows = 0.0028 + 0.0014 * sin(2 * pi * (t / 365 + 0.12))
    
    # ODEs
    dPO4_dt = -(PO4 / (e + PO4)) * a * sigma_t * PHYT * P / C_phyt +
      (beta * lambda * ((PHYT * P / C_phyt)^2 + omega * DET^2)) / 
      (mu^2 + (PHYT * P / C_phyt)^2 + omega * DET^2) * sigma_tz * ZOOP * P / C_zoop +
      gamma * d * sigma_tz * (ZOOP^3) / (pred^2 + ZOOP^2) * P / C_zoop +
      phi * sigma_t * DET + k * (1 - sigma_t) * (PO4_hypo - PO4) + PO4_exog - outflows * PO4
    
    dPHYT_dt = (PO4 / (e + PO4)) * a * sigma_t * PHYT - r * sigma_t * PHYT -
      (lambda * (PHYT * P / C_phyt)^2) / (mu^2 + (PHYT * P / C_phyt)^2 + omega * DET^2) * 
      sigma_tz * ZOOP - s * PHYT - outflows * PHYT
    
    dZOOP_dt = (alpha * lambda * ((PHYT * P / C_phyt)^2 + omega * DET^2)) /
      (mu^2 + (PHYT * P / C_phyt)^2 + omega * DET^2) * sigma_tz * ZOOP -
      d * sigma_tz * (ZOOP^3) / (pred^2 + ZOOP^2) - outflows * ZOOP
    
    dDET_dt = r * sigma_t * PHYT * P / C_phyt + 
      (((1 - alpha - beta) * (PHYT * P / C_phyt)^2 - (alpha + beta) * omega * DET^2) * 
         lambda) / (mu^2 + (PHYT * P / C_phyt)^2 + omega * DET^2) * sigma_tz * ZOOP * P / C_zoop -
      phi * sigma_t * DET - psi * DET + DET_exog - outflows * DET
    
    return(list(c(dPO4_dt, dPHYT_dt, dZOOP_dt, dDET_dt)))
  })
}


state <- c(PO4 = 1, PHYT = 0.1, ZOOP = 0.1, DET = 0.1)

parameters <- list(epsilon = 0.9, a = 1, d = 0.12, pred = 10, e = 5, k = 0.02, r = 0.04, s = 0.01,
                   alpha = 0.2, beta = 0.3, gamma = 0.3, lambda = 0.4, mu = 3, phi = 0.07, 
                   psi = 0.03, omega = 0.5, P = 1, C_phyt = 67, C_zoop = 35)

times <- seq(1, 365*5, by = 1)  # From day 1 to day 1825

library(deSolve)

# Initial values for PO4, PHYT, ZOOP, and DET
state <- c(PO4 = 1, PHYT = 0.1, ZOOP = 0.1, DET = 0.1)

# Solve the system
results <- ode(y = state, times = times, func = eutrophication_model, parms = parameters)

# Convert results to a data frame
results_df <- as.data.frame(results)
# 
# # Plotting
# par(mfrow = c(2, 2))  # Set up the plotting area to display 4 plots
# 
# plot(results_df$time, results_df$PO4, type = 'l', main = "PO4 over Time", xlab = "Time (days)", ylab = "PO4")
# plot(results_df$time, results_df$PHYT, type = 'l', main = "PHYT over Time", xlab = "Time (days)", ylab = "PHYT")
# plot(results_df$time, results_df$ZOOP, type = 'l', main = "ZOOP over Time", xlab = "Time (days)", ylab = "ZOOP")
# plot(results_df$time, results_df$DET, type = 'l', main = "DET over Time", xlab = "Time (days)", ylab = "DET")

library(ggplot2)
library(tidyr)  # For pivot_longer

#results_long <- pivot_longer(results_df, -time, names_to = "Variable", values_to = "Value")
# Replace pivot_longer with base R reshape
results_long <- data.frame(
  time = rep(results_df$time, times = 4),
  Variable = rep(c("PO4", "PHYT", "ZOOP", "DET"), each = length(results_df$time)),
  Value = c(results_df$PO4, results_df$PHYT, results_df$ZOOP, results_df$DET)
)

# ggplot(results_long, aes(x = time, y = Value, color = Variable)) +
#   geom_line() +
#   theme_minimal() +
#   labs(title = "Eutrophication Model Results over Time", x = "Time (days)", y = "Concentration") +
#   facet_wrap(~ Variable, scales = 'free_y')

# The calibration block below reads Calibration_data.csv from the working directory.
# Place the course calibration file there before running this section.
calibration_data <- read.csv("Calibration_data.csv")

library(tidyr)

#calibration_long <- pivot_longer(calibration_data, 
#                                 -time, 
#                                 names_to = "Variable", 
#                                 values_to = "Value")

# Replace pivot_longer with base R reshape
calibration_long <- data.frame(
  time = rep(calibration_data$time, times = ncol(calibration_data) - 1),
  Variable = rep(names(calibration_data)[-1], each = nrow(calibration_data)),
  Value = unlist(calibration_data[, -1])
)


results_long$Type <- "Simulated"
calibration_long$Type <- "Observed"
combined_data <- rbind(results_long, calibration_long)


library(ggplot2)

# Assuming combined_data contains both simulated and observed data in long format
# and it has a column 'Type' to differentiate between "Simulated" and "Observed"

ggplot(combined_data, aes(x = time, y = Value, color = Variable)) + 
  geom_line(data = filter(combined_data, Type == "Simulated"), aes(linetype = Type), size = 1) +  # Draw lines for Simulated data
  geom_point(data = filter(combined_data, Type == "Observed"), aes(shape = Type), size = 2) +   # Plot points for Observed data
  theme_minimal() +
  labs(title = "Eutrophication Model Results with Calibration Data", 
       x = "Time (days)", 
       y = "Concentration") +
  scale_shape_manual(values = c("Observed" = 19)) +  # Customize shapes for observed data
  scale_linetype_manual(values = c("Simulated" = "solid")) +  # Ensure simulated data has a solid line
  facet_wrap(~ Variable, scales = 'free_y')


library(ggplot2)
library(dplyr)

# Calculate the start time for the last 5th year, assuming the model runs for 5 years (1825 days)
start_last_year <- max(combined_data$time) - 365+1

# Filter the combined_data to include only data from the last year
last_year_data <- combined_data %>%
  filter(time > start_last_year)

# Plot the filtered data
ggplot(last_year_data, aes(x = time, y = Value, color = Variable)) + 
  geom_line(data = filter(last_year_data, Type == "Simulated"), size = 1) +  # Draw lines for Simulated data
  geom_point(data = filter(last_year_data, Type == "Observed"), size = 2, shape = 19) +   # Plot points for Observed data
  theme_minimal() +
  labs(title = "Eutrophication Model Results with Calibration Data (Last Year)", 
       x = "Time (days)", 
       y = "Concentration") +
  facet_wrap(~ Variable, scales = 'free_y')
############

#assessing goodness of fit
library(dplyr)
library(tidyr)

combined_last_year <- combined_data %>%
  filter(time > max(time) - 365)

# Separate observed and simulated data
observed_data <- combined_last_year %>% filter(Type == "Observed")
simulated_data <- combined_last_year %>% filter(Type == "Simulated")
# Assuming combined_last_year contains both observed and simulated data for the last year

# Join observed and simulated data based on time and Variable to ensure matching pairs
data_joined <- observed_data %>%
  inner_join(simulated_data, by = c("time", "Variable"), suffix = c("_obs", "_sim"))

# Function to calculate RMSE
calculate_rmse <- function(observed, simulated) {
  sqrt(mean((observed - simulated) ^ 2, na.rm = TRUE))
}

# Calculate RMSE for each Variable
rmse_results <- data_joined %>%
  group_by(Variable) %>%
  summarise(RMSE = calculate_rmse(Value_obs, Value_sim))

# Print the RMSE results
print(rmse_results)
rmse_table<-data.frame(t(rmse_results))
rmse_table$MC<-0
rmse_table[2,5]<-1

##########

