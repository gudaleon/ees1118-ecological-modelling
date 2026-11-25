library(deSolve)

# Lotka-Volterra with Holling Type II: G(R) = Gmax * R / (R + K)
lv_hollingII <- function(t, state, parms) {
  with(as.list(c(state, parms)), {
    
    # Holling Type II functional response using slide notation
    G <- Gmax * x / (x + Kfun)  # G(x)
    
    # System of equations
    dx_dt <- b * x * (1 - x / K) - G * y
    dy_dt <- e * G * y - d * y
    
    list(c(dx_dt, dy_dt))
  })
}

# Initial values
state <- c(x = 1, y = 0.5)

# Parameters (matching slide notation)
parms <- c(
  b = 1.0,       # intrinsic growth rate of producer
  K = 5.0,       # carrying capacity of producer
  Gmax = 1.5,    # maximum ingestion rate
  Kfun = 0.5,    # half-saturation constant
  e = 0.4,       # efficiency of predator assimilation
  d = 0.2        # predator loss rate
)

# Time range
times <- seq(0, 150, by = 0.1)

# Solve
out <- ode(y = state, times = times, func = lv_hollingII, parms = parms)

# Plot
matplot(out[, "time"], out[, 2:3], type = "l", lty = 1, col = c("forestgreen", "darkred"),
        xlab = "Time", ylab = "Biomass (mg C/L)")
legend("topright", legend = c("Producer (x)", "Grazer (y)"),
       col = c("forestgreen", "darkred"), lty = 1)
