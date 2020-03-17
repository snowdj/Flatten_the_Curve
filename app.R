

library(shiny)

library(tidyverse)
library(deSolve)


## Configuration
options(theme_set(theme_minimal(base_size = 15)))


# Function to compute the derivative of the ODE system
# -----------------------------------------------------------
#  t - time
#  y - current state vector of the ODE at time t
#  parms - Parameter vector used by the ODE system

sir <- function(t, y, parms, 
                social_dist_period, 
                reduction) {
    
    beta0 <- parms[1]
    gamma <- parms[2]
    
    # Reduce contact rate 
    beta_t <- if_else(t <= social_dist_period[1], 
                      beta0,
                      if_else(t <= social_dist_period[2], 
                              beta0 * reduction[1],
                              beta0 * reduction[2]
                      )
    )
    
    S <- y[1]
    I <- y[2]
    
    return(list(c(S = -beta_t * S * I, 
                  I =  beta_t * S * I - gamma * I)))
}


## we assume that many globals exist...!
solve_ode <- function(sdp, red, typ) {
    ode_solution <- lsoda(y = c(N - n_init, n_init), 
                          times = times, 
                          func  = sir, 
                          parms = c(beta, gamma), 
                          social_dist_period = sdp,
                          reduction = red) %>%
        as.data.frame() %>%
        setNames(c("t", "S", "I")) %>%
        mutate(beta = beta, 
               gama = gamma,
               R0 = N * beta / gamma, 
               s  = S / N, 
               i  = I / N, 
               type = typ)
    
    daily <- ode_solution %>%
        filter(t %in% seq(0, max_time, by = 1)) %>%
        mutate(C = if_else(row_number() == 1, 0, lag(S, k = 1) - S), 
               c = C / N)
    
    daily
}


## CONSTANTS

# Population size 
N <- 1e6 

# Rate at which person stays in the infectious compartment (disease specific and tracing specific)
gamma <- 1/5 

# Infectious contact rate - beta = R0/N*gamma and when R0  ~2.25 then  2.25/N*gamma
beta <- 4.5e-07 # CAUTION, depends on R0 which is defined below...

# R0 for the beta and gamma values
R0 <- beta*N/gamma

# Initial number of infected people
n_init <- 10

# Grid where to evaluate
max_time <- 150
times <- seq(0, max_time, by = 0.01)


## run only once
if (!exists("ode_solution_daily")) {
    ode_solution_daily <- solve_ode(
        sdp = c(0, max_time),     # social_dist_period
        red = c(0.9999, 0.9999),  # reduction
        typ = "without_intervention"
    )
}

# add results with intervention and plot
run <- function(sdp, red) {
    
    # solve with interventions
    ode_solution2_daily <- solve_ode(
        sdp = sdp,
        red = red,
        typ = "with_intervention"
    )
    
    # Combine the two solutions into one dataset
    ode_df <- rbind(ode_solution_daily, ode_solution2_daily)
    
    # Plot
    pp <- ggplot(ode_df, 
                 aes(x = t, 
                     y = 0, 
                     xend = t, 
                     yend = c, 
                     color = type)) + 
        geom_segment(alpha = 0.7) + 
        geom_line(aes(x = t, y = c)) + 
        labs(
            x = "Days", 
            y = NULL, 
            subtitle = "Daily new cases in % of the population", 
            caption  = "sdp: Start of social distance period\n
                  red: Reduction of 'contacts per day' during period") +
        scale_y_continuous(labels = scales::percent, 
                           limits = c(0, 0.05)) +
        theme(# axis.text.x = element_blank(),
            # axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
        scale_color_brewer(name = "Epidemic model:", type = "qual", palette = 1) +
        # sdp text
        geom_vline(xintercept = sdp, lty = 2, color = "darkgray") +
        geom_text(aes(x = sdp[1], 
                      y = 0.05, 
                      label = "sdp 1"),
                  vjust = 0,
                  color = "blue2") +
        geom_text(aes(x = sdp[2], 
                      y = 0.05, 
                      label = "sdp 2"),
                  vjust = 0,
                  color = "blue2") +
        # red text and arrows
        geom_text(aes(x = sdp[1] + (sdp[2] - sdp[1])/2, 
                      y = 0.048, 
                      label = "red 1"),
                  vjust = 0,
                  color = "brown4") +   
        geom_text(aes(x = sdp[2] + 20, 
                      y = 0.048, 
                      label = "red 2"),
                  vjust = 0,
                  color = "brown4") +   
        geom_segment(aes(
            x = sdp[1],
            y = 0.0475, 
            xend = sdp[2] * 0.995,
            yend = 0.0475
        ),
        size = 0.5, 
        color = "brown4",
        arrow = arrow(length = unit(2, "mm"))) +
        geom_segment(aes(
            x = sdp[2]*1.01,
            y = 0.0475, 
            xend = max_time,
            yend = 0.0475
        ),
        size = 0.5, 
        color = "brown4",
        arrow = arrow(length = unit(2, "mm"))) +
        # capacity
        geom_hline(yintercept = 0.02, lty = 3, color = "steelblue", size = 0.5) +
        geom_text(aes(x = 120, 
                      y = 0.021, 
                      label = "Healthcare system capacity"),
                  color ="steelblue", 
                  size = 3.5)
    print(pp)
}


ui <- fluidPage(

    titlePanel("Flatten the Curve"),
    
    sidebarLayout(
        sidebarPanel(
            # p("Change the reduction"),
            sliderInput("sdp.one",
                        "Start of first 'social distance period' (sdp 1)",
                        min   =  0,
                        max   = 50,
                        value = 30, 
                        step  = 10), 
            sliderInput("red.one",
                        "Reduction during first period (red 1)",
                        min = 0.4,
                        max = 0.9,
                        value = 0.6, 
                        step = 0.1), 
            br(), br(),
            sliderInput("sdp.two",
                        "Start of second 'social distance period' (sdp 2)",
                        min   =  50,
                        max   = 100,
                        value =  60, 
                        step  =  10), 
            sliderInput("red.two",
                        "Reduction during second period (red 2)",
                        min = 0.4,
                        max = 0.9,
                        value = 0.6, 
                        step = 0.1), 
            br(),
            hr(),            
            tags$small("Code on Github:"),
            br(),
            tags$a(href="https://github.com/tinu-schneider/Flatten_the_Curve", 
                   "github.com/tinu-schneider/Flatten_the_Curve"),
            br(), br(),
            tags$div("This work is licensed under a "),
            tags$a(href="http://creativecommons.org/licenses/by-sa/4.0/)", 
                   "Creative Commons Attribution-ShareAlike 4.0 International License"), 
            br(),
            img(src="license.png"),
            br()
        ),

        mainPanel(
           plotOutput("chart", height = "525px"), 
           br(),
           hr(),
           h4("Initial code, mathematical model and idea:"),
           tags$a(href="https://staff.math.su.se/hoehle/blog/2020/03/16/flatteningthecurve.html", 
                  "Michael Höhle, 'Flatten the COVID-19 curve'") 
        )
    )
)


server <- function(input, output) {

    output$chart <- renderPlot({
        run(sdp = c(input$sdp.one, input$sdp.two), 
            red = c(input$red.one, input$red.two))
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
