library(shiny)
library(ggplot2)
library(dplyr)
library(broom)

# ---------------------------
# 1. Load NHANES RData
# ---------------------------
load("/Users/zhuokeyue/Downloads/nhanes.RData")  # loads nhanes_dpq

# ---------------------------
# 2. PHQ-9 recoding function
# ---------------------------
dep2score <- function(x) {
  case_match(as.integer(x),
             1 ~ 0,
             2 ~ 1,
             3 ~ 2,
             4 ~ 3,
             .default = NA_integer_)
}

# ---------------------------
# 3. Clean PHQ-9 scoring
# ---------------------------
nhanes_dpq <- nhanes_dpq |>
  mutate(across(starts_with("DPQ"), dep2score)) |>
  rowwise() |>
  mutate(
    DepScore   = sum(c_across(DPQ010:DPQ090), na.rm = TRUE),
    Depressed  = DepScore >= 10
  ) |>
  ungroup()

# ---------------------------
# 4. Convert variables
# ---------------------------
nhanes_dpq <- nhanes_dpq %>%
  mutate(
    RIAGENDR = factor(RIAGENDR, levels = c(1,2), labels = c("Male","Female")),
    DMDEDUC2 = factor(DMDEDUC2)
  )

# ---------------------------
# 5. Covariates
# ---------------------------
covariates_list <- c("RIDAGEYR", "RIAGENDR", "DMDEDUC2")

# ---------------------------
# 6. UI
# ---------------------------
ui <- fluidPage(
  titlePanel("NHANES PHQ-9 Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("gender", "Select gender:", 
                  choices = levels(nhanes_dpq$RIAGENDR)),
      
      checkboxGroupInput("covariates", "Adjust for:", 
                         choices = covariates_list),
      
      radioButtons("model_type", "Model type:", 
                   choices = c("Unadjusted", "Adjusted"))
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Raw Differences",
                 plotOutput("dist_plot"),
                 tableOutput("summary_table")),
        
        tabPanel("Adjusted Differences",
                 plotOutput("coef_plot"),
                 tableOutput("model_table"))
      )
    )
  )
)

# ---------------------------
# 7. SERVER
# ---------------------------
server <- function(input, output, session) {
  
  # Filtered data based on gender
  filtered_data <- reactive({
    nhanes_dpq %>% filter(RIAGENDR == input$gender)
  })
  
  # ----------------------------------------
  # Raw Distribution Plot  (PHQ-9 by gender)
  # ----------------------------------------
  output$dist_plot <- renderPlot({
    ggplot(nhanes_dpq, aes(x = RIAGENDR, y = DepScore)) +
      geom_boxplot(fill = "skyblue") +
      theme_minimal() +
      labs(title = "PHQ-9 Score Distribution by Gender",
           x = "Gender", y = "PHQ-9 Total Score")
  })
  
  # ----------------------------------------
  # Summary statistics
  # ----------------------------------------
  output$summary_table <- renderTable({
    nhanes_dpq %>%
      group_by(RIAGENDR) %>%
      summarise(
        Mean = mean(DepScore, na.rm = TRUE),
        SD = sd(DepScore, na.rm = TRUE),
        Median = median(DepScore, na.rm = TRUE),
        n = n()
      )
  })
  
  # ----------------------------------------
  # Regression model
  # ----------------------------------------
  model <- reactive({
    
    formula_text <- 
      if (input$model_type == "Unadjusted" || length(input$covariates) == 0) {
        "DepScore ~ RIAGENDR"
      } else {
        paste("DepScore ~ RIAGENDR +", paste(input$covariates, collapse=" + "))
      }
    
    lm(as.formula(formula_text), data = nhanes_dpq)
  })
  
  # Model summary table
  output$model_table <- renderTable({
    tidy(model(), conf.int = TRUE)
  })
  
  # Coefficient Plot
  output$coef_plot <- renderPlot({
    tidy(model(), conf.int = TRUE) %>%
      ggplot(aes(x = term, y = estimate)) +
      geom_point(size = 3, color = "darkred") +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                    width = 0.2, color = "darkred") +
      theme_minimal() +
      labs(title = "Regression Coefficients",
           x = "", y = "Estimate")
  })
}

# ---------------------------
# 8. Run App
# ---------------------------
shinyApp(ui, server)
