#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
library(shiny)
library(ggplot2)
library(dplyr)
library(broom)
library(survey)  # for proper weighting / design-based analyses

# Example: load data
# phq <- haven::read_xpt("P_DPQ.xpt")
# demo <- haven::read_xpt("DEMO_.xpt")  # choose relevant cycle
# df <- phq %>% left_join(demo, by = "SEQN")

# compute total PHQ-9 score
df <- df %>%
  mutate(across(starts_with("DPQ0"), ~ ifelse(. %in% c(7, 9), NA, .))) %>%
  rowwise() %>%
  mutate(PHQ9_total = sum(c_across(starts_with("DPQ0")), na.rm = FALSE)) %>%
  ungroup()

# define survey design (example — adapt to your cycle / weight variable)
nhanes_design <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC2YR,
  data = df,
  nest = TRUE
)

ui <- fluidPage(
  titlePanel("NHANES PHQ-9: Depression by Race"),
  sidebarLayout(
    sidebarPanel(
      selectInput("race_var", "Race/ethnicity:", choices = unique(df$RIDRETH1)),
      checkboxGroupInput("covariates", "Adjust for:", 
                         choices = c("age", "gender", "income_poverty_ratio")),
      radioButtons("model_type", "Model type:", choices = c("Unadjusted", "Adjusted"))
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

server <- function(input, output) {
  filtered <- reactive({
    df %>% filter(RIDRETH1 == input$race_var)
  })
  
  output$dist_plot <- renderPlot({
    ggplot(filtered(), aes(x = RIDRETH1, y = PHQ9_total)) +
      geom_boxplot() +
      theme_minimal()
  })
  
  output$summary_table <- renderTable({
    filtered() %>% summarise(
      Mean = mean(PHQ9_total, na.rm = TRUE),
      SD = sd(PHQ9_total, na.rm = TRUE),
      Median = median(PHQ9_total, na.rm = TRUE),
      n = n()
    )
  })
  
  model <- reactive({
    f <- if (input$model_type == "Unadjusted") {
      as.formula("PHQ9_total ~ RIDRETH1")
    } else {
      as.formula(paste("PHQ9_total ~ RIDRETH1 +",
                       paste(input$covariates, collapse = " + ")))
    }
    svyglm(f, design = nhanes_design)
  })
  
  output$model_table <- renderTable({
    broom::tidy(model(), conf.int = TRUE)
  })
  
  output$coef_plot <- renderPlot({
    broom::tidy(model(), conf.int = TRUE) %>%
      ggplot(aes(x = term, y = estimate)) +
      geom_point() +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
      theme_minimal()
  })
}

shinyApp(ui, server)

