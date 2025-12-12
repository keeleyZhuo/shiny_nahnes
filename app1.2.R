
# ------------------------------------------------------------
# Packages
# ------------------------------------------------------------
library(shiny)
library(dplyr)
library(ggplot2)
library(survey)
library(broom)
library(forcats)
library(purrr)
library(janitor)
library(nhanesA)

# ------------------------------------------------------------
# 1. Cycle definitions
# ------------------------------------------------------------
cycle_years <- list(
  "D" = "2005-2006",
  "E" = "2007-2008",
  "F" = "2009-2010",
  "G" = "2011-2012",
  "H" = "2013-2014",
  "I" = "2015-2016",
  "J" = "2017-2020"
)
cycles <- names(cycle_years)

# ------------------------------------------------------------
# 2. Helper: recode DPQ items (NHANES → PHQ-9 scores)
# ------------------------------------------------------------
dep2score <- function(x) {
  case_when(
    x == "Not at all" ~ 0L,
    x == "Several days" ~ 1L,
    x == "More than half the days" ~ 2L,
    x == "Nearly every day" ~ 3L,
    TRUE ~ NA_integer_
  )
}

# ------------------------------------------------------------
# 3. Load and standardize one cycle
# ------------------------------------------------------------
load_cycle <- function(cycle_letter) {
  cat("Loading cycle", cycle_letter, "...
")
  demo <- tryCatch(nhanes(paste0("DEMO_", cycle_letter)), error = function(e) NULL)
  dpq  <- tryCatch(nhanes(paste0("DPQ_",  cycle_letter)), error = function(e) NULL)
  if (is.null(demo) || is.null(dpq)) {
    cat("Cycle", cycle_letter, "not available.
")
    return(NULL)
  }
  merged <- merge(demo, dpq, by = "SEQN", all = FALSE) %>%
    mutate(across(everything(), as.character))
  merged$cycle <- cycle_years[[cycle_letter]]
  merged
}

# ------------------------------------------------------------
# 4. Load all cycles + stack into one dataset
# ------------------------------------------------------------
raw_all <- map_df(cycles, load_cycle)
cat("All cycles successfully combined.
")

# ------------------------------------------------------------
# 5. Clean names
# ------------------------------------------------------------
cleaned <- raw_all %>% janitor::clean_names()

# ------------------------------------------------------------
# 6. Recode PHQ-9 items and compute total score
# ------------------------------------------------------------
phq_vars <- paste0("dpq0", seq(10, 90, 10))  # dpq010 ... dpq090
cleaned <- cleaned %>%
  mutate(across(all_of(phq_vars), dep2score)) %>%
  mutate(phq9_total = rowSums(across(all_of(phq_vars)), na.rm = FALSE))

# ------------------------------------------------------------
# 7. Harmonize demographics, severity, and year
# ------------------------------------------------------------
harmonized <- cleaned %>%
  mutate(
    race_eth = case_when(
      ridreth1 == "Non-Hispanic White" ~ "Non-Hispanic White",
      ridreth1 == "Non-Hispanic Black" ~ "Non-Hispanic Black",
      ridreth1 %in% c("Mexican American","Other Hispanic") ~ "Hispanic",
      TRUE ~ "Other"
    ),
    race_eth = factor(race_eth,
                      levels = c("Non-Hispanic White","Non-Hispanic Black","Hispanic","Other")),
    gender = case_when(
      riagendr == "Male" ~ "Male",
      riagendr == "Female" ~ "Female",
      TRUE ~ NA_character_
    ),
    gender = factor(gender, levels = c("Male","Female")),
    age = suppressWarnings(as.numeric(ridageyr)),
    pir = suppressWarnings(as.numeric(indfmpir)),
    wtmec2yr = suppressWarnings(as.numeric(wtmec2yr)),
    sdmvpsu  = suppressWarnings(as.numeric(sdmvpsu)),
    sdmvstra = suppressWarnings(as.numeric(sdmvstra)),
    year = as.integer(substr(cycle, 1, 4)),
    phq9_severity = cut(
      phq9_total,
      breaks = c(-Inf, 4, 9, 14, 19, Inf),
      labels = c("Minimal (0–4)","Mild (5–9)","Moderate (10–14)",
                 "Moderately severe (15–19)","Severe (20–27)")
    )
  ) %>%
  filter(!is.na(age), age >= 18)

# ------------------------------------------------------------
# 8. Shiny app
# ------------------------------------------------------------
make_formula <- function(response, covars, include_race = TRUE) {
  rhs <- character(0)
  if (include_race) rhs <- c(rhs, "race_eth")
  if (length(covars)) rhs <- c(rhs, covars)
  as.formula(paste(response, "~", paste(rhs, collapse = " + ")))
}

ui <- fluidPage(
  titlePanel("NHANES PHQ-9 Disparities Explorer (2005–2020)"),
  sidebarLayout(
    sidebarPanel(
      h4("Filters and settings"),
      selectInput("years", "Survey cycles",
                  choices = sort(unique(harmonized$cycle)),
                  multiple = TRUE, selected = sort(unique(harmonized$cycle))),
      selectInput("raceRef","Reference race/ethnicity",
                  choices = c("Non-Hispanic White","Non-Hispanic Black","Hispanic","Other"),
                  selected = "Non-Hispanic White"),
      checkboxInput("useWeights","Use survey weights", TRUE),
      checkboxGroupInput("covars","Adjust for",
                         choices = c("Age"="age","Gender"="gender","Income (PIR)"="pir","Year"="year"),
                         selected = c("age","gender","pir","year")),
      helpText("Models are design-based (NHANES weights) when enabled. Coefficients reflect associations, not causal effects.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Distributions",
                 fluidRow(
                   column(6, plotOutput("distPlot")),
                   column(6, plotOutput("severityPlot"))
                 ),
                 br(),
                 h5("Summary by race/ethnicity"),
                 tableOutput("summaryTable")),
        tabPanel("Trends",
                 plotOutput("trendPlot", height = "500px")),
        tabPanel("Models",
                 plotOutput("forestPlot", height = "500px"),
                 br(),
                 h5("Race coefficients (unadjusted vs adjusted)"),
                 tableOutput("modelTable")),
        tabPanel("Items",
                 plotOutput("itemHeat", height = "520px"),
                 br(),
                 h5("Adjusted item-level coefficients with BH q-values"),
                 tableOutput("itemTable"))
      )
    )
  )
)

server <- function(input, output, session) {
  dat_f <- reactive({
    harmonized %>%
      filter(cycle %in% input$years) %>%
      filter(!is.na(phq9_total), !is.na(race_eth)) %>%
      mutate(race_eth = fct_relevel(race_eth, input$raceRef))
  })
  
  des_r <- reactive({
    df <- dat_f()
    if (nrow(df) == 0) return(svydesign(ids = ~1, data = df))
    if (input$useWeights &&
        all(c("wtmec2yr","sdmvstra","sdmvpsu") %in% names(df))) {
      svydesign(ids = ~sdmvpsu, strata = ~sdmvstra,
                weights = ~wtmec2yr, nest = TRUE,
                data = df)
    } else {
      svydesign(ids = ~1, data = df)
    }
  })
  
  # A. Raw differences
  output$distPlot <- renderPlot({
    df <- dat_f(); validate(need(nrow(df) > 0, "No data available"))
    ggplot(df, aes(x = phq9_total, fill = race_eth)) +
      geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.25, bins = 28) +
      geom_density(alpha = 0.35) +
      scale_x_continuous(limits = c(0, 27), breaks = seq(0, 27, 3)) +
      labs(x = "PHQ-9 total (0–27)", y = "Density", fill = "Race/ethnicity",
           title = "PHQ-9 distributions by race/ethnicity") +
      theme_minimal()
  })
  
  output$severityPlot <- renderPlot({
    df <- dat_f(); validate(need(nrow(df) > 0, "No data available"))
    df %>%
      filter(!is.na(phq9_severity)) %>%
      count(race_eth, phq9_severity) %>%
      group_by(race_eth) %>%
      mutate(p = n / sum(n)) %>%
      ungroup() %>%
      ggplot(aes(x = race_eth, y = p, fill = phq9_severity)) +
      geom_col(position = "fill") +
      scale_y_continuous(labels = scales::percent) +
      labs(x = "Race/ethnicity", y = "Percent", fill = "Severity",
           title = "PHQ-9 severity categories by race/ethnicity") +
      theme_minimal()
  })
  
  # B. Trends over time
  output$trendPlot <- renderPlot({
    df <- dat_f(); validate(need(nrow(df) > 0, "No data available"))
    df %>%
      group_by(cycle, race_eth) %>%
      summarise(mean_phq9 = mean(phq9_total, na.rm = TRUE), .groups = "drop") %>%
      ggplot(aes(x = cycle, y = mean_phq9, color = race_eth, group = race_eth)) +
      geom_point() + geom_line() +
      labs(x = "Survey cycle", y = "Mean PHQ-9 total",
           color = "Race/ethnicity",
           title = "Mean PHQ-9 total over survey cycles by race/ethnicity") +
      theme_minimal()
  })
  
  # C. Models: Forest plot and table
  model_results <- reactive({
    df <- dat_f()
    des <- des_r()
    covars <- input$covars
    race_ref <- input$raceRef
    
    # Unadjusted model
    formula_unadj <- make_formula("phq9_total", character(0), include_race = TRUE)
    fit_unadj <- svyglm(formula_unadj, design = des)
    tidy_unadj <- broom::tidy(fit_unadj) %>% filter(grepl("race_eth", term))
    
    # Adjusted model
    formula_adj <- make_formula("phq9_total", covars, include_race = TRUE)
    fit_adj <- svyglm(formula_adj, design = des)
    tidy_adj <- broom::tidy(fit_adj) %>% filter(grepl("race_eth", term))
    
    # Combine
    combined <- tidy_unadj %>%
      select(term, estimate_unadj = estimate, std.error_unadj = std.error, p.value_unadj = p.value) %>%
      left_join(tidy_adj %>%
                  select(term, estimate_adj = estimate, std.error_adj = std.error, p.value_adj = p.value),
                by = "term")
    
    combined %>% mutate(
      term = gsub("race_eth", "", term),
      term = factor(term, levels = levels(harmonized$race_eth))
    )
  })
  
  output$forestPlot <- renderPlot({
    res <- model_results(); validate(need(nrow(res) > 0, "No model results"))
    res_long <- res %>%
      pivot_longer(cols = c(estimate_unadj, estimate_adj), names_to = "model", values_to = "estimate") %>%
      mutate(model = recode(model, estimate_unadj = "Unadjusted", estimate_adj = "Adjusted"))
    
    ggplot(res_long, aes(x = term, y = estimate, color = model)) +
      geom_point(position = position_dodge(width = 0.5), size = 3) +
      geom_errorbar(aes(ymin = estimate - 1.96 * std.error_adj, ymax = estimate + 1.96 * std.error_adj),
                    position = position_dodge(width = 0.5), width = 0.2) +
      labs(x = "Race/ethnicity", y = "Coefficient (PHQ-9 total)", color = "Model",
           title = "Race/ethnicity coefficients from linear models") +
      theme_minimal()
  })
  
  output$modelTable <- renderTable({
    res <- model_results(); validate(need(nrow(res) > 0, "No model results"))
    res %>%
      mutate(
        estimate_unadj = round(estimate_unadj, 3),
        std.error_unadj = round(std.error_unadj, 3),
        p.value_unadj = signif(p.value_unadj, 3),
        estimate_adj = round(estimate_adj, 3),
        std.error_adj = round(std.error_adj, 3),
        p.value_adj = signif(p.value_adj, 3),
        term = as.character(term)
      ) %>%
      rename(
        "Race/Ethnicity" = term,
        "Estimate (Unadjusted)" = estimate_unadj,
        "Std. Error (Unadjusted)" = std.error_unadj,
        "P-value (Unadjusted)" = p.value_unadj,
        "Estimate (Adjusted)" = estimate_adj,
        "Std. Error (Adjusted)" = std.error_adj,
        "P-value (Adjusted)" = p.value_adj
      )
  })
  
  # D. Item-level heatmap and table
  item_results <- reactive({
    df <- dat_f()
    des <- des_r()
    covars <- input$covars
    
    phq_items <- paste0("dpq0", seq(10, 90, 10))
    
    # Fit models for each item
    res_list <- map(phq_items, function(item) {
      formula <- make_formula(item, covars, include_race = TRUE)
      fit <- svyglm(formula, design = des)
      tidy_fit <- broom::tidy(fit) %>% filter(grepl("race_eth", term))
      tidy_fit$item <- item
      tidy_fit
    })
    
    res_df <- bind_rows(res_list) %>%
      mutate(
        term = gsub("race_eth", "", term),
        item = factor(item, levels = phq_items)
      )
    
    # Adjust p-values for multiple testing (BH)
    res_df <- res_df %>% group_by(term) %>%
      mutate(q.value = p.adjust(p.value, method = "BH")) %>% ungroup()
    
    res_df
  })
  
  output$itemHeat <- renderPlot({
    res <- item_results(); validate(need(nrow(res) > 0, "No item results"))
    ggplot(res, aes(x = term, y = term, fill = estimate)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(x = "PHQ-9 Item", y = "Race/Ethnicity", fill = "Coefficient",
           title = "Item-level race/ethnicity coefficients") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$itemTable <- renderTable({
    res <- item_results(); validate(need(nrow(res) > 0, "No item results"))
    res %>%
      select(item, term, estimate, std.error, p.value, q.value) %>%
      rename(
        "PHQ-9 Item" = item,
        "Race/Ethnicity" = term,
        "Estimate" = estimate,
        "Std. Error" = std.error,
        "P-value" = p.value,
        "BH q-value" = q.value
      )
  })
}

shinyApp(ui, server)
