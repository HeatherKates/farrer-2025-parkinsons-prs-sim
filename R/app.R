library(shiny)
library(readxl)
library(ggplot2)
library(DT)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  titlePanel("Interactive GWAS Polygenic Risk Score (PRS) Simulator"),
  tabsetPanel(
    tabPanel("Simulation",
             sidebarLayout(
               sidebarPanel(
                 h4("1. Upload data"),
                 HTML('
<p>Your file can include as many columns as you want. The columns used for simulation will be available to select after your data is loaded. The <b>first sheet</b> will be used.</p>
'),
                 fileInput("xlsx_file", "Upload Summary Statistics (.xlsx)", accept = ".xlsx"),
                 uiOutput("choose_cols"),
                 h4("2. Filter data to be used for simulation"),
                 HTML('
<p>Your filtered table will update based on the filters set here. When your filters are ready, click <code>Filter Data</code>.</p>
'),
                 textAreaInput("rsid_drop", "Exclude rsIDs (paste, separated by comma/space/newline):", height = "60px"),
                 uiOutput("filter_beta"),
                 uiOutput("filter_af"),
                 uiOutput("factor_key"),
                 uiOutput("filter_subset_var"),
                 actionButton("apply_filter", "Filter Data", class = "btn-success"),
                 br(),
                 downloadButton("download_filtered", "Download Filtered Data"),
                 h4("3. Simulate genotypes and calculate PRS"),
                 HTML('
<p>We simulate genotypes for <i>n</i> synthetic individuals assuming Hardy-Weinberg equilibrium using the effect allele frequencies. For each SNP:</p>
<ul>
<li>Genotype is drawn from a binomial distribution with 2 trials and success probability = Effect Allele Frequency.</li>
<li>Each simulated individual\'s PRS is then computed as: <code>sum(Beta * Genotype)</code></li>
</ul>
'),
                 numericInput("n_inds", "Simulated Population Size", value = 100000, min = 10, step = 100),
                 actionButton("run_sim", "Run Simulation", class = "btn-primary")
               ),
               mainPanel(
                 h4("Preview of Uploaded Data"),
                 DTOutput("table_raw"),
                 hr(),
                 h4("Filtered Data"),
                 DTOutput("table_filtered"),
                 hr(),
                 h4("PRS Simulation Summary"),
                 verbatimTextOutput("stats_out"),
                 plotOutput("prs_plot"),
                 h4("OR Simulation Summary"),
                 verbatimTextOutput("or_stats_out"),
                 plotOutput("or_plot")
               )
             )
    ),
    tabPanel("About",
             h3("About this Application"),
             HTML('
<p>
This app simulates polygenic risk scores (PRS) from GWAS summary statistics provided by the user.<br>
Please see the sidebar for instructions.<br>
</p>
<p>
<b>Developed by:</b> Dr. Heather Kates, UF Health Cancer Center BCB-SR<br>
<b>Contact:</b> hkates@ufl.edu
</p>
')
    )
  )
)

server <- function(input, output, session) {
  dat <- reactive({
    req(input$xlsx_file)
    as.data.frame(read_excel(input$xlsx_file$datapath, sheet = 1))
  })
  
  output$choose_cols <- renderUI({
    req(dat())
    cn <- colnames(dat())
    tagList(
      selectInput("rsid_col", "rsID Column",
                  choices = cn,
                  selected = if (!is.null(input$rsid_col) && input$rsid_col %in% cn) input$rsid_col
                  else grep("rsid", cn, ignore.case = TRUE, value = TRUE)[1]
      ),
      selectInput("beta_col", "Beta Column",
                  choices = cn,
                  selected = if (!is.null(input$beta_col) && input$beta_col %in% cn) input$beta_col
                  else grep("beta", cn, ignore.case = TRUE, value = TRUE)[1]
      ),
      selectInput("af_col", "Effect Allele Frequency Column",
                  choices = cn,
                  selected = if (!is.null(input$af_col) && input$af_col %in% cn) input$af_col
                  else grep("freq", cn, ignore.case = TRUE, value = TRUE)[1]
      ),
      selectInput("subset_var", "Subset SNPs (optional):",
                  choices = setdiff(cn, c(input$rsid_col, input$beta_col, input$af_col)),
                  selected = if (!is.null(input$subset_var) && input$subset_var %in% cn) input$subset_var
                  else if ("Known or Novel" %in% cn) "Known or Novel"
                  else setdiff(cn, c(input$rsid_col, input$beta_col, input$af_col))[1]
      )
    )
  })
  
  output$table_raw <- renderDT({
    req(dat())
    DT::datatable(
      dat(),
      rownames=FALSE, options = list(pageLength=10, scrollX=TRUE)
    )
  })
  
  output$filter_beta <- renderUI({
    req(dat(), input$beta_col)
    d <- dat()
    beta_vals <- suppressWarnings(as.numeric(d[[input$beta_col]]))
    beta_vals <- beta_vals[!is.na(beta_vals) & is.finite(beta_vals)]
    if (length(beta_vals) > 0) {
      rng <- range(beta_vals)
      sliderInput("beta_slider", "Beta Range:", min = rng[1], max = rng[2], value = rng)
    } else {
      helpText("No numeric values detected in Beta column.")
    }
  })
  
  output$filter_af <- renderUI({
    req(dat(), input$af_col)
    d <- dat()
    af_vals <- suppressWarnings(as.numeric(d[[input$af_col]]))
    af_vals <- af_vals[!is.na(af_vals) & is.finite(af_vals)]
    if (length(af_vals) > 0) {
      rng <- range(af_vals)
      sliderInput("af_slider", "Effect Allele Frequency Range:",
                  min = rng[1], max = rng[2], value = rng)
    } else {
      helpText("No numeric values detected in Effect Allele Frequency column.")
    }
  })
  
  output$factor_key <- renderUI({
    req(dat(), input$subset_var)
    vals <- dat()[[input$subset_var]]
    if (!is.numeric(vals)) {
      fac <- as.factor(vals)
      levels_str <- paste(sprintf("<b>%d</b> = %s", seq_along(levels(fac)), levels(fac)), collapse = "<br>")
      HTML(sprintf("<b>Level code mapping for '%s':</b><br>%s", input$subset_var, levels_str))
    }
  })
  
  output$filter_subset_var <- renderUI({
    req(dat(), input$subset_var)
    vals <- dat()[[input$subset_var]]
    asnum <- as.numeric(as.factor(vals))
    rng <- range(asnum, na.rm = TRUE)
    sliderInput("subset_var_slider",
                sprintf("Filter %s (numeric code):", input$subset_var),
                min = rng[1], max = rng[2], value = rng, step = 1)
  })
  output$download_filtered <- downloadHandler(
    filename = function() {
      paste0("filtered_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      fd <- tryCatch(filtered_data(), error=function(e) data.frame())
      if (is.null(fd) || nrow(fd) == 0) {
        write.csv(data.frame(), file, row.names = FALSE)
      } else {
        write.csv(fd, file, row.names = FALSE)
      }
    }
  )
  
  filter_done <- reactiveVal(FALSE)
  
  filtered_data <- eventReactive(input$apply_filter, {
    req(dat(), input$rsid_col, input$beta_col, input$af_col)
    d <- dat()
    d <- d[complete.cases(d[, c(input$rsid_col, input$beta_col, input$af_col)]), , drop = FALSE]
    # Remove rsIDs to drop
    to_drop <- unique(unlist(strsplit(input$rsid_drop, "[, \n]+")))
    if (length(to_drop) > 0 && any(to_drop != "")) {
      d <- d[!(as.character(d[[input$rsid_col]]) %in% to_drop), , drop = FALSE]
    }
    # Filter by beta/freq
    beta_vals <- suppressWarnings(as.numeric(d[[input$beta_col]]))
    af_vals <- suppressWarnings(as.numeric(d[[input$af_col]]))
    keep <- rep(TRUE, nrow(d))
    if (!is.null(input$beta_slider) && length(beta_vals) == nrow(d)) {
      keep <- keep & beta_vals >= input$beta_slider[1] & beta_vals <= input$beta_slider[2]
    }
    if (!is.null(input$af_slider) && length(af_vals) == nrow(d)) {
      keep <- keep & af_vals >= input$af_slider[1] & af_vals <= input$af_slider[2]
    }
    # Subset SNPs (always as numeric code)
    if (!is.null(input$subset_var) && input$subset_var %in% names(d) && !is.null(input$subset_var_slider)) {
      vals <- d[[input$subset_var]]
      codes <- as.numeric(as.factor(vals))
      keep <- keep & codes >= input$subset_var_slider[1] & codes <= input$subset_var_slider[2]
    }
    d_filtered <- d[keep, , drop = FALSE]
    if (nrow(d_filtered) > 0) filter_done(TRUE) else filter_done(FALSE)
    d_filtered
  })
  
  observe({
    if (isTRUE(filter_done())) enable("run_sim") else disable("run_sim")
  })
  observeEvent(input$xlsx_file, { filter_done(FALSE) })
  observeEvent(list(input$rsid_col, input$beta_col, input$af_col), { filter_done(FALSE) })
  observeEvent(dat(), {
    cn <- colnames(dat())
    updateSelectInput(session, "rsid_col", choices = cn,
                      selected = grep("rsid", cn, ignore.case = TRUE, value = TRUE)[1])
    updateSelectInput(session, "beta_col", choices = cn,
                      selected = grep("beta", cn, ignore.case = TRUE, value = TRUE)[1])
    updateSelectInput(session, "af_col", choices = cn,
                      selected = grep("freq", cn, ignore.case = TRUE, value = TRUE)[1])
    subset_choices <- setdiff(cn, c(input$rsid_col, input$beta_col, input$af_col))
    updateSelectInput(session, "subset_var", choices = subset_choices,
                      selected = if ("Known or Novel" %in% cn) "Known or Novel" else subset_choices[1])
  })
  
  output$table_filtered <- renderDT({
    filtered <- filtered_data()
    req(filtered)
    if (nrow(filtered) == 0) {
      DT::datatable(data.frame())
    } else {
      DT::datatable(filtered[, c(input$rsid_col, input$beta_col, input$af_col), drop = FALSE],
                    rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
    }
  })
  
  simulate_and_get_prs <- function(prs_table, rsid_col, beta_col, af_col, n_individuals) {
    prs_table <- prs_table[complete.cases(prs_table[, c(rsid_col, beta_col, af_col)]), , drop = FALSE]
    AF <- suppressWarnings(as.numeric(prs_table[[af_col]]))
    Beta <- suppressWarnings(as.numeric(prs_table[[beta_col]]))
    keep <- !is.na(AF) & !is.na(Beta) & AF >= 0 & AF <= 1
    AF <- AF[keep]; Beta <- Beta[keep]; rsids <- prs_table[[rsid_col]][keep]
    genos <- sapply(AF, function(freq) rbinom(n_individuals, 2, freq))
    colnames(genos) <- rsids
    PRS <- as.vector(genos %*% Beta)
    data.frame(PRS = PRS)
  }
  simres <- eventReactive(input$run_sim, {
    req(filtered_data())
    simulate_and_get_prs(filtered_data(), input$rsid_col, input$beta_col, input$af_col, as.integer(input$n_inds))
  })
  output$stats_out <- renderPrint({
    req(simres())
    x <- simres()$PRS
    cat(sprintf("Mean PRS: %.3f", mean(x)), "\n")
    cat(sprintf("SD PRS: %.3f", sd(x)), "\n")
    cat(sprintf("Min PRS: %.3f", min(x)), "\n")
    cat(sprintf("Max PRS: %.3f", max(x)), "\n")
  })
  output$prs_plot <- renderPlot({
    req(simres())
    x <- simres()$PRS
    m <- mean(x)
    s <- sd(x)
    d <- density(x)
    y_max <- max(d$y)
    x_pos <- m + 0.02 * (max(x) - min(x))
    y_pos <- y_max * 0.7
    ggplot(data.frame(PRS = x), aes(x = PRS)) +
      geom_density(fill = "#3182bd", alpha = 0.4) +
      geom_vline(xintercept = m, linetype = "dashed", color = "black", linewidth = 1) +
      annotate(
        "text", x = x_pos, y = y_pos,
        label = paste0("Mean = ", sprintf("%.2f", m), "\nSD = ", sprintf("%.2f", s)),
        hjust = 0, vjust = 0, color = "black", size = 4
      ) +
      labs(title = "Simulated PRS Distribution", x = "PRS", y = "Density") +
      theme_minimal()
  })
  output$or_stats_out <- renderPrint({
    req(simres())
    or <- exp(simres()$PRS)
    cat(sprintf("Mean OR: %.3f", mean(or)), "\n")
    cat(sprintf("SD OR: %.3f", sd(or)), "\n")
    cat(sprintf("Min OR: %.3f", min(or)), "\n")
    cat(sprintf("Max OR: %.3f", max(or)), "\n")
  })
  output$or_plot <- renderPlot({
    req(simres())
    or <- exp(simres()$PRS)
    m <- mean(or)
    s <- sd(or)
    d <- density(or)
    y_max <- max(d$y)
    x_pos <- m + 0.02 * (max(or) - min(or))
    y_pos <- y_max * 0.7
    ggplot(data.frame(OR = or), aes(x = OR)) +
      geom_density(fill = "#de2d26", alpha = 0.4) +
      geom_vline(xintercept = m, linetype = "dashed", color = "black", linewidth = 1) +
      annotate(
        "text", x = x_pos, y = y_pos,
        label = paste0("Mean = ", sprintf("%.2f", m), "\nSD = ", sprintf("%.2f", s)),
        hjust = 0, vjust = 0, color = "black", size = 4
      ) +
      labs(title = "Simulated Odds Ratio (OR) Distribution",
           x = "Odds Ratio (exp(PRS))", y = "Density") +
      theme_minimal()
  })
}

shinyApp(ui, server)