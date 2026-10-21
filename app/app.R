## LC-MS Run Time & Cost Calculator
## Shinylive-compatible (uses only shiny + bslib)

library(shiny)
library(bslib)

# ---- helpers ---------------------------------------------------------------

fmt_dur <- function(min) {
  min <- round(min)
  d <- min %/% 1440
  h <- (min %% 1440) %/% 60
  m <- min %% 60
  paste0(if (d > 0) paste0(d, " d ") else "", h, " h ", m, " min")
}

fmt_money <- function(x) paste0("$", formatC(x, format = "f", digits = 2, big.mark = ","))

# Build one polarity's injection order. Blanks are spread evenly through the
# run (first at the start, last at the end once there are 2 or more).
build_sequence <- function(n, reps, n_blanks, qc_start, qc_end, qc_every) {
  blank_slots <- if (n_blanks > 0) round(seq(0, n, length.out = n_blanks)) else integer(0)
  type <- character(0); label <- character(0)
  qc_i <- 0; blk_i <- 0

  add <- function(t, l) { type <<- c(type, t); label <<- c(label, l) }
  add_blanks <- function(g) {
    for (j in seq_len(sum(blank_slots == g))) {
      blk_i <<- blk_i + 1
      add("Blank", sprintf("Blank_%02d", blk_i))
    }
  }

  add_blanks(0)
  for (j in seq_len(qc_start)) add("QC", sprintf("QC_start_%02d", j))
  for (i in seq_len(n)) {
    add(rep("Sample", reps),
        if (reps > 1) sprintf("S%03d_rep%d", i, seq_len(reps)) else sprintf("S%03d", i))
    if (i < n) {
      if (qc_every > 0 && i %% qc_every == 0) {
        qc_i <- qc_i + 1
        add("QC", sprintf("QC_%02d", qc_i))
      }
      add_blanks(i)
    }
  }
  for (j in seq_len(qc_end)) add("QC", sprintf("QC_end_%02d", j))
  add_blanks(n)

  data.frame(Type = type, Name = label, stringsAsFactors = FALSE)
}

num <- function(id, label, value, min = 0, step = 1) {
  numericInput(id, label, value = value, min = min, step = step)
}

inj_types <- c("Sample", "QC", "Blank")

# ---- UI --------------------------------------------------------------------

ui <- page_sidebar(
  title = "LC-MS Run Time & Cost Calculator",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  sidebar = sidebar(
    width = 330,
    accordion(
      open = c("Samples & method"),
      accordion_panel(
        "Samples & method",
        num("n_samples", "Number of samples", 48, min = 1),
        num("reps", "Injections per sample (technical reps)", 1, min = 1),
        num("gradient", "Gradient / method run time (min)", 15, step = 0.5),
        num("buffer", "Buffer per injection (min)", 1, step = 0.25),
        helpText("Buffer = autosampler draw, needle wash, re-equilibration, ",
                 "or any dead time between injections.")
      ),
      accordion_panel(
        "Polarity",
        radioButtons("polarity", NULL,
          choices = c("Single polarity" = "single",
                      "Pos & neg (separate injections)" = "both",
                      "Polarity switching (one injection)" = "switching"),
          selected = "both"),
        conditionalPanel(
          "input.polarity == 'single'",
          radioButtons("single_pol", NULL, c("Pos", "Neg"), inline = TRUE)
        ),
        conditionalPanel(
          "input.polarity == 'both'",
          radioButtons("pol_order", "Run order",
            choices = c("All pos, then all neg" = "batch",
                        "Alternate pos/neg per injection" = "alternate")),
          conditionalPanel(
            "input.pol_order == 'batch'",
            num("switchover", "Source switchover / stabilisation (min)", 0, step = 5)
          ),
          helpText("Every sample, QC and blank is injected once per polarity.")
        )
      ),
      accordion_panel(
        "Blanks & QCs",
        num("n_blanks", "Total blanks (per polarity)", 3),
        checkboxInput("blk_short", "Blanks use a different (shorter) method", FALSE),
        conditionalPanel(
          "input.blk_short",
          num("blk_time", "Blank method run time (min)", 5, step = 0.5)
        ),
        hr(),
        num("qc_start", "QCs at start (conditioning)", 3),
        num("qc_every", "QC every N samples (0 = none)", 10),
        num("qc_end", "QCs at end", 1),
        checkboxInput("qc_diff", "QCs use a different method", FALSE),
        conditionalPanel(
          "input.qc_diff",
          num("qc_time", "QC method run time (min)", 15, step = 0.5)
        )
      ),
      accordion_panel(
        "Instrument cost",
        num("rate", "Instrument rate ($ per hour)", 0, step = 5),
        num("fee", "Fixed booking / setup fee ($)", 0, step = 10),
        checkboxInput("round_up", "Bill in whole hours (round up)", FALSE)
      ),
      accordion_panel(
        "Extras & start time",
        num("setup", "One-off pre-sequence equilibration (min)", 0, step = 5),
        dateInput("start_date", "Start date", value = Sys.Date()),
        textInput("start_time", "Start time (HH:MM, 24 h)", "09:00")
      )
    )
  ),

  layout_columns(
    fill = FALSE,
    value_box(title = "Total run time", value = textOutput("total_time"),
              theme = "primary"),
    value_box(title = "Total injections", value = textOutput("total_inj"),
              p(textOutput("inj_note", inline = TRUE)), theme = "secondary"),
    value_box(title = "Estimated finish", value = textOutput("finish"),
              theme = "success"),
    value_box(title = "Instrument cost", value = textOutput("cost"),
              p(textOutput("cost_note", inline = TRUE)), theme = "warning")
  ),

  layout_columns(
    col_widths = c(5, 7),
    card(
      card_header("Breakdown"),
      tableOutput("breakdown"),
      plotOutput("bar", height = "220px")
    ),
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center",
        "Injection sequence",
        downloadButton("dl", "CSV", class = "btn-sm")
      ),
      div(style = "max-height: 520px; overflow-y: auto;",
          tableOutput("seq_table"))
    )
  )
)

# ---- server ----------------------------------------------------------------

server <- function(input, output, session) {

  clean <- function(x, default = 0) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x < 0) default else x
  }

  start_dt <- reactive({
    tm <- if (grepl("^\\d{1,2}:\\d{2}$", input$start_time)) input$start_time else "09:00"
    as.POSIXct(paste(input$start_date, tm), format = "%Y-%m-%d %H:%M", tz = "UTC")
  })

  seq_df <- reactive({
    req(input$n_samples >= 1)
    base <- build_sequence(
      n = floor(input$n_samples), reps = max(1, floor(clean(input$reps, 1))),
      n_blanks = floor(clean(input$n_blanks)),
      qc_start = floor(clean(input$qc_start)), qc_end = floor(clean(input$qc_end)),
      qc_every = floor(clean(input$qc_every))
    )

    grad <- clean(input$gradient)
    method_time <- c(
      Sample = grad,
      Blank  = if (isTRUE(input$blk_short)) clean(input$blk_time) else grad,
      QC     = if (isTRUE(input$qc_diff))   clean(input$qc_time)  else grad
    )
    base$Minutes <- unname(method_time[base$Type]) + clean(input$buffer)

    gap_row <- function(type, name, minutes)
      data.frame(Polarity = "", Type = type, Name = name, Minutes = minutes,
                 stringsAsFactors = FALSE)

    # Expand by polarity
    if (input$polarity == "single") {
      df <- cbind(Polarity = input$single_pol, base, stringsAsFactors = FALSE)
    } else if (input$polarity == "switching") {
      df <- cbind(Polarity = "Pos/Neg", base, stringsAsFactors = FALSE)
    } else if (input$pol_order == "alternate") {
      df <- base[rep(seq_len(nrow(base)), each = 2), ]
      df <- cbind(Polarity = rep(c("Pos", "Neg"), nrow(base)), df, stringsAsFactors = FALSE)
      df$Name <- paste0(df$Name, "_", tolower(df$Polarity))
    } else {
      pos <- cbind(Polarity = "Pos", base, stringsAsFactors = FALSE)
      neg <- cbind(Polarity = "Neg", base, stringsAsFactors = FALSE)
      pos$Name <- paste0(pos$Name, "_pos"); neg$Name <- paste0(neg$Name, "_neg")
      sw <- clean(input$switchover)
      df <- if (sw > 0)
        rbind(pos, gap_row("Switchover", "Polarity switchover", sw), neg)
      else rbind(pos, neg)
    }

    setup <- clean(input$setup)
    if (setup > 0) df <- rbind(gap_row("Equilibration", "Pre-sequence equilibration", setup), df)

    df$Position <- seq_len(nrow(df))
    df$Start <- format(start_dt() + 60 * cumsum(c(0, head(df$Minutes, -1))),
                       "%a %d %b %H:%M", tz = "UTC")
    rownames(df) <- NULL
    df[, c("Position", "Polarity", "Type", "Name", "Minutes", "Start")]
  })

  total_min <- reactive(sum(seq_df()$Minutes))
  n_inj     <- reactive(sum(seq_df()$Type %in% inj_types))

  billed_hours <- reactive({
    h <- total_min() / 60
    if (isTRUE(input$round_up)) ceiling(h) else h
  })
  cost <- reactive(billed_hours() * clean(input$rate) + clean(input$fee))

  output$total_time <- renderText(fmt_dur(total_min()))
  output$total_inj  <- renderText(n_inj())
  output$inj_note   <- renderText({
    s <- sum(seq_df()$Type == "Sample")
    paste0(s, " sample injections (", round(100 * s / n_inj()), "%)")
  })
  output$finish <- renderText(
    format(start_dt() + 60 * total_min(), "%a %d %b, %H:%M", tz = "UTC")
  )
  output$cost <- renderText(
    if (clean(input$rate) == 0 && clean(input$fee) == 0) "Enter a rate" else fmt_money(cost())
  )
  output$cost_note <- renderText({
    if (clean(input$rate) == 0 && clean(input$fee) == 0) return("Set $/h under Instrument cost")
    paste0(fmt_money(cost() / floor(input$n_samples)), " per sample · ",
           formatC(billed_hours(), format = "f", digits = if (isTRUE(input$round_up)) 0 else 2),
           " h billed")
  })

  breakdown <- reactive({
    df <- seq_df()
    types <- intersect(c(inj_types, "Switchover", "Equilibration"), unique(df$Type))
    out <- data.frame(
      Type = types,
      Injections = sapply(types, function(t) if (t %in% inj_types) sum(df$Type == t) else NA),
      Minutes = sapply(types, function(t) sum(df$Minutes[df$Type == t])),
      stringsAsFactors = FALSE
    )
    out$Time <- sapply(out$Minutes, fmt_dur)
    out$`% of run` <- sprintf("%.1f%%", 100 * out$Minutes / sum(out$Minutes))
    if (clean(input$rate) > 0)
      out$Cost <- fmt_money(out$Minutes / 60 * clean(input$rate))
    rownames(out) <- NULL
    out
  })

  output$breakdown <- renderTable(
    breakdown()[, setdiff(names(breakdown()), "Minutes")],
    na = "–", striped = TRUE, width = "100%", digits = 0
  )

  output$bar <- renderPlot({
    b <- breakdown()
    par(mar = c(4, 7, 1, 1))
    cols <- c(Sample = "#2C3E50", QC = "#18BC9C", Blank = "#95A5A6",
              Switchover = "#E74C3C", Equilibration = "#F39C12")[b$Type]
    barplot(rev(b$Minutes / 60), names.arg = rev(b$Type), horiz = TRUE, las = 1,
            col = rev(cols), border = NA, xlab = "Hours")
  })

  output$seq_table <- renderTable(
    transform(seq_df(), Minutes = sprintf("%.2f", Minutes)),
    striped = TRUE, width = "100%"
  )

  output$dl <- downloadHandler(
    filename = function() paste0("LCMS_sequence_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(seq_df(), file, row.names = FALSE)
  )
}

shinyApp(ui, server)
