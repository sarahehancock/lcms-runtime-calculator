## LC-MS Run Time, Cost & Solvent Calculator
## Shinylive-compatible (uses only shiny + bslib, plus jsonlite which ships with shiny)

library(shiny)
library(bslib)

# ---- built-in LC methods ---------------------------------------------------
# Each method is a pump timetable: time (min), %B, flow (mL/min).
# Composition and flow change linearly between rows, as on the LC pump.
# For a step change, add two rows close together in time (e.g. 25 and 25.1).

builtin_methods <- list(
  "Lipidomics (RP, 30 min)" = data.frame(
    time = c(0, 1.5, 4, 5, 8, 11, 14, 18, 21, 25, 25.1, 30),
    B    = c(32, 32, 45, 52, 58, 66, 70, 75, 97, 97, 32, 32),
    flow = 0.26
  ),
  # BEH Z-HILIC method from Girel et al. 2023, J Chromatogr A 1697:463994.
  # "gradient at 350 uL/min and equilibration stage at 450 uL/min" (Discussion).
  "Metabolomics (Z-HILIC, Girel 2023)" = data.frame(
    time = c(0, 2, 20, 20.5, 24, 24.1, 27, 27.1, 34.5),
    B    = c(0, 0, 72, 100, 100, 0, 0, 0, 0),
    flow = c(0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.35, 0.45, 0.45)
  )
)

# Column used with each built-in method: length & ID (mm), porosity, and column
# volumes for pre-sequence equilibration. Selecting a method fills these in.
builtin_columns <- list(
  # Waters ACQUITY UPLC CSH C18, 2.1 x 100 mm, 1.7 um (gradient from Hu et al. 2008,
  # J Proteome Res 7:4982, originally on an Ascentis Express C8 2.1 x 150 mm)
  "Lipidomics (RP, 30 min)"            = c(len = 100, id = 2.1, porosity = 0.65, n_cv = 20),
  # Waters Atlantis Premier BEH Z-HILIC, 2.1 x 150 mm, 1.7 um; 50 column volumes (paper)
  "Metabolomics (Z-HILIC, Girel 2023)" = c(len = 150, id = 2.1, porosity = 0.65, n_cv = 50)
)

# ---- helpers ---------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

fmt_dur <- function(min) {
  min <- round(min)
  d <- min %/% 1440
  h <- (min %% 1440) %/% 60
  m <- min %% 60
  paste0(if (d > 0) paste0(d, " d ") else "", h, " h ", m, " min")
}

fmt_money <- function(x) paste0("$", formatC(x, format = "f", digits = 2, big.mark = ","))

fmt_vol <- function(ml) {
  if (ml >= 1000) sprintf("%.2f L", ml / 1000)
  else if (ml >= 10) sprintf("%.0f mL", ml)
  else sprintf("%.2f mL", ml)
}

method_to_text <- function(df) {
  paste(sprintf("%g, %g, %g", df$time, df$B, df$flow), collapse = "\n")
}

# Parse "time, %B, flow" lines. Returns list(ok, df, msg).
parse_method <- function(txt) {
  lines <- trimws(strsplit(txt %||% "", "\n")[[1]])
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  if (length(lines) < 2) return(list(ok = FALSE, msg = "Enter at least two rows: time, %B, flow"))
  rows <- lapply(lines, function(l) suppressWarnings(as.numeric(strsplit(l, "[,;\t ]+")[[1]])))
  bad <- which(vapply(rows, function(r) length(r) != 3 || anyNA(r), logical(1)))
  if (length(bad)) return(list(ok = FALSE, msg = paste0(
    "Row ", bad[1], " needs three numbers (time, %B, flow): \"", lines[bad[1]], "\"")))
  df <- as.data.frame(do.call(rbind, rows))
  names(df) <- c("time", "B", "flow")
  if (any(diff(df$time) < 0)) return(list(ok = FALSE, msg = "Times must be in increasing order"))
  if (any(df$B < 0 | df$B > 100)) return(list(ok = FALSE, msg = "%B must be between 0 and 100"))
  if (any(df$flow < 0)) return(list(ok = FALSE, msg = "Flow cannot be negative"))
  if (max(df$time) <= 0) return(list(ok = FALSE, msg = "Method must be longer than 0 min"))
  list(ok = TRUE, df = df, msg = NULL)
}

# Solvent used by one run of a timetable (mL). %B and flow are linear within
# each segment, so the integral of flow x %B is computed exactly.
gradient_volumes <- function(df) {
  if (df$time[1] > 0) df <- rbind(df[1, ], df)
  if (df$time[1] > 0) df$time[1] <- 0
  if (nrow(df) < 2) return(c(A = 0, B = 0, total = 0))
  i  <- seq_len(nrow(df) - 1)
  dt <- diff(df$time)
  F1 <- df$flow[i]; F2 <- df$flow[i + 1]
  B1 <- df$B[i] / 100; B2 <- df$B[i + 1] / 100
  total <- sum(dt * (F1 + F2) / 2)
  vb    <- sum(dt * (F1 * B1 / 3 + F1 * B2 / 6 + F2 * B1 / 6 + F2 * B2 / 3))
  c(A = total - vb, B = vb, total = total)
}

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

MANUAL <- "Manual entry"   # LC method option: type run time, flow and average %B

# Browser storage for custom methods (per browser; use Export/Import to move them)
storage_js <- tags$script(HTML("
  $(document).on('shiny:connected', function() {
    var s = '';
    try { s = localStorage.getItem('lcms_custom_methods') || ''; } catch (e) {}
    Shiny.setInputValue('stored_methods', s, {priority: 'event'});
  });
  Shiny.addCustomMessageHandler('save_methods', function(x) {
    try { localStorage.setItem('lcms_custom_methods', x); } catch (e) {}
  });
"))

# ---- UI --------------------------------------------------------------------

ui <- page_navbar(
  title = "LC-MS Run Calculator",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  header = storage_js,
  fillable = FALSE,

  sidebar = sidebar(
    width = 330,
    accordion(
      open = c("Samples & method"),
      accordion_panel(
        "Samples & method",
        selectInput("method", "LC method", choices = c(names(builtin_methods), MANUAL)),
        conditionalPanel(
          sprintf("input.method != '%s'", MANUAL),
          checkboxInput("sync_gradient", "Set run time from LC method", TRUE)
        ),
        conditionalPanel(
          sprintf("input.method == '%s'", MANUAL),
          layout_columns(
            num("man_flow", "Flow (mL/min)", 0.3, step = 0.01),
            num("man_B", "Average %B", 50, step = 5)
          ),
          helpText("Solvent use = run time × flow, split by the average %B. ",
                   "Set the column under Column equilibration.")
        ),
        num("n_samples", "Number of samples", 48, min = 1),
        num("reps", "Injections per sample (technical reps)", 1, min = 1),
        num("gradient", "Gradient / method run time (min)", 30, step = 0.5),
        num("buffer", "Buffer per injection (min)", 1, step = 0.25),
        helpText("Buffer = autosampler draw, needle wash, or any dead time ",
                 "between injections.")
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
        "Solvents",
        num("margin", "Safety margin (%)", 10, step = 5),
        num("prime", "Prime / purge volume per channel (mL)", 20, step = 5),
        num("bottle", "Bottle size (mL)", 1000, min = 1, step = 250),
        checkboxInput("idle_flow",
          "Count flow during buffer and switchover time (at initial conditions)",
          TRUE)
      ),
      accordion_panel(
        "Instrument cost",
        num("rate", "Instrument rate ($ per hour)", 0, step = 5),
        num("fee", "Fixed booking / setup fee ($)", 0, step = 10),
        checkboxInput("round_up", "Bill in whole hours (round up)", FALSE)
      ),
      accordion_panel(
        "Column equilibration",
        radioButtons("eq_mode", "Pre-sequence equilibration",
          choices = c("Column volumes" = "cv", "Fixed time" = "min", "None" = "none"),
          selected = "cv", inline = TRUE),
        conditionalPanel(
          "input.eq_mode == 'cv'",
          layout_columns(
            num("col_len", "Column length (mm)", 100, min = 1, step = 50),
            num("col_id", "Column ID (mm)", 2.1, min = 0.1, step = 0.1)
          ),
          layout_columns(
            num("n_cv", "Column volumes", 20, step = 5),
            num("porosity", "Porosity", 0.65, step = 0.05)
          ),
          helpText("Filled in from the selected LC method; saving a method stores these too. ",
                   "20-50 column volumes is typical; HILIC usually needs the high end. ",
                   "Porosity ~0.65 for fully porous particles, ~0.5 for core-shell. ",
                   "Runs at the method's initial %B and flow.")
        ),
        conditionalPanel(
          "input.eq_mode == 'min'",
          num("setup", "Equilibration time (min)", 30, step = 5)
        ),
        uiOutput("eq_summary")
      ),
      accordion_panel(
        "Start time",
        dateInput("start_date", "Start date", value = Sys.Date()),
        textInput("start_time", "Start time (HH:MM, 24 h)", "09:00")
      )
    )
  ),

  # ---- Tab 1: run time & cost ----
  nav_panel(
    "Run time & cost",
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
  ),

  # ---- Tab 2: solvents ----
  nav_panel(
    "Solvents",
    layout_columns(
      fill = FALSE,
      value_box(title = "Solvent A needed", value = textOutput("solv_a"),
                p(textOutput("solv_a_note", inline = TRUE)), theme = "primary"),
      value_box(title = "Solvent B needed", value = textOutput("solv_b"),
                p(textOutput("solv_b_note", inline = TRUE)), theme = "info"),
      value_box(title = "Per injection", value = textOutput("solv_per"),
                p(textOutput("solv_per_note", inline = TRUE)), theme = "secondary")
    ),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Solvent breakdown"),
        tableOutput("solv_table"),
        helpText("Every injection is assumed to run the selected gradient, including ",
                 "blanks and QCs. Column equilibration, buffer and switchover time ",
                 "are counted at the method's initial %B and flow.")
      ),
      card(
        card_header("Gradient"),
        plotOutput("grad_plot", height = "260px")
      )
    ),
    conditionalPanel(
      sprintf("input.method != '%s'", MANUAL),
      card(
      card_header("Method editor"),
      layout_columns(
        col_widths = c(5, 7),
        div(
          textAreaInput("method_text", "Timetable: time (min), %B, flow (mL/min) - one row per line",
                        value = method_to_text(builtin_methods[[1]]), rows = 13, width = "100%"),
          uiOutput("method_msg")
        ),
        div(
          helpText("Edits are used in the calculation straight away. %B and flow change ",
                   "linearly between rows, as on the pump. For a step change, use two rows ",
                   "close together (e.g. 25, 97 then 25.1, 32). Times before the first row ",
                   "are held at the first row's conditions."),
          textInput("method_name", "Method name", value = names(builtin_methods)[1]),
          div(class = "d-flex flex-wrap gap-2 mb-3",
              actionButton("save_method", "Save method", class = "btn-primary btn-sm"),
              actionButton("delete_method", "Delete method", class = "btn-outline-danger btn-sm")),
          helpText("Saved methods are kept in this browser. Built-in methods can't be ",
                   "overwritten - save edits under a new name. Use export/import to move ",
                   "methods to another computer or share them."),
          div(class = "d-flex flex-wrap gap-2 align-items-start",
              downloadButton("export_methods", "Export my methods", class = "btn-sm"),
              fileInput("import_methods", NULL, accept = ".json",
                        buttonLabel = "Import...", placeholder = "methods .json"))
        )
      )
      )
    )
  )
)

# ---- server ----------------------------------------------------------------

server <- function(input, output, session) {

  clean <- function(x, default = 0) {
    if (is.null(x) || length(x) == 0 || is.na(x) || x < 0) default else x
  }

  # ---- LC methods ----
  # Custom methods (timetables) and their column settings, stored in the browser
  # as JSON: {"methods": {name: [rows]}, "columns": {name: {len, id, porosity, n_cv}}}
  custom      <- reactiveVal(list())
  custom_cols <- reactiveVal(list())
  all_methods <- reactive(c(builtin_methods, custom()))
  all_cols    <- reactive(c(builtin_columns, custom_cols()))

  methods_from_json <- function(txt) {
    x <- tryCatch(jsonlite::fromJSON(txt, simplifyDataFrame = TRUE), error = function(e) NULL)
    if (!is.list(x)) return(list(methods = list(), columns = list()))
    if (!is.null(x[["methods"]])) { meths <- x[["methods"]]; cols <- x[["columns"]] %||% list() }
    else { meths <- x; cols <- list() }                     # older methods-only format
    meths <- lapply(meths, function(d) {
      d <- tryCatch(as.data.frame(d), error = function(e) NULL)
      if (is.null(d) || !all(c("time", "B", "flow") %in% names(d))) return(NULL)
      d <- d[, c("time", "B", "flow")]
      d[] <- lapply(d, as.numeric)
      d
    })
    meths <- meths[!vapply(meths, is.null, logical(1))]
    cols <- lapply(cols, function(cc) {
      v <- suppressWarnings(as.numeric(unlist(cc)[c("len", "id", "porosity", "n_cv")]))
      if (anyNA(v)) return(NULL)
      setNames(v, c("len", "id", "porosity", "n_cv"))
    })
    cols <- cols[!vapply(cols, is.null, logical(1))]
    keep <- setdiff(names(meths), c(names(builtin_methods), MANUAL))
    list(methods = meths[keep], columns = cols[intersect(names(cols), keep)])
  }

  methods_json <- function(pretty = FALSE) {
    as.character(jsonlite::toJSON(
      list(methods = custom(), columns = lapply(custom_cols(), as.list)),
      digits = NA, auto_unbox = TRUE, pretty = pretty))
  }

  store_custom <- function(x, cols, select = NULL) {
    custom(x); custom_cols(cols[intersect(names(cols), names(x))])
    session$sendCustomMessage("save_methods", methods_json())
    updateSelectInput(session, "method", choices = c(names(all_methods()), MANUAL),
                      selected = select %||% input$method)
  }

  observeEvent(input$stored_methods, {
    if (nzchar(input$stored_methods)) {
      x <- methods_from_json(input$stored_methods)
      custom(x$methods); custom_cols(x$columns)
      updateSelectInput(session, "method", choices = c(names(all_methods()), MANUAL),
                        selected = input$method)
    }
  }, once = TRUE)

  observeEvent(input$method, {
    m <- all_methods()[[input$method]]
    req(m)
    updateTextAreaInput(session, "method_text", value = method_to_text(m))
    updateTextInput(session, "method_name", value = input$method)
    cc <- all_cols()[[input$method]]
    if (!is.null(cc)) {
      updateNumericInput(session, "col_len",  value = cc[["len"]])
      updateNumericInput(session, "col_id",   value = cc[["id"]])
      updateNumericInput(session, "porosity", value = cc[["porosity"]])
      updateNumericInput(session, "n_cv",     value = cc[["n_cv"]])
    }
  })

  is_manual <- reactive(identical(input$method, MANUAL))

  cur_method <- reactive({
    if (is_manual()) {
      rt <- clean(input$gradient)
      if (rt <= 0) return(list(ok = FALSE, msg = "Enter a run time above 0 min"))
      B  <- min(100, clean(input$man_B))
      return(list(ok = TRUE, msg = NULL,
                  df = data.frame(time = c(0, rt), B = B, flow = clean(input$man_flow))))
    }
    parse_method(input$method_text)
  })

  output$method_msg <- renderUI({
    m <- cur_method()
    if (!m$ok) return(div(class = "text-danger small", m$msg))
    saved <- all_methods()[[input$method]]
    changed <- is.null(saved) || nrow(saved) != nrow(m$df) ||
      any(abs(as.matrix(saved) - as.matrix(m$df)) > 1e-9)
    div(class = "small text-muted",
        sprintf("%d rows, %g min", nrow(m$df), max(m$df$time)),
        if (changed) span(class = "text-warning ms-2", "Unsaved edits"))
  })

  observeEvent(input$save_method, {
    m <- cur_method()
    if (!m$ok) return(showNotification(m$msg, type = "error"))
    nm <- trimws(input$method_name)
    if (!nzchar(nm)) return(showNotification("Give the method a name", type = "error"))
    if (nm %in% c(names(builtin_methods), MANUAL))
      return(showNotification("That's a built-in method name - choose a new name to save your version.",
                              type = "warning"))
    x <- custom(); x[[nm]] <- m$df
    cols <- custom_cols()
    cols[[nm]] <- c(len = clean(input$col_len), id = clean(input$col_id),
                    porosity = clean(input$porosity), n_cv = clean(input$n_cv))
    store_custom(x, cols, select = nm)
    showNotification(paste0("Saved \"", nm, "\" with the current column settings"), type = "message")
  })

  observeEvent(input$delete_method, {
    nm <- input$method
    if (nm %in% c(names(builtin_methods), MANUAL))
      return(showNotification("Built-in methods can't be deleted", type = "warning"))
    x <- custom(); x[[nm]] <- NULL
    store_custom(x, custom_cols(), select = names(builtin_methods)[1])
    showNotification(paste0("Deleted \"", nm, "\""), type = "message")
  })

  output$export_methods <- downloadHandler(
    filename = function() paste0("lcms_methods_", Sys.Date(), ".json"),
    content  = function(file) writeLines(methods_json(pretty = TRUE), file)
  )

  observeEvent(input$import_methods, {
    txt <- paste(readLines(input$import_methods$datapath, warn = FALSE), collapse = "\n")
    x <- methods_from_json(txt)
    if (!length(x$methods)) return(showNotification("No valid methods found in that file", type = "error"))
    y <- custom(); y[names(x$methods)] <- x$methods
    cols <- custom_cols(); cols[names(x$columns)] <- x$columns
    store_custom(y, cols, select = names(x$methods)[1])
    showNotification(sprintf("Imported %d method(s)", length(x$methods)), type = "message")
  })

  observe({
    req(isTRUE(input$sync_gradient), !is_manual())
    m <- cur_method(); req(m$ok)
    updateNumericInput(session, "gradient", value = max(m$df$time))
  })

  # ---- column equilibration ----
  # Column volume (mL) = pi * r^2 * L * porosity, with r and L in cm.
  eq <- reactive({
    m  <- cur_method()
    F0 <- if (m$ok) m$df$flow[1] else 0
    if (input$eq_mode == "cv") {
      v0  <- pi * (clean(input$col_id) / 20)^2 * (clean(input$col_len) / 10) * clean(input$porosity)
      vol <- clean(input$n_cv) * v0
      min <- if (F0 > 0) vol / F0 else 0
      list(min = min, v0 = v0, F0 = F0)
    } else if (input$eq_mode == "min") {
      list(min = clean(input$setup), v0 = NA, F0 = F0)
    } else list(min = 0, v0 = NA, F0 = F0)
  })

  output$eq_summary <- renderUI({
    e <- eq()
    if (input$eq_mode == "none") return(NULL)
    if (e$F0 <= 0) return(div(class = "small text-danger",
      "The method's initial flow is 0, so equilibration can't be calculated."))
    txt <- if (input$eq_mode == "cv")
      sprintf("1 column volume = %.0f µL. %g CV = %s at %g mL/min → %s.",
              1000 * e$v0, clean(input$n_cv), fmt_vol(e$min * e$F0), e$F0, fmt_dur(e$min))
    else sprintf("%s at %g mL/min = %s.", fmt_dur(e$min), e$F0, fmt_vol(e$min * e$F0))
    div(class = "small text-muted", txt)
  })

  # ---- sequence & run time ----
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

    setup <- eq()$min
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

  # ---- solvents ----
  solv <- reactive({
    m <- cur_method()
    validate(need(m$ok, m$msg))
    df  <- m$df
    per <- gradient_volumes(df)
    n   <- n_inj()

    idle_min <- if (isTRUE(input$idle_flow)) {
      clean(input$buffer) * n + sum(seq_df()$Minutes[seq_df()$Type == "Switchover"])
    } else 0
    F0 <- df$flow[1]; B0 <- df$B[1] / 100
    at_initial <- function(min) c(A = min * F0 * (1 - B0), B = min * F0 * B0)
    idle   <- at_initial(idle_min)
    eq_min <- eq()$min
    equil  <- at_initial(eq_min)

    run    <- per[c("A", "B")] * n
    sub    <- run + idle + equil
    margin <- sub * clean(input$margin) / 100
    prime  <- c(A = clean(input$prime), B = clean(input$prime))
    list(per = per, n = n, run = run, idle = idle, idle_min = idle_min,
         equil = equil, eq_min = eq_min,
         margin = margin, prime = prime, total = sub + margin + prime)
  })

  bottles <- function(ml) {
    b <- max(1, clean(input$bottle, 1000))
    k <- ceiling(ml / b)
    sprintf("%d × %s bottle%s", k, fmt_vol(b), if (k == 1) "" else "s")
  }

  output$solv_a      <- renderText(fmt_vol(solv()$total[["A"]]))
  output$solv_b      <- renderText(fmt_vol(solv()$total[["B"]]))
  output$solv_a_note <- renderText(bottles(solv()$total[["A"]]))
  output$solv_b_note <- renderText(bottles(solv()$total[["B"]]))
  output$solv_per    <- renderText(fmt_vol(solv()$per[["total"]]))
  output$solv_per_note <- renderText(sprintf("A %s · B %s per run",
    fmt_vol(solv()$per[["A"]]), fmt_vol(solv()$per[["B"]])))

  output$solv_table <- renderTable({
    s <- solv()
    rows <- list(
      c("Gradient, per injection",                  s$per[["A"]],   s$per[["B"]]),
      c(sprintf("Gradient × %d injections", s$n),   s$run[["A"]],   s$run[["B"]]),
      c(sprintf("Column equilibration (%s)", fmt_dur(s$eq_min)),
                                                    s$equil[["A"]], s$equil[["B"]]),
      c(sprintf("Buffer & switchover at initial conditions (%s)", fmt_dur(s$idle_min)),
                                                    s$idle[["A"]],  s$idle[["B"]]),
      c(sprintf("Safety margin (%g%%)", clean(input$margin)),
                                                    s$margin[["A"]], s$margin[["B"]]),
      c("Prime / purge",                            s$prime[["A"]], s$prime[["B"]]),
      c("Total to prepare",                         s$total[["A"]], s$total[["B"]])
    )
    out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    names(out) <- c("Component", "A (mL)", "B (mL)")
    out[[2]] <- as.numeric(out[[2]]); out[[3]] <- as.numeric(out[[3]])
    out$`Total (mL)` <- out[[2]] + out[[3]]
    out
  }, striped = TRUE, width = "100%", digits = 1)

  output$grad_plot <- renderPlot({
    m <- cur_method()
    validate(need(m$ok, m$msg))
    d <- m$df
    if (d$time[1] > 0) d <- rbind(transform(d[1, ], time = 0), d)
    par(mar = c(4, 4, 1, 4.5))
    plot(d$time, d$B, type = "n", ylim = c(0, 100), las = 1,
         xlab = "Time (min)", ylab = "%B")
    polygon(c(d$time, rev(d$time)), c(d$B, rep(0, nrow(d))),
            col = adjustcolor("#3498DB", 0.25), border = NA)
    lines(d$time, d$B, lwd = 2, col = "#2C3E50")
    par(new = TRUE)
    fmax <- max(d$flow) * 1.25
    plot(d$time, d$flow, type = "l", lty = 2, lwd = 2, col = "#E67E22",
         axes = FALSE, xlab = "", ylab = "", ylim = c(0, if (fmax > 0) fmax else 1))
    axis(4, las = 1)
    mtext("Flow (mL/min)", side = 4, line = 3)
    legend("topleft", c("%B", "Flow"), col = c("#2C3E50", "#E67E22"),
           lty = c(1, 2), lwd = 2, bty = "n", cex = 0.9)
  })
}

shinyApp(ui, server)
