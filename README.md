# LC-MS Run Time, Cost & Solvent Calculator

A browser-based Shiny app for planning LC-MS sequences. It estimates total run time, finish time, instrument cost and solvent use from your samples, QCs, blanks, polarity and LC method, and generates the injection order as a CSV.

**Open the app:** https://sarahehancock.github.io/lcms-runtime-calculator/

The app runs entirely in your browser via [Shinylive](https://posit-dev.github.io/r-shinylive/) (webR) — no server, no installation, and nothing you enter leaves your computer. The first load takes ~10–30 s while R downloads; later visits are faster.

---

## How to use it

All settings are in the sidebar on the left, grouped into collapsible sections. Results update as you type. The app has two tabs along the top: **Run time & cost** and **Solvents**.

### 1. Samples & method

| Setting | What it does |
|---|---|
| **LC method** | Choose a built-in method, one you've saved, or **Manual entry** (see below). |
| **Set run time from LC method** | When ticked, the run time is taken from the end of the method's timetable. Untick to type your own. |
| **Number of samples** | Biological samples in the batch. |
| **Injections per sample** | Technical replicates (e.g. 2 = each sample injected twice). |
| **Gradient / method run time** | Length of one injection's LC method, in minutes. |
| **Buffer per injection** | Extra time per injection for autosampler draw, needle wash or any dead time between runs. |

**Manual entry:** to use a method that isn't in the list, choose *Manual entry*. This clears the **Method editor** on the **Solvents** tab (the sidebar has a link to it), where you type the pump timetable — see [Adding or editing a method](#adding-or-editing-a-method). Results update as you type; save it with a name to add it to the LC method list.

### 2. Polarity

- **Single polarity** — one injection per sample (pos or neg).
- **Pos & neg (separate injections)** — every sample, QC and blank is injected once in each polarity. Choose whether to run all pos then all neg (with optional source switchover/stabilisation time between them) or to alternate pos/neg per injection.
- **Polarity switching** — one injection per sample, both polarities acquired in the same run.

### 3. Blanks & QCs

- **Total blanks** — spread evenly through the sequence, with the first at the start and the last at the end. Optionally give blanks a shorter method time.
- **QCs at start** — conditioning QCs before the first sample.
- **QC every N samples** — interval QCs (0 = none).
- **QCs at end** — closing QCs after the last sample.

Counts are per polarity, so with separate pos & neg injections they are doubled.

### 4. Solvents

| Setting | What it does |
|---|---|
| **Safety margin** | Extra % added to the calculated volumes. |
| **Prime / purge volume** | Volume per channel used to prime the pump before the run. |
| **Bottle size** | Used to report how many bottles of A and B you need. |
| **Count flow during buffer and switchover** | Includes solvent pumped between injections, at the method's starting %B and flow. |

### 5. Column equilibration

Time to equilibrate the column before the first injection. It's added to the start of the sequence and counted in run time, cost and solvent use.

- **Column volumes** (default) — enter column length and internal diameter, porosity, and the number of column volumes (20–50 is typical; HILIC usually needs the high end). The app shows the column volume and the resulting time at the method's starting flow.
  Column volume = π × (ID/2)² × length × porosity. Use porosity ~0.65 for fully porous particles and ~0.5 for core-shell.
- **Fixed time** — enter minutes directly.
- **None** — no pre-sequence equilibration.

Each built-in or saved method fills in its own column settings when selected.

### 6. Instrument cost

Enter the hourly instrument rate and any fixed booking fee. Tick **Bill in whole hours** to round up to the next hour. The cost box shows the total, cost per sample and hours billed.

### 7. Start time

Set the start date and time (24 h, e.g. `09:00`) to get an estimated finish time and a start time for each injection.

---

## Reading the results

### Run time & cost tab

- **Summary boxes** — total run time, total injections, estimated finish and instrument cost.
- **Breakdown** — time (and cost) split between samples, QCs, blanks, polarity switchover and equilibration.
- **Injection sequence** — the full run order with polarity, sample names and estimated start times. Click **CSV** to download it as a starting point for your instrument sequence.

### Solvents tab

- **Solvent A / Solvent B needed** — total volume to prepare, including margin and priming, plus the number of bottles.
- **Per injection** — volume of A and B used by one run of the gradient.
- **Solvent breakdown** — volume per injection, all injections, equilibration, buffer/switchover, safety margin and priming.
- **Gradient** — plot of %B (solid) and flow (dashed) over the method.

Every injection, including blanks and QCs, is assumed to run the selected gradient.

---

## LC methods

### Built-in methods

| Method | Column | Flow | Source |
|---|---|---|---|
| Lipidomics (RP, 30 min) | ACQUITY UPLC CSH C18, 2.1 × 100 mm, 1.7 µm; 20 column volumes | 0.26 mL/min | Gradient from Hu et al. (2008) *J Proteome Res* 7:4982–4991, [doi:10.1021/pr800373m](https://doi.org/10.1021/pr800373m) |
| Metabolomics (Z-HILIC, Girel 2023) | Atlantis Premier BEH Z-HILIC, 2.1 × 150 mm, 1.7 µm; 50 column volumes | 0.35 mL/min gradient, 0.45 mL/min equilibration | Girel et al. (2023) *J Chromatogr A* 1697:463994, [doi:10.1016/j.chroma.2023.463994](https://doi.org/10.1016/j.chroma.2023.463994) |

Built-in methods can't be overwritten or deleted, but you can edit one and save it under a new name.

### Adding or editing a method

The **Method editor** is at the bottom of the **Solvents** tab.

1. Select a method to start from, or choose **Manual entry** to start with an empty timetable.
2. Enter the pump timetable, one row per line: `time (min), %B, flow (mL/min)`. For example:
   ```
   0, 32, 0.26
   1.5, 32, 0.26
   4, 45, 0.26
   ...
   25, 97, 0.26
   25.1, 32, 0.26
   30, 32, 0.26
   ```
   - %B and flow change **linearly** between rows, as on the pump, so a gradient segment is just two rows.
   - For a **step change** (e.g. returning to starting conditions for re-equilibration), use two rows close together, e.g. `25, 97` then `25.1, 32`.
   - Times before the first row are held at the first row's conditions.
   - Commas, tabs, semicolons or spaces all work as separators. Lines starting with `#` are ignored.
3. Set the column details under **Column equilibration** in the sidebar.
4. Type a **Method name** and click **Save method**. The timetable and column settings are saved together.

Edits are used in the calculations straight away, even before saving; "Unsaved edits" appears under the box until you save.

### Saving and sharing methods

Saved methods are stored **in your browser only** — they won't appear on another computer or in another browser, and clearing your browser data removes them.

- **Export my methods** downloads all your saved methods as a `.json` file. Keep it as a backup or send it to colleagues.
- **Import...** loads methods from an exported `.json` file. Methods with the same name are replaced.

---

## For developers

### Repository layout

```
app/app.R                       # the Shiny app
.github/workflows/deploy.yml    # builds with shinylive and deploys to GitHub Pages
```

### Deploying to GitHub Pages

1. Push this repository to GitHub (default branch `main`).
2. In the repo, go to **Settings → Pages** and set **Source** to **GitHub Actions**.
3. Push a commit (or run the workflow manually from the **Actions** tab).
4. The app will be live at `https://<username>.github.io/<repo-name>/`.

Every push to `main` rebuilds and redeploys the app.

### Running locally

```r
# As a normal Shiny app
shiny::runApp("app")

# Or build and preview the Shinylive version
install.packages(c("shinylive", "httpuv"))
shinylive::export("app", "site")
httpuv::runStaticServer("site")
```

Opening `site/index.html` directly from disk will not work — it must be served over HTTP.

To add or change a built-in method, edit `builtin_methods` and `builtin_columns` near the top of `app/app.R`.
