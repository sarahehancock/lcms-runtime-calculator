# LC-MS Run Time & Cost Calculator

A browser-based Shiny app (via [Shinylive](https://posit-dev.github.io/r-shinylive/)) that estimates total LC-MS sequence time and instrument cost from:

- number of samples and technical replicates
- gradient time and per-injection buffer time
- total blanks and QCs (start, interval, end), with optional separate method times
- single polarity, separate pos & neg injections, or polarity switching
- instrument rate ($/h), booking fee, and optional whole-hour billing

It also generates the full injection order with estimated start times, downloadable as CSV.

The app runs entirely in the browser (webR) — no server needed.

## Repository layout

```
app/app.R                       # the Shiny app
.github/workflows/deploy.yml    # builds with shinylive and deploys to GitHub Pages
```

## Deploying to GitHub Pages

1. Push this repository to GitHub (default branch `main`).
2. In the repo, go to **Settings → Pages** and set **Source** to **GitHub Actions**.
3. Push a commit (or run the workflow manually from the **Actions** tab).
4. The app will be live at `https://<username>.github.io/<repo-name>/`.

The first load takes ~10–30 s while webR and packages download; later loads are cached by the browser.

## Running locally

```r
# As a normal Shiny app
shiny::runApp("app")

# Or build and preview the Shinylive version
install.packages(c("shinylive", "httpuv"))
shinylive::export("app", "site")
httpuv::runStaticServer("site")
```

Opening `site/index.html` directly from disk will not work — it must be served over HTTP.
