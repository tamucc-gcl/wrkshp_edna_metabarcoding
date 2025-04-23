# ── Begin renv autoloading ───────────────────────────────────────────────────────
source("renv/activate.R")
# ── End renv autoloading ─────────────────────────────────────────────────────────

if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  setHook("rstudio.sessionInit", function(newSession) {
    if (newSession && rstudioapi::isAvailable()) {
      later::later(function() {
        files_to_open <- c(
          "analysis/data_cleanup.qmd"
        )
        for (f in files_to_open) {
          if (file.exists(f)) {
            try(rstudioapi::navigateToFile(f, line = -1L, column = -1L), silent = TRUE)
          }
        }
      }, delay = 1) # delay 1 second
    }
  }, action = "append")
}

