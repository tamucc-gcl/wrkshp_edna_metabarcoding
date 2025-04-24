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

# only run on Colab
if (nzchar(Sys.getenv("COLAB"))) {
  # paths under Drive
  gdrive <- "/content/edna_libraries/library"
  lib    <- file.path(gdrive, "library")
  cache  <- file.path(gdrive, "cache")

  dir.create(lib,   recursive = TRUE, showWarnings = FALSE)
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)

  # tell renv to use Drive
  Sys.setenv(
    RENV_PATHS_LIBRARY = lib,
    RENV_PATHS_CACHE   = cache
  )
  .libPaths(lib)
}
