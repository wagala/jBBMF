#===============================================================================
# Create an optional progress bar
#===============================================================================

# Creates a progress bar when the progress package is available.
# Returns NULL otherwise, allowing the sampler to run without displaying
# progress information.

.progress_bar_or_null <- function(total) {
  if (requireNamespace("progress", quietly = TRUE)) {
    progress::progress_bar$new(
      format = "working [:bar] :percent in :elapsed",
      total = total,
      clear = FALSE,
      width = 80
    )
  } else {
    NULL
  }
}