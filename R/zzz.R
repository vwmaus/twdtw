.onLoad <- function(libname, pkgname) {
  # Prevent NOTEs from R CMD check
  requireNamespace("proxy", quietly = TRUE)

  # Register twdtw distance function
  if("twdtw" %in% names(proxy::pr_DB$get_entries())) {
    proxy::pr_DB$delete_entry("twdtw")
  }
  proxy::pr_DB$set_entry(FUN = function(x, y, ...) twdtw(x, y, ...), names = c("twdtw"), distance = TRUE, loop = FALSE)

  invisible(NULL)

}
