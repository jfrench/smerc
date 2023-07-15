if (requireNamespace("sf", quietly = TRUE) &
    requireNamespace("sp", quietly = TRUE)) {
  nysp <- sf::as_Spatial(get0("smerc_nysf_internal",
                              envir = asNamespace("smerc")))
} else {
  message("The sf and sp packages must be installed to produce this data. Returning a NULL object.")
  nysp <- NULL
}
