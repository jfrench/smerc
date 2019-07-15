.onAttach <- function(...) {
    date <- date()
    x <- regexpr("[0-9]{4}", date)
    yr <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)

	greet <- paste("# This research was partially supported under NSF grant 1463642")
	packageStartupMessage(greet)
}

