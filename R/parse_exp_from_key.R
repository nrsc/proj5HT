parse_exp_from_key <- function(key) {
  m <- sub(".*\\|exp=([^|]*)\\|.*", "\\1", key)
  m <- trimws(m)

  if (m == "" || is.na(m)) {
    "standard puff"
  } else if (grepl("wash-in", m, ignore.case = TRUE)) {
    "wash-in"
  } else if (grepl("wash-out|wout", m, ignore.case = TRUE)) {
    "wash-out"
  } else {
    m
  }
}
