# Ignore checkout line endings while retaining every other source byte.
canonical_source_md5 <- function(path) {
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  if (length(bytes) > 1L) {
    crlf <- c(bytes[-length(bytes)] == as.raw(13L) &
                bytes[-1L] == as.raw(10L), FALSE)
    bytes <- bytes[!crlf]
  }
  canonical <- tempfile()
  on.exit(unlink(canonical), add = TRUE)
  writeBin(bytes, canonical)
  unname(tools::md5sum(canonical))
}
