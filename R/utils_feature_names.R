### ------------------------ Feature Name Helpers ------------------------ ###
#' Check that feature names can be encoded
#'
#' Aborts if any name is missing or empty, if names are duplicated, or if a
#' feature is named "class", which is reserved for the outcome column.
#'
#' @param x Character vector of feature names.
#' @param call Environment of the calling function, used in error messages.
#' @return Invisibly TRUE if all names are valid, otherwise aborts with an error
#'   of class `tentacles_invalid_feature_names`.
#' @keywords internal
validate_feature_names <- function(x, call = parent.frame()) {
  missing_pos <- which(is.na(x) | x == "")
  if (length(missing_pos) > 0) {
    cli::cli_abort(
      c("Feature names must not be missing or empty.",
        "x" = "Found at {cli::qty(length(missing_pos))}position{?s} {missing_pos}."),
      class = "tentacles_invalid_feature_names", call = call
    )
  }

  dups <- unique(x[duplicated(x)])
  if (length(dups) > 0) {
    cli::cli_abort(
      c("Feature names must be unique.",
        "x" = "Duplicated name{?s}: {.val {dups}}."),
      class = "tentacles_invalid_feature_names", call = call
    )
  }

  if ("class" %in% x) {
    cli::cli_abort(
      c("No feature can be named {.val class}.",
        "i" = "{.val class} is reserved for the outcome column."),
      class = "tentacles_invalid_feature_names", call = call
    )
  }

  invisible(TRUE)
}

#' Encode feature names as syntactic ASCII names
#'
#' Names such as "HLA-B" or "C18:1" break formula-based selectors and models,
#' and some engines (C5.0, LightGBM) reject non-ASCII or special characters.
#' Non-ASCII bytes and invalid characters become ".", and collisions get the
#' suffixes ".1", ".2", ... Names that are already syntactic and ASCII are
#' kept as they are and take priority in collisions.
#'
#' @param x Character vector of original feature names.
#' @param call Environment of the calling function, used in error messages.
#' @return A named character vector mapping the internal names (names) to the
#'   original names (values), in the order of `x`.
#' @keywords internal
encode_feature_names <- function(x, call = parent.frame()) {
  validate_feature_names(x, call = call)

  internal <- make.names(iconv(enc2utf8(x), "UTF-8", "ASCII", sub = "."))
  # Unchanged names come first, so make.unique() only suffixes changed ones
  o <- order(internal != x)
  internal[o] <- make.unique(internal[o])

  stats::setNames(x, internal)
}

#' Decode internal feature names back to the original names
#'
#' @param x Character vector of internal names.
#' @param map Named character vector returned by `encode_feature_names()`.
#' @return `x` with internal names replaced by the original names. Names not
#'   found in `map` (e.g. "class") are left unchanged.
#' @keywords internal
decode_feature_names <- function(x, map) {
  idx <- match(x, names(map))
  found <- !is.na(idx)
  x[found] <- unname(map[idx[found]])
  x
}
