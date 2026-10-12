is_ascii <- function(x) !grepl("[^ -~]", x)

test_that("encode_feature_names maps syntactic ASCII names to the originals", {
  map <- encode_feature_names(tricky_names)

  expect_type(map, "character")
  expect_equal(unname(map), unname(tricky_names))
  expect_equal(make.names(names(map)), names(map))
  expect_false(anyDuplicated(names(map)) > 0)
  expect_true(all(is_ascii(names(map))))
})

test_that("encode_feature_names keeps syntactic names and suffixes the changed ones", {
  map <- encode_feature_names(c("HLA-B", "HLA.B", "C18:1", "C18.1", "A2M"))

  expect_equal(names(map), c("HLA.B.1", "HLA.B", "C18.1.1", "C18.1", "A2M"))
})

test_that("encode_feature_names replaces non-ASCII characters", {
  # Greek alpha and beta, built from code points to keep this file ASCII.
  # Each takes two bytes in UTF-8, so each becomes two dots, plus one for "-".
  x <- paste0(intToUtf8(c(945, 946), multiple = TRUE), "-tocopherol")
  map <- encode_feature_names(x)

  expect_equal(names(map), c("...tocopherol", "...tocopherol.1"))
  expect_equal(unname(map), x)
})

test_that("encode_feature_names handles latin1 and invalid UTF-8 strings", {
  # "cafe" with an e-acute, marked as latin1
  latin1 <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
  Encoding(latin1) <- "latin1"
  # A byte that is not valid UTF-8, with unknown encoding (e.g. a file read
  # without the right `encoding`)
  invalid <- rawToChar(as.raw(c(0x67, 0xff, 0x31)))

  map <- encode_feature_names(c(latin1, invalid))

  expect_false(anyNA(names(map)))
  expect_equal(make.names(names(map)), names(map))
  expect_true(all(is_ascii(names(map))))
  expect_identical(decode_feature_names(names(map), map), c(latin1, invalid))
})

test_that("decode_feature_names restores originals and leaves other names unchanged", {
  map <- encode_feature_names(c("HLA-B", "HLA.B", "C18:1"))

  expect_equal(
    decode_feature_names(c("C18.1", "HLA.B.1", "class", "HLA.B"), map),
    c("C18:1", "HLA-B", "class", "HLA.B")
  )
})

test_that("encode_feature_names aborts on missing or empty names", {
  expect_error(encode_feature_names(c("A2M", NA)), class = "tentacles_invalid_feature_names")
  expect_error(encode_feature_names(c("A2M", "")), class = "tentacles_invalid_feature_names")
  expect_error(
    encode_feature_names(c("A2M", NA, "")),
    regexp = "positions 2 and 3",
    class = "tentacles_invalid_feature_names"
  )
})

test_that("encode_feature_names aborts on duplicated names, naming them", {
  expect_error(
    encode_feature_names(c("HLA-B", "A2M", "HLA-B")),
    regexp = "HLA-B",
    class = "tentacles_invalid_feature_names"
  )
})

test_that("encode_feature_names aborts when a feature is named 'class'", {
  expect_error(
    encode_feature_names(c("A2M", "class")),
    class = "tentacles_invalid_feature_names"
  )
})
