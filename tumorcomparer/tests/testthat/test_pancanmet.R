# All tests are done on files in package using system.file()

context("ID Conversion")

test_that("convertChemIds", {
    results <- convertChemIds("6305", "PubChem CID", "ChEBI")
    expect_identical(results, c("CHEBI:16828", "CHEBI:57912"))
})
