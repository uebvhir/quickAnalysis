
context("quickAnalysus test")
test_that("quickAnalysis exports", {

  ex <- quickAnalysis(expres = exprs(sampleSet),
                      groupingVar = pData(sampleSet)$group,
                      min = "Case", sust = "Control", pvalThreshold = 0.25)

  expect_that(ex$genes, is_a("character"))
  expect_that(ex$resT, is_a("data.frame"))
  expect_that(ex$topTable, is_a("data.frame"))
})


test_that("quickAnalysis top Tables", {

  ex <- quickAnalysis(expres = exprs(sampleSet),
                      groupingVar = pData(sampleSet)$group,
                      min = "Case", sust = "Control", pvalThreshold = 0.25)

  expect_equivalent(sign(ex$resT$teststat)[1:5], sign(ex$topTable$t)[1:5])
  expect_equal(sort(rownames(ex$resT)), sort(rownames(ex$topTable)))
})

