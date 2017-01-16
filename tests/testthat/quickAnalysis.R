
context("quickAnalysis test")
test_that("quickAnalysis exports", {

  ex <- quickAnalysis(expres = exprs(sampleSet),
                      groupingVar = pData(sampleSet)$group,
                      min = "Case", sust = "Control", pvalThreshold = 0.25)

  expect_that(ex$genes, is_a("character")) ## l'objecte genes es tracta d'un vector caracter
  expect_that(ex$resT, is_a("data.frame")) ## l'objecte resT, es tracta d'un data.frame
  expect_that(ex$topTable, is_a("data.frame")) ## l'objecte topTable, es tracta d'un data.frame
})


test_that("quickAnalysis top Tables", {

  ex <- quickAnalysis(expres = exprs(sampleSet),
                      groupingVar = pData(sampleSet)$group,
                      min = "Case", sust = "Control", pvalThreshold = 0.25)

  expect_equivalent(sign(ex$resT$teststat)[1:5], sign(ex$topTable$t)[1:5]) ## Surt el mateix top5 amb els dos metodes
  expect_equal(sort(rownames(ex$resT)), sort(rownames(ex$topTable))) ## Hi ha els mateixos gens en ambdos anàlisis
  expect_equal((ex$topTable$meanCase - ex$topTable$meanControl), ex$topTable$logFC  ) ## la resta de mitjanes correspon al logFC
  expect_equal(ex$topTable$P.Value, sort(ex$topTable$P.Value)) ## comprovem que estiguin endreçats per p.valor
})

