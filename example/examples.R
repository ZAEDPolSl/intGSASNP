test <- ORAtest(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = TRUE,
  adjust.method = "BH"
)

test <- CERNOtest(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = FALSE,
  adjust.method = "BH"
)

test <- MAGENTAtest(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = FALSE,
  cutoff = 0.05,
  adjust.method = "BH",
  permutations = 100
)

test <- GSEASNPtest(
  data = example,
  pathways = KEGGhsa,
  adjust.method = "BH",
  permutations = 100
)

test <- GSEAtest_serial(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = FALSE,
  adjust.method = "BH",
  permutations = 10,
  permutation.method = "entrez"
)

test <- GSEAtest_parallel(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = FALSE,
  adjust.method = "BH",
  permutations = 10,
  permutation.method = "entrez"
)

test <- iGSEA4GWAStest_serial(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = FALSE,
  adjust.method = "BH",
  permutations = 10
)

test <- iGSEA4GWAStest_parallel(
  data = example,
  pathways = KEGGhsa,
  int.method = "min",
  adjust.int = FALSE,
  int.cor.LD = FALSE,
  adjust.method = "BH",
  permutations = 10
)
