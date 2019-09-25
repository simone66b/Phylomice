source("kvals.R")
setup <- list(
  miceMethod = c("phnorm"),
  treeSize = c(64),
  traitSize = c(5),
  treeShape = c("pureBirth"),
  lambdaTransform = c(1),
  missDataMech = c("mar"),
  missDataLevel = c(3)
)
print(system.time({
	runSimulations(setup)
}))
