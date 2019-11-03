M <- t(read.csv(file="Projet M1 AD 1920.csv", header=FALSE, sep=","))

# our function
acp_package <- pcafunction(M, axisNumber=2)

#explain function
explainPCA(acp_package)
