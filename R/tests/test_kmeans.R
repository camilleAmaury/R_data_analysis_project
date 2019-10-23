library(NantesStatisticalAnalysis)

test_kmeans <- function(){
  M <- matrix(c(1:3, c(7,3,5), c(5,11,30), c(20,15,1), c(20,40,19)), nrow=5, byrow=T)
  print(M)
  print(kmeansfunction(M, k=2, initialization = "kmeans++"))
  print(kmeansfunction(M, k=2))
}

test_kmeans()
