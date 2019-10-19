kmeansfunction <- function(M, k=-1, d="euclidian", initialization="random", precision=0.001, logs=FALSE){
  # k = -1 to test all cluster number and choose the best cluster number --> to implement


  if(k > dim(M)[1]){
    stop(paste0("Impossible to create ", k, " clusters since number of individuals is ", dim(M)[1]))
  }

  # Function which determines the first iteration centers
  kcenters_initial <- function(M, k, initialization){
    centers = list()
    if(initialization=="random"){
      indice <- sample(1:dim(M)[1], k, replace=F)
      centers <- list("individuals"=indice, "coords"=M[indice,])
    }else if(initialization=="kmeans++"){
      # choose a random s1 point
      Di <- c(sample(1:dim(M)[1], 1))
      # choose the next points
      for(l in 2:k){
        # compute euclidian distance for every other point, to choose s_l
        M_step <- 1:dim(M)[1]
        M_step <- M_step[M_step != Di]
        # foreach point not already selected, compute the probability to choose the point to the next center point
        Pi <- c()
        for(i in M_step){
          dist_compare <- c()
          # foreach different center existing point
          for(j in Di){
            dist_compare <- c(dist_compare, sqrt(sum((M[j,] - M[i,])**2)))
          }
          # get the minimum distance
          Pi <- c(Pi, min(dist_compare))
        }
        # uniform probability
        proba <- runif(1)
        # cumultaive probability of all probabilities computed
        sum_Pi <- sum(Pi)
        for(i in 2:length(Pi)){
          Pi[i] <- Pi[i-1] + Pi[i]
        }
        Pi <- Pi / sum_Pi
        # get the point regarding to probabilities
        for(i in 1:length(Pi)){
          if(proba <= Pi[i]){
            Di <- c(Di, i)
            break
          }
        }
      }
      centers <- list("individuals"=Di, "coords"=M[Di,])
    }
    return(centers)
  }

  # Function which determines the next iteration centers
  kcenters_next <- function(M, k, clusters){
    gravity_center_vector <- function(M){
      return(if(is.matrix(M)) colMeans(M, na.rm = FALSE, dims = 1) else M)
    }
    # process
    centers = matrix(0, nrow=k, ncol=dim(M)[2])
    keys <- names(clusters)
    for(i in 1:length(keys)){
      cluster <- M[clusters[[keys[i]]],]
      centers[i,] <- gravity_center_vector(cluster)
    }
    return(list("coords"=centers))
  }

  # Function which assign membership of a k-cluster for every individuals
  kclusters <- function(M, distance, centers){
    clusters <- list()
    if(distance == "euclidian"){
      for(i in 1:dim(M)[1]){
        # computing the nearest center
        bestCenter = 1
        bestDist <- 1000000
        for(j in 1:dim(centers[["coords"]])[1]){
          centerPoint <- centers[["coords"]][j,]
          dist <- sqrt(sum((centerPoint - M[i,])**2))
          if(dist < bestDist){
            bestCenter <- j
            bestDist <- dist
          }
        }
        # assigning the individual
        clusters[[toString(bestCenter)]] <- c(clusters[[toString(bestCenter)]], i)
      }
    }else if(distance == "mahalanobis"){
      S <- solve(cov(M))
      for(i in 1:dim(M)[1]){
        # computing the nearest center
        bestCenter = 1
        bestDist <- 1000000
        for(j in 1:dim(centers[["coords"]])[1]){
          centerPoint <- centers[["coords"]][j,]
          dist <- sqrt(t(centerPoint - M[i,]) %*% S %*% (centerPoint - M[i,]))
          if(dist < bestDist){
            bestCenter <- j
            bestDist <- dist
          }
        }
        # assigning the individual
        clusters[[toString(bestCenter)]] <- c(clusters[[toString(bestCenter)]], i)
      }
    }
    return(clusters)
  }

  # Function which determines if two matrix are equal, regarding to their rows and a specified precision rate
  equal_matrix_round <- function(U, V, precision){
    # Function which determines if two vectors are equal, regarding to a specified precision rate
    equal_vectors_round <- function(u, v, precision){
      if(length(u) != length(v)){
        stop(paste0("Impossible to compare a vector of length ", length(u), " to a vector of length ", length(v)))
      }
      equal = TRUE
      for(i in 1:length(u)){
        cond1 <- u[i] <= v[i] + (if(u[i] >= 0) precision else (-precision))
        cond2 <- u[i] >= v[i] - (if(u[i] >= 0) precision else (- precision))
        if(!(cond1 && cond2)){
          equal = FALSE
          break
        }
      }
      return(equal)
    }
    if(dim(U)[1] != dim(V)[1]){
      stop(paste0("Impossible to compare a matrix of ", dim(U)[1], " rows to a matrix of ", dim(V)[1], " rows"))
    }
    equal = TRUE
    for(i in 1:dim(U)[1]){
      if(!equal_vectors_round(U[i,], V[i,], precision)){
        equal = FALSE
        break
      }
    }
    return(equal)
  }

  # Process
  last_centers <- kcenters_initial(M, k, initialization)
  clusters <- kclusters(M, d, last_centers)
  if(logs){
    print("initial centers : ")
    print(last_centers[["coords"]])
    print("clusters [individuals of M] : ")
    print(clusters)
  }
  new_centers = list("coords"=matrix(0, nrow=k, ncol=dim(M)[2]))
  step <- 1
  # while our precision changes between two iterations is higher than precision, doing another step
  while(!equal_matrix_round(new_centers[["coords"]], last_centers[["coords"]], precision)){
    if(step == 1){
      new_centers <- kcenters_next(M, k, clusters)
    }else{
      last_centers <- new_centers
      new_centers <- kcenters_next(M, k, clusters)
    }
    clusters <- kclusters(M, d, new_centers)
    if(logs){
      print(paste0("\n # ---------------> step ", toString(step), " <--------------- # \n"))
      print("actual centers : ")
      print(new_centers[["coords"]])
      print("clusters [individuals of M] : ")
      print(clusters)
    }
    step <- step + 1
  }

  # reformat cluster in a unique vector
  new_clusters <- 1:dim(M)[1]
  for(i in 1:length(names(clusters))){
    for(j in 1:length(clusters[[names(clusters)[i]]])){
      new_clusters[clusters[[names(clusters)[i]]][j]] = as.integer(names(clusters)[i])
    }
  }

  return(list("clusters"=new_clusters, "centers"=new_centers[["coords"]]))
}
M <- matrix(c(1:3, c(7,3,5), c(5,11,30), c(20,15,1), c(20,40,19)), nrow=5, byrow=T)
print(M)
print(kmeansfunction(M, k=2, initialization = "kmeans++"))
print(kmeansfunction(M, k=2))

# --> our algorithm is either implemented with the Forgy or Lloyd algorithm of the stats kmeans function
# --> https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kmeans.html kmeans documentation, for other algorithm
#library("stats")
#print(kmeans(M, 2, algorithm = "Forgy")$cluster)
