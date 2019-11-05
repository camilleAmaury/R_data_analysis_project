kmeansfunction <- function(M, k=-1, d="euclidian", initialization="random", precision=0.001, logs=FALSE, seed = -1, kMax = -1){
  if(seed != -1){
    set.seed(seed)
  }

  if(k > dim(M)[1]){
    stop(paste0("Impossible to create ", k, " clusters since number of individuals is ", dim(M)[1]))
  }

  selectBestKvalue <- function(M, kMax = 2, plot = F){
    silhouette_imp <- function(M, kMax, plot){
      matches <- function(i, l){
        res <- c()
        for(j in 1:length(l)){
          if(i == l[j]){
            res <- c(res, j)
          }
        }
        return(res)
      }
      silhouette_ind <- function(i, M, clusters){
        index_my_cluster = matches(clusters[i], clusters)
        a <- 0

        if(length(index_my_cluster)>1){
          for(k in 1:length(index_my_cluster)){
            if(i != index_my_cluster[k]){
              distance <- sqrt(sum((M[index_my_cluster[k],]-M[i,])**2))
              a <- a + distance
            }
          }
          a <- a / (length(index_my_cluster)-1)
        }else{
          a <- 0
        }

        #print(paste0("Distance moyenne du point n",i," à son groupe :", a))

        b <- c()
        for(j in 1:max(clusters)){
          if(j != clusters[i]){
            index_cluster_j = matches(j, clusters)
            bp <- 0
            for(k in 1:length(index_cluster_j)){
              if(i != index_cluster_j[k]){
                bp <- bp + sqrt(sum((M[index_cluster_j[k],]-M[i,])**2))
              }
            }
            bp <- bp / length(index_cluster_j)
            #print(paste0("Distance moyenne du point n",i," au groupe ", j, " :", bp))
            b <- c(b, bp)
          }
        }


        return((min(b)-a)/max(a, min(b)))
      }

      res <- c()

      for(km in 2:kMax){
        kme <- kmeansfunction(M, km, precision = 0.01, seed = 200)
        ssil <- 0
        for(k in 1:km){

          ssili <- 0
          Ik <- matches(k, kme$clusters)
          for(ik in Ik){
            ssili <- ssili + silhouette_ind(ik, M, kme$clusters)
          }
          ssili <- ssili/length(Ik)
          ssil <- ssil + ssili
        }
        ssil <- ssil/km
        res <- c(res, ssil)
      }

      return(res)
    }

    argmax <- function(l){
      max <- 0
      arg <- 1
      for(i in 1:length(l)){
        if(l[i] > max){
          max <- l[i]
          arg <- i
        }
      }
      return(arg)
    }

    softmax <- function(l){
      theta <- c()
      print(l)
      for(i in l){
        theta <- c(theta, exp(i)/sum(exp(l)))
      }
      return(theta)
    }

    res <- silhouette_imp(M, kMax, plot)
    i <- c()
    j <- 2
    for(val in res[2:kMax-1]){
      i <- c(i, (1-(j/nrow(M)))*val)
      j <- j + 1
    }
    return(argmax(softmax(i))+1)
  }

  # k_means for a k fixed
  kmeans_process <- function(M, k=-1, d="euclidian", initialization="random", precision=0.001, logs=FALSE){
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

        return(sqrt(sum((v-u)**2)) <= precision)
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

    return(list("clusters"=new_clusters, "centers"=new_centers[["coords"]], "k"=k))
  }

  # k = -1 to test all cluster number and choose the best cluster number --> to implement
  choosenK <- list()
  if(k == -1){
    if(kMax == -1){
      kMax = min(c(10, nrow(M)))
    }
    return(kmeansfunction(M, selectBestKvalue(M, kMax, logs), d, initialization, precision, logs))
    # choose the best k by elbow method

  }else{
    choosenK <- kmeans_process(M, k, d, initialization, precision, logs)
  }

  M <- cbind(M, choosenK$clusters)

  return(choosenK)
}



generate_dataset <- function(nbCluster, dim=2, nbIndividu=1000, distanceToKernel=1){
  borne_inf <- 0
  borne_sup <- 10
  c <- c()
  for(i in 1:nbCluster){
    cluster <- c(sample(borne_sup, dim, replace=FALSE))
    c <- c(c, cluster)
    for(j in 1:nbIndividu){
      c <- c(c, cluster[1]+(distanceToKernel*runif(1, -1, 1)), cluster[2]+(distanceToKernel*runif(1, -1, 1)))
    }
  }
  data <- matrix(c, ncol = dim, byrow = T)
  plot(data)
  return(data)
}



# --> our algorithm is either implemented with the Forgy or Lloyd algorithm of the stats kmeans function
# --> https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kmeans.html kmeans documentation, for other algorithm
#library("stats")
#print(kmeans(M, 2, algorithm = "Forgy")$cluster)
