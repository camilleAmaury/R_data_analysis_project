
pcafunction <- function(M, center=TRUE, scale=TRUE, bias=FALSE, Q=diag(dim(M)[2]), D=(1/dim(M)[1])*diag(dim(M)[1]), axisMethod="kaiser", axisNumber=-1){
  # returned object :
  base_result = list("base"=M)


  # --------------------------------------------------------- BASE COMPUTATION FOR PCA --------------------------------------------------------- #
  # fonction qui calcule la matrice de covariance de la matrice initiale et renvoit toutes les étapes précédentes de calcul
  covariance_matrix <- function(M){
    # fonction qui retourne la matrice des données centrées et le centre de gravité de la matrice initiale
    centered_matrix <- function(M){
      # fonction qui calcule la matrice G, centre de gravité de la matrice passée en paramètre
      gravity_center_matrix <- function(M){
        moy_T <- colMeans(M, na.rm = FALSE, dims = 1)
        M_moy_t <- t(matrix(rep(moy_T, dim(M)[1]), ncol = dim(M)[1]))
        return(M_moy_t)
      }
      G <- gravity_center_matrix(M)
      return(list("G"=G, "M_centered"=M-G))
    }
    Mc <- centered_matrix(M)
    M_temp <- t(Mc[["M_centered"]]) %*% Mc[["M_centered"]]
    S_bias <- (1 / dim(M)[1]) * M_temp
    S <- (1 / (dim(M)[1]-1)) * M_temp
    return(append(Mc, list("S"=S, "S_bias"=S_bias)))
  }
  if(scale || center){
    base_result <- append(base_result, covariance_matrix(M))
  }
  # fonction de calcul de la matrice de correlation de la matrice initiale et renvoit toutes les étapes précédentes de calcul
  correlation_matrix <- function(S, M_centered){
    # fonction qui calcule la matrice diagonale des (1/écarts-type) des colonnes
    diag_ecart_type <- function(S){
      return(list("D1_div_sd"= (1/sqrt(diag(S))) * diag(dim(S)[1])))
    }
    result <- diag_ecart_type(S)

    # fonction qui calcule la matrice centrée réduite
    centered_reduced_matrix <- function(M_centered, D1_S){
      return(list("M_scale"= M_centered %*% D1_S))
    }
    result <- append(result, centered_reduced_matrix(M_centered, result[["D1_div_sd"]]))

    # calcul de la matrice de correlation
    R_temp <- t(result[["M_scale"]]) %*% result[["M_scale"]]
    R <- (1 / dim(M)[1]) * R_temp
    return(append(result, list("R"=R)))
  }
  if(scale){
    base_result <- append(base_result, correlation_matrix(base_result[["S_bias"]], base_result[["M_centered"]]))
  }
  if(bias){
    base_result[["S"]] <- base_result[["S_bias"]]
  }
  base_result[["S_bias"]] <- NULL


  # --------------------------------------------------------- MAIN PCA COMPONENT --------------------------------------------------------- #
  component_pca <- function(X, Q, D, AxisMethod, axisNumber){
    number_axes_conservated <- function(eigenvalues, method){
      i <- 0
      if(method=="kaiser"){
        for(val in eigenvalues){
          if(val > 1){
            i <- i + 1
          }
        }
      }

      if(method=="elbow"){
        lastval <- eigenvalues[1]
        i <- i + 1
        for(val in eigenvalues[2:length(eigenvalues)]){
          if((lastval - val)/lastval <= 0.5){
            lastval <- val
            i <- i + 1
          }
        }
      }

      return(i)
    }
    # building PCA main component
    pca_component <- list("Si"=t(X) %*% D %*% X %*% Q)

    # eigen values and vectors
    eig <- eigen(pca_component[["Si"]])
    eig[["vectors"]] <- round(eig[["vectors"]], 6)
    eig[["values"]] <- round(eig[["values"]], 6)

    # choose the right number of axis
    axis_conservated <- if(axisNumber == -1) number_axes_conservated(eig[["values"]], axisMethod) else axisNumber
    eig[["vectors"]] <- matrix(eig[["vectors"]][,1:axis_conservated], ncol=axis_conservated)
    eig[["values"]] <- eig[["values"]][1:axis_conservated]
    pca_component <- append(pca_component, eig)

    # calculs des Fi et Gi
    individuals_variables <- function(X, axis_conservated, eigenValues, eigenVectors){
      # individuals
      Fi <- matrix(0, nrow=dim(X)[1], ncol=axis_conservated)
      for(i in 1:axis_conservated){
          Fi[,i] <- X[i,] %*% eigenVectors[,i]
      }
      # variables
      Gi <- matrix(0, nrow=dim(eigenVectors)[1], ncol=axis_conservated)
      for(i in 1:axis_conservated){
        Gi[,i] <- sqrt(eigenValues[i])*eigenVectors[,i]
      }
      return(list("Fi"=Fi, "Gi"=Gi))
    }
    pca_component <- append(pca_component,
                            individuals_variables((if(scale) base_result[["M_scale"]] else base_result["M_centered"]),
                                                  axis_conservated, pca_component[["values"]], pca_component[["vectors"]])
    )
    # Valeurs propres et vecteurs propres (axes principaux)
    #pca_component <- append(pca_component, eigen(pca_component[["Si"]]))
    # Facteurs principaux Uk
    #pca_component <- append(pca_component, list("PrincipalFactors"=Q%%pca_component$vectors))
    return(pca_component)
  }
  X <- if(scale) base_result[["M_scale"]] else (if(center) base_result[["M_centered"]] else M)
  pca_res <- component_pca(X, Q, D, AxisMethod, axisNumber)
  li_res <- list("base_component"=base_result, "PCA_component"=pca_res)


  # --------------------------------------------------------- PCA INTERPRETATION RESULT --------------------------------------------------------- #
  pca_inter <- list()
  inertia <- function(eigenvalues){
    return(sum(eigenvalues))
  }

  cumulative_information_caught <- function(eigenvalues, num_axes = 2, percent = T){
    return(sum(eigenvalues[1:num_axes])/inertia(eigenvalues)((1-(percent1))+(percent1)*100))
  }

  information_caught <- function(eigenvalues, index_axe, percent = T){
    return(eigen_values[index_axe]/inertia(eigenvalues)((1-(percent1))+(percent1)*100))
  }


  contribution_abs <- function(individual, main_comp, F, eigenvalues){
    return((F[individual, main_comp]**2)/(length(eigenvalues)*eigenvalues[main_comp]))
  }

  contribution_rel <- function(individual, main_comp, F, eigenvalues){
    return(F[individual, main_comp]**2/rowSums(F**2)[individual])
  }

  cumulative_contribution_abs <- function(individual, main_comps, F, eigenvalues){
    return((F[individual, 1:main_comps]**2)/(length(eigenvalues)*eigenvalues[1:main_comps]))
  }

  cumulative_contribution_rel <- function(individual, main_comps, F, eigenvalues){
    return(rowSums(F[individual, 1:main_comps]**2)/rowSums(F**2)[individual])
  }

  # process here


  li_res <- append(li_res, pca_inter)
  return(li_res)
}

M <- matrix(c(1:3, c(7,3,5), c(5,11,30), c(20,15,1), c(20,40,19)), nrow=5, byrow=T)
print(pcafunction(M))
