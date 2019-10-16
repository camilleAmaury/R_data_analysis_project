
pca <- function(M, center=TRUE, scale=TRUE, bias=FALSE, AxisNb=2, Q=diag(dim(M)[2]), D=(1/dim(M)[1])diag(dim(M)[1])){
  # returned object :
  base_result = list("base"=M)

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
    M_temp <- t(Mc[["M_centered"]]) %% Mc[["M_centered"]]
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
      return(list("M_scale"= M_centered %% D1_S))
    }
    result <- append(result, centered_reduced_matrix(M_centered, result[["D1_div_sd"]]))

    # calcul de la matrice de correlation
    R_temp <- t(result[["M_scale"]]) %% result[["M_scale"]]
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
  # component PCA
  component_pca <- function(X, Q, D, AxisNb){
    # building PCA main component
    pca_component <- list("Si"=t(X) %% D %% X %% Q)
    eig <- eigen(pca_component[["Si"]])
    eig[["vectors"]] <- round(eig[["vectors"]], 6)
    eig[["values"]] <- round(eig[["values"]], 6)
    pca_component <- append(pca_component, eig)
    # calculs des Fi et Gi
    individuals_variables <- function(X, eigenValues, eigenVectors){
      return(list("Fi"=X%%eigenVectors, "Gi"=sqrt(eigenValues)eigenVectors))
    }
    pca_component <- append(pca_component,
                            individuals_variables((if(scale) base_result[["M_scale"]] else base_result["M_centered"]),
                                                  pca_component[["values"]], pca_component[["vectors"]])
    )
    # Valeurs propres et vecteurs propres (axes principaux)
    pca_component <- append(pca_component, eigen(pca_component[["Si"]]))
    # Facteurs principaux Uk
    pca_component <- append(pca_component, list("PrincipalFactors"=Q%%pca_component$vectors))
    return(pca_component)
  }
  X <- if(scale) base_result[["M_scale"]] else (if(center) base_result[["M_centered"]] else M)
  pca_res <- component_pca(X, Q, D, AxisNb)


  li_res <- list("base_component"=base_result, "PCA_component"=pca_res)
  return(li_res)
}



inertia <- function(eigenvalues){
  return(sum(eigenvalues))
}

cumulative_information_caught <- function(eigenvalues, num_axes = 2, percent = T){
  return(sum(eigenvalues[1:num_axes])/inertia(eigenvalues)((1-(percent1))+(percent1)100))
}

information_caught <- function(eigenvalues, index_axe, percent = T){
  return(eigen_values[index_axe]/inertia(eigenvalues)((1-(percent1))+(percent1)100))
}
