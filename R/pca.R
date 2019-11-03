
pcafunction <- function(M, center=TRUE, scale=TRUE, bias=FALSE, Q=diag(dim(M)[2]), D=(1/dim(M)[1])*diag(dim(M)[1]), axisMethod="elbow", axisNumber=-1){
  # returned object :
  base_result = list("base"=as.matrix(M))
  rownames = c()
  for(i in 1:dim(M)[1]){
    rownames = c(rownames, paste0("i", i))
  }
  base_result = append(base_result, list("rownames"=rownames))
  colnames = c()
  if(is.null(names(M))){
    for(i in 1:dim(M)[2]){
      colnames = c(colnames, paste0("v", i))
    }
  }else{
    colnames = names(M)
  }
  base_result = append(base_result, list("colnames"=colnames))
  M <- as.matrix(M)
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
  component_pca <- function(X, Q, D, AxisMethod, axisNumber, base_result){
    cumulative_information_caught <- function(eigenvalues){
      inertia <- function(eigenvalues){
        return(sum(eigenvalues))
      }
      information_caught <- function(eigenvalues){
        return(eigenvalues/inertia(eigenvalues)*100)
      }
      res <- list("inertia"=inertia(eigenvalues),"axisIC"=information_caught(eigenvalues))
      return(append(res, list("cumulative_axisIC"=cumsum(res$axisIC))))
    }
    number_axes_conservated <- function(inertia, global_inertia, method, size=length(inertia)){
      res <- c()
      if(method=="kaiser"){
        res <- inertia[inertia > 100/global_inertia]
      }

      if(method=="elbow"){
        res <- inertia[inertia > 100/size]
      }
      return(length(res))
    }
    # building PCA main component
    pca_component <- list("Si"=t(X) %*% D %*% X %*% Q, "X"=X)

    # eigen values and vectors
    eig <- eigen(pca_component$Si)
    eig[["values"]] <- round(eig$values[round(eig$values, 5) != 0], 6)
    eig[["vectors"]] <- round(eig$vectors[,1:length(eig$values)], 6)

    # rank
    pca_component[["rank"]] <- length(eig$values)

    # choose the right number of axis
    pca_component <- append(pca_component,cumulative_information_caught(eig$values))
    axis_conservated <- if(axisNumber == -1) number_axes_conservated(pca_component$axisIC, pca_component$inertia, axisMethod) else axisNumber
    pca_component <- append(pca_component, append(list("axisNumber"=axis_conservated),eig))

    # calculs des Fi et Gi
    individuals_variables <- function(X, eigenValues, eigenVectors){
      # individuals / principal component
      Fi <- matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
      for(i in 1:length(eigenValues)){
          Fi[,i] <- X %*% eigenVectors[,i]
      }
      # variables
      Gi <- matrix(0, nrow=dim(eigenVectors)[1], ncol=dim(X)[2])
      for(i in 1:length(eigenValues)){
        Gi[,i] <- sqrt(eigenValues[i])*eigenVectors[,i]
      }
      return(list("Fi"=Fi, "Gi"=Gi))
    }
    pca_component <- append(pca_component,
                            individuals_variables((if(scale) base_result$M_scale else base_result$M_centered), pca_component$values, pca_component$vectors)
    )
    # Valeurs propres et vecteurs propres (axes principaux)
    #pca_component <- append(pca_component, eigen(pca_component[["Si"]]))
    # Facteurs principaux Uk
    #pca_component <- append(pca_component, list("PrincipalFactors"=Q%%pca_component$vectors))
    return(pca_component)
  }
  X <- if(scale) base_result$M_scale else (if(center) base_result$M_centered else M)
  pca_res <- component_pca(X, Q, D, AxisMethod, axisNumber, base_result)
  li_res <- list("base_component"=base_result, "PCA_component"=pca_res)


  # --------------------------------------------------------- PCA INTERPRETATION RESULT --------------------------------------------------------- #
  interpretation_pca <- function(Fi, Gi, Q, D, eigenvalues, axisNumber){
    pca_inter <- list()

    contribution_rel_individuals <- function(Fi, axisNumber){
      qlt <- matrix(0, nrow=dim(Fi)[1], ncol=axisNumber)
      for(column in 1:axisNumber){
        for(row in 1:dim(Fi)[1]){
          Xi_square <- sum(Fi[row,]**2)
          qlt[row,column] <- ((Fi[row,column]**2 * 100) / Xi_square) * sign(Fi[row,column])
        }
      }
      return(qlt)
    }

    contribution_abs_individuals <- function(Fi, D, eigenvalues, axisNumber){
      ctr <- matrix(0, nrow=dim(Fi)[1], ncol=axisNumber)
      for(column in 1:axisNumber){
        for(row in 1:dim(Fi)[1]){
          P_i <- D[row,row]
          ctr[row,column] <- (Fi[row,column]**2  * P_i * 100)/ eigenvalues[column]
        }
      }
      return(ctr)
    }
    # row contribution
    pca_inter <- append(pca_inter, list("qlt_Fi"=contribution_rel_individuals(Fi, axisNumber), "ctr_Fi"=contribution_abs_individuals(Fi, D, eigenvalues, axisNumber)))

    contribution_rel_variables <- function(Gi, axisNumber){
      qlt <- matrix(0, nrow=dim(Gi)[1], ncol=axisNumber)
      for(column in 1:axisNumber){
        for(row in 1:dim(Gi)[1]){
          Xj_square <- sum(Gi[row,]**2)
          qlt[row,column] <- ((Gi[row,column]**2 * 100) / Xj_square) * sign(Gi[row,column])
        }
      }
      return(qlt)
    }

    contribution_abs_variables <- function(Gi, Q, eigenvalues, axisNumber){
      ctr <- matrix(0, nrow=dim(Gi)[1], axisNumber)
      for(column in 1:axisNumber){
        for(row in 1:dim(Gi)[1]){
          S_i <- Q[row,row]
          ctr[row,column] <- (Gi[row,column]**2  * S_i * 100) / eigenvalues[column]
        }
      }
      return(ctr)
    }
    # columncontribution
    pca_inter <- append(pca_inter, list("qlt_Gi"=contribution_rel_variables(Gi, axisNumber), "ctr_Gi"=contribution_abs_variables(Gi, Q, eigenvalues, axisNumber)))

    return(pca_inter)
  }

  li_res <- append(li_res, list("PCA_interpretation"=interpretation_pca(li_res$PCA_component$Fi, li_res$PCA_component$Gi, Q, D, li_res$PCA_component$values, li_res$PCA_component$axisNumber)))
  return(li_res)
}

M <- matrix(c(1:3, c(7,3,5), c(5,11,30), c(20,15,1), c(20,40,19)), nrow=5, byrow=T)
print(pcafunction(M))
