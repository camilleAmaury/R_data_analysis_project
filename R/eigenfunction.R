eigen_function_hh <- function(M){
  norm <- function(v){
    return(sqrt(sum(v*2)))
  }
  QRDecompHH <- function(A){
    T <- A
    I <- diag(dim(A)[1])

    vi <- t(A[,1]) - sign(A[1,1])norm(A[,1])I[,1]


    hi <- I-(t(vi)%%vi)matrix(rep(2/vi%%t(vi), dim(A)[1]dim(A)[1]), nrow=dim(A)[1], ncol=dim(A)[1])


    his <- list(hi)
    Ak <- hi%%A
    for(i in 2:dim(Ak)[1]){
      if(dim(A)[1] > 2){
        A <- Ak[2:dim(Ak)[1], 2:dim(Ak)[2]]
        #print(A)
        #print(i)
        I <- diag(dim(A)[1])

        vi <- t(A[,1]) - sign(A[1,1])norm(A[,1])I[,1]
        #print(vi)

        hi <- I-(t(vi)%%vi)matrix(rep(2/vi%%t(vi), dim(A)[1]dim(A)[1]), nrow=dim(A)[1], ncol=dim(A)[1])
        #print(hi)

        his[[i]] <- hi
        Ak <- hi%%A
        #print(Ak)
      }
    }

    Q <- diag(dim(T)[1])
    R <- diag(dim(T)[1])
    for(i in 1:length(his)){
      I = diag(dim(T)[1])

      c <- his[[i]]
      I[(dim(I)[1]-dim(c)[1]+1):dim(I)[1], (dim(I)[1]-dim(c)[1]+1):dim(I)[1]] <- c
      his[[i]] <- I
      Q <- Q%%I

    }


    for(i in 1:length(his)){
      R <- R%%his[[dim(T)[1]-i]]
    }
    return(list(R%%T, Q))
  }

  eigen_v <- function(A, accuracy){
    Ao <- A
    Uo <- diag(dim(A)[1])
    for(i in 1:accuracy){
      Akk <- QRDecompHH(Ao)
      Ao <- Akk[[1]]%%Akk[[2]]
      Uo <- Uo%%Akk[[2]]
    }
    return(list(Uo, Ao))
  }

  dimM = dim(M)[1]
  eigen_v_ret <- eigen_v(M, 100)


  return(list("eigen_values"=diag(eigen_v_ret[[2]]), "eigen_vectors"=eigen_v_ret[[1]]))
}
