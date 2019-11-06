# --------------------------------------------------------- Functions which can interpret and draw results --------------------------------------------------------- #
explainPCA <- function(pca, individualRate=-1, variableRate=-1){

  library(ade4)

  if(individualRate==-1){
    individualRate=(1/dim(pca$base_component$base)[1]*100)
  }
  if(variableRate==-1){
    variableRate=(1/dim(pca$base_component$base)[2]*100)
  }
  # create a tableau for each axis choosen, and determines the most contributing variables/individuals
  for(i in 1:pca$PCA_component$axisNumber){
    li <- list("indi_neg"=c(), "indi_pos"=c(), "var_neg"=c(), "var_pos"=c())
    indi <- which(pca$PCA_interpretation$ctr_Fi[,i] > individualRate)
    var <- which(pca$PCA_interpretation$ctr_Gi[,i] > variableRate)
    chaine_indi_min = ""
    chaine_indi_max = ""
    for(j in indi){
      if(pca$PCA_interpretation$qlt_Fi[j,i] < 0){
        li[["indi_neg"]] <- c(li[["indi_neg"]], j)
        chaine_indi_min = paste0(chaine_indi_min, pca$base_component$rownames[j], sep=",")
      }else{
        li[["indi_pos"]] <- c(li[["indi_pos"]], j)
        chaine_indi_max = paste0(chaine_indi_max, pca$base_component$rownames[j], sep=",")
      }
    }
    chaine_var_min = ""
    chaine_var_max = ""
    for(j in var){
      if(pca$PCA_interpretation$qlt_Gi[j,i] < 0){
        li[["var_neg"]] <- c(li[["var_neg"]], j)
        chaine_var_min = paste0(chaine_var_min, pca$base_component$colnames[j], sep=",")
      }else{
        li[["var_pos"]] <- c(li[["var_pos"]], j)
        chaine_var_max =  paste0(chaine_var_max, pca$base_component$colnames[j], sep=",")
      }
    }
    resultats <- data.frame("negatively"=c(chaine_indi_min, chaine_var_min),"positively"=c(chaine_indi_max, chaine_var_max),
                            row.names=c("individuals","variables"))

    string = ""
    for(col in 1:length(names(resultats))){
      for(row in 1:length(row.names(resultats))){
        string = cat(string, "\n", "\n", "The ", row.names(resultats)[row], " which have ", names(resultats)[col], " contributed on axis ", i, " are : ")
        string = cat(string,"\n", as.matrix(resultats)[row,col],"\n", "-------------")
      }
    }
    print(string)
    #plot correlation circle
    s.corcircle(pca$PCA_component$Gi[li[["var_neg"]],], xax = 1, yax = 2, label = pca$base_component$colnames[li[["var_neg"]]], clabel = 1, grid = TRUE, sub = paste0("Axis ", i," : variable negative contribution"), csub = 1, possub = "bottomleft", cgrid = 0, fullcircle = TRUE, box = FALSE, add.plot = FALSE)
    #plot correlation circle
    s.corcircle(pca$PCA_component$Gi[li[["var_pos"]],], xax = 1, yax = 2, label = pca$base_component$colnames[li[["var_pos"]]], clabel = 1, grid = TRUE, sub = paste0("Axis ", i," : variable positive contribution"), csub = 1, possub = "bottomleft", cgrid = 0, fullcircle = TRUE, box = FALSE, add.plot = FALSE)
    #plot cloud of point
    plot(pca$PCA_component$Fi[li[["indi_neg"]],], type="p",  sub = paste0("Axis ", i," : individual negative contribution"))
    text(pca$PCA_component$Fi[li[["indi_neg"]],], label = pca$base_component$colnames[li[["indi_neg"]]], cex= 0.7, pos=3)
    #plot cloud of point
    plot(pca$PCA_component$Fi[li[["indi_pos"]],], type="p",   sub = paste0("Axis ", i," : individual positive contribution"))
    text(pca$PCA_component$Fi[li[["indi_pos"]],],label = pca$base_component$colnames[li[["indi_pos"]]], cex= 0.7, pos=3)
  }
}


explainDataset <- function(M, colnames, bidimensionnal=FALSE){

  library(e1071)

  # analysis monodimensionnal
  if(!bidimensionnal){
    # paramètre de position
    meanVar <- apply(M,2,mean)
    modeVar <- apply(M,2,mode)
    medianVar <- apply(M,2,median)
    dataframe_position <- data.frame("mean" = meanVar, "median"=medianVar, "mode"=modeVar)

    # paramètre de dispersion
    quantileVar <- apply(M,2,quantile)
    sdVar <- apply(M,2,sd)
    etendueVar <- apply(quantileVar,2,function(x){return(x[5] + if(sign(x[1]) == sign(x[5])) - x[1] else x[1])})
    interquantileVar <- apply(quantileVar,2,function(x){return(x[4] + if(sign(x[2]) == sign(x[4])) - x[2] else x[2])})
    dataframe_dispersion <- data.frame("standart deviation" = sdVar, "scope"=etendueVar, "interquantile deviation"=interquantileVar)

    # paramètre de forme
    kurtosisVar <- apply(M,2,kurtosis)
    skewnessVar <- apply(M,2,skewness)
    dataframe_forme <- data.frame("kurtosis" = kurtosisVar, "skewness"=skewnessVar)

    #
    for(i in 1:ncol(M)){
      hist(M[,i])
      boxplot(M[,i])
    }
    # position
    print("# ------------------------ # Comparison of position # ------------------------ # ")

    for(i in 1:length(meanVar)){
      mean_med_comparison <- if(meanVar[i] == medianVar[i]) "Centred Variable" else (if(meanVar[i] * 1.25 >= medianVar[i] && meanVar[i] * 0.75 <= medianVar[i]) "Not centred but no outliers" else "Outliers detected")
      print(cat("v", i, " : \n", "   mean : ", meanVar[i], "\n   median : ", medianVar[i], "\n   mode : ", modeVar[i], "\n   Comparison between mean and median : ", mean_med_comparison))
      print("------------------------")
    }

    # dispersion
    print("# ------------------------ # Comparison of dispersion # ------------------------ # ")
    for(i in 1:length(sdVar)){
      print(cat("v", i, " : \n", "   standart deviation : ", sdVar[i], "\n   scope (outliers can affect) : ", etendueVar[i], "\n   interquantile deviation : ", interquantileVar[i]))
      print("------------------------")
    }

    # forme
    print("# ------------------------ # Comparison of form # ------------------------ # ")
    for(i in 1:length(sdVar)){
      kurtosis_interpretation <- if(kurtosisVar[i] == 0) "flattened gauss distrubution" else (if(kurtosisVar[i] > 0) "more concentrated than Gauss" else "more flattened gauss distribution")
      skewness_interpretation <- if(skewnessVar[i] == 0) "Symmetric distribution" else (if(skewnessVar[i] > 0) "Converging to the right Distribution" else "Converging to the left Distribution")
      print(cat("v", i, " : \n", "   kurtosis rate : ", kurtosisVar[i], "\n   interpretation : ", kurtosis_interpretation, "\n   skewness rate : ", skewnessVar[i], "\n   interpretation : ", skewness_interpretation))
      print("------------------------")
    }

  }else{
    # Test du chi2 pour chaque paire de variable
    print("# ------------------------ # X2 # ------------------------ # ")
    for(i in 1:(ncol(M)-1)){
      newM <- abs(M[,i:i+1])
      khi <- chisq.test(newM)
      print(khi)
      print(paste0("Variables[", i, ",", i+1,"], Considérons un seuil à 0.05, on peut alors ", if(khi[["p.value"]] > 0.05) "rejeter l'hypothèse nulle d'indépendance" else "rejeter l'hypothèse d'indépendance"))
      print("------------------------")
    }

    meanVar <- apply(M,2,mean)
    modeVar <- apply(M,2,mode)
    medianVar <- apply(M,2,median)
    dataframe_position <- data.frame("mean" = meanVar, "median"=medianVar, "mode"=modeVar)

    for(i in 1:(ncol(M)-1)){
      linear_reg <- lm(M[,i] ~ M[,i+1])
      print(summary(linear_reg))
      plot(M[,i] ~ M[,i+1])
      abline(linear_reg, col = "red")
    }
  }
}

explainKmeans <- function(M, kmeans){
  plot(x = M[,1], y = M[,2], col=sample(c("antiquewhite2", "antiquewhite4", "aquamarine1", "aquamarine4", "azure4",
                                          "blue2", "brown1", "brown4", "chartreuse", "chartreuse4", "chocolate1",
                                          "cyan3", "darkgoldenrod1", "darkmagenta", "darkolivegreen1", "deeppink1",
                                          "gray8"))[kmeans$clusters])
  text(x = M[,1], y = M[,2], label = c(1:length(M[,1])), cex= 0.7, pos=3)
}
