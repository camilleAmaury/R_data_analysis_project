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
        li[["indi_pos"]] <- c(li[["indi_neg"]], j)
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
