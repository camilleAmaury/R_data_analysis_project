library(ade4)

M <- t(read.csv(file="dataset.csv", header=FALSE, sep=","))

# our function
acp_package <- pcafunction(M, axisNumber=2)
ax <- acp_package$PCA_component$axisNumber

# function which compare each value of a two matrix and tell if equal or not
matequal <- function(x, y, precision=0.001, samesign=TRUE){
   return(dim(x) == dim(y) && all(x + (if(samesign) -y else y) <= precision))
}


# --------------------------------------------------------- BASE COMPUTATION FOR PCA --------------------------------------------------------- #


# centered matrix test
print("centered matrix test : ")
center_R <- apply(M,2,scale,center=TRUE,scale=FALSE)
center_Package <- acp_package$base_component$M_centered
# matrix print
print(center_R[1:2,1:2])
print(center_Package[1:2,1:2])
# result should be TRUE
print(matequal(center_R, center_Package))


# Covariance matrix test
print("Covariance matrix test : ")
# R's covariance matrix is not biased (so 1/n-1)
cov_R <- cov(M)
# By default, our pca function gives the non-biaised matrix too
cov_Package <- acp_package$base_component$S
# matrix print
print(cov_R[1:2,1:2])
print(cov_Package[1:2,1:2])
# result should be TRUE
print(matequal(cov_R, cov_Package))


# Centered and reduced matrix test
print("Centered and reduced matrix test : ")
#  R's scale matrix implements a n-1 formula
scale_R <- apply(M,2,scale,center=TRUE,scale=TRUE)
# Our function doesn't apply the n-1 on the scale matrix
scale_Package <- acp_package$base_component$M_scale
# matrix print : as we ca see result are slighty different due to the n-1 formula
print(scale_R[1:2,1:2])
print(scale_Package[1:2,1:2])
# result should be FALSE
print(matequal(scale_R, scale_Package))


# Correlation matrix test
print("Correlation matrix test : ")
cor_R <- cor(M)
cor_Package <- acp_package$base_component$R
# matrix print
print(cor_R[1:2,1:2])
print(cor_Package[1:2,1:2])
# result should be TRUE
print(matequal(cor_R, cor_Package))


# --------------------------------------------------------- MAIN PCA COMPONENT --------------------------------------------------------- #


acp_R <- dudi.pca(M, scannf=FALSE, nf=2)


# X acp's triplet matrix test
print("X acp's triplet matrix test : ")
# the R's pca function is by default specified with scale=TRUE, center=TRUE --> so the X of (X,Q,D) matrix should be Xcr (scale matrix)
X_R <- acp_R$tab
# our pca's function has the same behavior
X_package <- acp_package$PCA_component$X
# matrix print
print(X_R[1:2,1:2])
print(X_package[1:2,1:2])
# result should be TRUE
print(matequal(X_R, X_package))


# Eigen values test
print("Eigen values test : ")
# dudi.pca package gives us all values > 0 (regarding to the rank)
eigval_R <- acp_R$eig
# our function gives only the selected number of axis/values regading to a specified number or the elbow/kaiser method
eigval_package <- acp_package$PCA_component$values
# vectors print
print(eigval_R)
print(eigval_package)
# result should be TRUE
print(all(round(eigval_R,3) == round(eigval_package,3)))


# Principal axes test (eigen vectors)
print("Principal axes test (eigen vectors) : ")
eigvec_R <- acp_R$c1
eigvec_package <- acp_package$PCA_component$vectors
# matrix print
print(eigvec_R[1:5,1:ax])
print(eigvec_package[1:5,1:ax])
# result should be TRUE
print(matequal(eigvec_R, eigvec_package[,1:ax], samesign=FALSE))


# Conservated axis number
print("Conservated axis number : ")
# Since dudi.pca doesn't implements automatic axis determination, we decided to compare manually with our function
ax_R <- acp_R$nf
ax_package <- acp_package$PCA_component$axisNumber
# Axis number print
print(paste0("R : ", ax_R, ", package : ", ax_package))
# result should be TRUE
print(ax_R == ax_package)


# Principal component (Fi)
print("Principal component (Fi) : ")
Fi_R <- as.matrix(acp_R$li)
Fi_package <- acp_package$PCA_component$Fi
# matrix print
print(Fi_R[1:5,1:2])
print(Fi_package[1:5,1:2])
# result should be TRUE
print(matequal(Fi_R, Fi_package[,1:ax], samesign=FALSE))


# Variable projection (Gi)
print("Variable projection (Gi) : ")
Gi_R <- as.matrix(acp_R$co)
Gi_package <- acp_package$PCA_component$Gi
# matrix print
print(Gi_R[1:5,1:2])
print(Gi_package[1:5,1:2])
# result should be TRUE
print(matequal(Gi_R, Gi_package[,1:ax], samesign=FALSE))


# Inertia of the incidence matrix
print("Inertia of the incidence matrix : ")
# dudi.pca doesn't give Inertia by default
inertia_R <- sum(acp_R$eig)
inertia_package <- acp_package$PCA_component$inertia
# Global Inertia print
print(paste0("R : ", inertia_R, ", package : ", inertia_package))
# result should be TRUE
print(round(inertia_R,3) ==  round(inertia_package,3))


R_inertia_row <- inertia.dudi(acp_R, row.inertia=TRUE)


# Information caught per axis
print("Information caught per axis : ")
# inertia.dudi strangely gives an inertia and cum column which are false (not percentage and not [0,1] numbers)
ic_R <- as.vector(as.matrix(R_inertia_row$tot.inertia)[,1])/10  #acp_R$eig/sum(acp_R$eig)*100
ic_package <- acp_package$PCA_component$axisIC
# vectors print
print(ic_R[1:5])
print(ic_package[1:5])
# result should be TRUE
print(all(round(ic_R,3) == round(ic_package,3)))


# Cumulative information caught per axis
print("Cumulative information caught per axis : ")
# dudi.pca doesn't give cumulative information caught by default
cumulative_ic_R <- as.vector(as.matrix(R_inertia_row$tot.inertia)[,3]) #cumsum(acp_R$eig/sum(acp_R$eig)*100)
cumulative_ic_package <- acp_package$PCA_component$cumulative_axisIC
# vectors print
print(cumulative_ic_R[1:5])
print(cumulative_ic_package[1:5])
# should be true
print(all(round(cumulative_ic_R,3) == round(cumulative_ic_package,3)))


# --------------------------------------------------------- PCA INTERPRETATION RESULT --------------------------------------------------------- #


# Relative contribution (qlt) for individuals
print("Relative contribution (qlt) for individuals : ")
qlt_indi_R <- as.matrix(R_inertia_row$row.rel)
qlt_indi_package <- acp_package$PCA_interpretation$qlt_Fi
# matrix print
print(qlt_indi_R[1:5,])
print(qlt_indi_package[1:5,])
# should be true
print(matequal(qlt_indi_R, qlt_indi_package[,1:ax],  samesign = FALSE))


# Absolute contribution (ctr) for individuals
print("Absolute contribution (ctr) for individuals : ")
ctr_indi_R <- as.matrix(R_inertia_row$row.abs)
ctr_indi_package <- acp_package$PCA_interpretation$ctr_Fi
# matrix print
print(ctr_indi_R[1:5,])
print(ctr_indi_package[1:5,])
# should be true
print(matequal(ctr_indi_R, ctr_indi_package[,1:ax]))


R_inertia_col <- inertia.dudi(acp_R, col.inertia=TRUE)


# Relative contribution (qlt) for variables
print("Relative contribution (qlt) for variables : ")
qlt_var_R <- as.matrix(R_inertia_col$col.rel)
qlt_var_package <- acp_package$PCA_interpretation$qlt_Gi
# matrix print
print(qlt_var_R[1:5,])
print(qlt_var_package[1:5,])
# should be true
print(matequal(qlt_var_R, qlt_var_package[,1:ax], samesign = FALSE))


# Absolute contribution (ctr) for variables
print("Absolute contribution (ctr) for variables : ")
ctr_var_R <- as.matrix(R_inertia_col$col.abs)
ctr_var_package <- acp_package$PCA_interpretation$ctr_Gi
# matrix print
print(ctr_var_R[1:5,])
print(ctr_var_package[1:5,])
# should be true
print(matequal(ctr_var_R, ctr_var_package[,1:ax]))

print(acp_package$base_component$colnames)

# test pca dudi
s.corcircle(Gi_R,xax=1,yax=2)
s.corcircle(Gi_package,xax=1,yax=2)
