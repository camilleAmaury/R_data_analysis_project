library(ade4)

M <- t(read.csv(file="Projet M1 AD 1920.csv", header=FALSE, sep=","))
acp_package <- pcafunction(M, method="kaiser")

# function which compare each value of a two matrix and tell if equal or not
matequal <- function(x, y, precision=3)
  dim(x) == dim(y) && all(round(x,precision) == round(y,precision))

# centered matrix test
center_R <- apply(M,2,scale,center=TRUE,scale=FALSE)
center_Package <- acp_package$base_component$M_centered
# result should be TRUE
print(matequal(center_R, center_Package))

# Covariance test --> R covariance matrix is not biased (so 1/n-1)
cov_R <- cov(M)
cov_Package <- acp_package$base_component$S
# result should be TRUE
print(matequal(cov_R, cov_Package))

# centered reduced matrix test
scale_R <- apply(M,2,scale,center=TRUE,scale=TRUE)
scale_Package <- acp_package$base_component$M_scale
# result should be TRUE
print(matequal(scale_R, scale_Package))

# correlation test
cor_R <- cor(M)
cor_Package <- acp_package$base_component$R
# result should be TRUE
print(matequal(cor_R, cor_Package))

# acp by default scale=TRUE, center=TRUE --> so the X of (X,Q,D) matrix should be Xcr (scale matrix)
# our pca has the same behavior
acp_R <- dudi.pca(M, scannf=FALSE)
X_R <- acp_R$tab
X_package <- acp_package$PCA_component$X
print(matequal(X_R, X_package))

# eigen values
eigval_R <- acp_R$eig
eigval_package <- acp_package$PCA_component$values
print(all(round(eigval_R,3) == round(eigval_package,3)))

# Principal component
Fi_R <- acp_R$li
Fi_package <- acp_package$PCA_component$Fi
print(matequal(Fi_R, Fi_package))

# axis number
ax_R <- acp_R$nf
ax_package <- acp_package$PCA_component$axisNumber
print(paste0("R : ", ax_R, ", package : ", ax_package))

# inertia
inertia_R <- acp_R$eig/sum(acp_R$eig)*100
inertia_package <- acp_package$PCA_component$inertia
print(all(round(inertia_R,3) == round(inertia_package,3)))

#for(i in 1:dim(X_R)[1]){for(j in 1:dim(X_R)[2]){if(round(X_R[i,j],3) != round(X_package[i,j],3)){print(paste0(i, ",", j))}}}

#inertia.dudi(acp_R,row.inertia=TRUE)
#inertia.dudi(acp_R,col.inertia=TRUE)
#names(inertierow) ; inertierow$TOT
