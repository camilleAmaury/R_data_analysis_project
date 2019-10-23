library(ade4)

M <- t(read.csv(file="Projet M1 AD 1920.csv", header=FALSE, sep=","))
acp <- pcafunction(M)

# function which compare each value of a two matrix and tell if equal or not
matequal <- function(x, y, precision=3)
  is.matrix(x) && is.matrix(y) && dim(x) == dim(y) && all(round(cov_R,precision) == round(cov_Package,precision))

# centered matrix test
center_R <- apply(M,2,scale,center=TRUE,scale=FALSE)
center_Package <- acp$base_component$M_centered
# result should be TRUE
print(matequal(center_R, center_Package))

# Covariance test --> R covariance matrix is not biased (so 1/n-1)
cov_R <- cov(M)
cov_Package <- acp$base_component$S
# result should be TRUE
print(matequal(cov_R, cov_Package))

# centered reduced matrix test
scale_R <- apply(M,2,scale,center=TRUE,scale=TRUE)
scale_Package <- acp$base_component$M_scale
# result should be TRUE
print(matequal(scale_R, scale_Package))

# correlation test
cor_R <- cor(M)
cor_Package <- acp$base_component$R
# result should be TRUE
print(matequal(cor_R, cor_Package))
