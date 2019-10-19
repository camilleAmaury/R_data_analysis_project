# R_data_analysis_project

## Global information
<br/>

### Authors :
<br/>
- <b>Aniss Ilyes BENTEBIB</b> <i>(1st year of Master degree [computer science] related to Natural Language Processing)</i>.<br/>
- <b>Camille-Amaury JUGE</b> <i>(1st year of Master degree [computer science] related to Data Science)</i>.<br/>
- This package is under the license : <i>CC-BY-4.0</i>, free to use except for commercial.<br/><br/>

### Objectives :
<br/>
- see the file nammed <i>"Exercise_Lise_Bellanger_Univ_Nantes.pdf"</i>.<br/>
- This project is related to the 1st year master Degree <i>"Analyse de de données"</i> held by <b>Lise BELLANGER-HUSI</b> in 2019.<br/>
- see the final report of this project nammed <i>"JUGE_BENTEBIB_ADD_Project.pdf"</i>.<br/><br/>

### Language, Tools used :
<br/>
- R<br/>
- Anaconda 3 Environment<br/>
- Trello and Git<br/>
- Jupyter Notebook <br/><br/>

## Documentation
<br/>

### Installation and beginning
<br/>
<br/>
Download the package on ... (coming soon).
<br/>
<br/>
After the package is installed, just use it as a normal library in R.
<br/>
<br/>

```R
library(nantes_statistical_analysis)
```
<br/>

### Principal component analysis 
<br/>
<br/>
The pca function have three main categories :<br/>
- The base objects (matrix and calculation to prepare our pca)<br/>
- The pca's axes, factors, projection of variables and individuals and more ...<br/>
- The pca's interpreted datas (inertia, cumulative percentages, correlation)<br/><br/>

In the future, the pca's object could be used in other functions to draw multiple plots ... (coming soon)
<br/>

```R
#By Default
nantes_pca <- pca(M, center=TRUE, scale=TRUE, bias=FALSE, Q=diag(dim(M)[2]), D=(1/dim(M)[1])*diag(dim(M)[1]), axisMethod="kaiser", axisNumber=-1)
```

<br/>
<br/>

```R
#Examples

#Centered pca
nantes_pca <- pca(M, scale=FALSE)
#biased-Centered pca
nantes_pca <- pca(M, scale=FALSE, bias=TRUE)
#Centered-reduced pca, if scale is true, whatever center is specified it will be also true
# note that bias is useless in a scale pca
nantes_pca <- pca(M, scale=TRUE)
# specified axis number, won't take in consideration the determination algorithm for the number of axis
nantes_pca <- pca(M, scale=TRUE, axisNumber=2)
# auto determine the number of axis by the kaiser or elbow method
nantes_pca <- pca(M, scale=TRUE, axisMethod="kaiser")
nantes_pca <- pca(M, scale=TRUE, axisMethod="elbow")
# for a totally different pca, just set all to FALSE, and you can change the individuals metric matrix D, and variables metric matrix Q
```

<br/>
This will return a list of pca components including :
<br/>

```R
# Base components
nantes_pca$base_component$base #the main matrix
nantes_pca$base_component$G #the gravity center of your matrix
nantes_pca$base_component$M_centered #the centered matrix
nantes_pca$base_component$S #the covariance matrix (if bias=TRUE S main diagonal is the biaised variance 1/n)
# Depending on the type of PCA you aim to do, the following objects won't appear
# --> "scale=TRUE" will enable thos following objects to appear
nantes_pca$base_component$D1_div_sd #the matrix which contains on its diagonal the (1/standard deviation) of the column 
nantes_pca$base_component$M_scale #the centered and reduced matrix
nantes_pca$base_component$R #the correlation matrix
```

<br/>
<br/>

```R
# PCA components
nantes_pca$PCA_component$Si #the incidence matrix
nantes_pca$PCA_component$values #eigen values of the incidence matrix
nantes_pca$PCA_component$vectors #eigen vectors of the incidence matrix
nantes_pca$PCA_component$Fi #Individuals matrix projection
nantes_pca$PCA_component$Gi #Variables matrix projection
```

<br/>

### Kmeans 

<br/>
<br/>
The categorical algorithm for classification K-means  :<br/>
- Returns the cluster as a unique vector with every individuals as an indice, and his value his cluster number<br/>
- Returns the clusters' centers point<br/>
- Implements multiple methods (distance, initialization) <br/><br/>

In the future, the K-means' object could be used in other functions to draw multiple plots ... (coming soon)
<br/>

```R
#By Default
nantes_kmeans <- kmeansfunction(M, k=-1, d="euclidian", initialization="random", precision=0.001, logs=FALSE)
```

<br/>
This will return a list compunded of clusters' centers point and cluster's membership of individuals :
<br/>

```R
# returns
nantes_pca$centers #the clusters' center points as a matrix
nantes_pca$clusters #cluster's membership of individuals as a vector
```

<br/>
<br/>

```R
#Examples

# initialization algorithm :
# random will take "k" points randomly in all individuals
nantes_kmeans <- kmeansfunction(M, initialization="random")
# we advize you to choose this one rather than the last one, see why here : https://en.wikipedia.org/wiki/K-means%2B%2B
# choose the first point randomly and then compute a probability for every other point regarding to the euclidian distance
nantes_kmeans <- kmeansfunction(M, initialization="kmeans++")

# distance choosen :
# basic euclidian distance between center and a point
nantes_kmeans <- kmeansfunction(M, d="euclidian")
# Distance of mahalanobis, which is sometimes more efficient because it doesn't take an account of outliers
nantes_kmeans <- kmeansfunction(M, d="mahalanobis")

# precision parameters
# precision limit between two iterations , note that if they are more than 1000 iterations it will stops even precision is not respected
nantes_kmeans <- kmeansfunction(M, precision=0.001)

# logs parameters
# if you want to print all steps (clusters and center points) else FALSE
nantes_kmeans <- kmeansfunction(M, logs=TRUE)

# k parameters
# the cluster number, note that if you let it to -1, it will automatically do 1:k clusters and choose the best one
nantes_kmeans <- kmeansfunction(M, k=-1) # k will be choosen within the algorithm
nantes_kmeans <- kmeansfunction(M, k=2) # k is specified, this will only create two clusters
```

<br/>

### Eigen reproduced function

<br/>
<br/>
We decided to work additionally on implementing our own Eigen function  :<br/>
- returns the same parameters than R's eigen function<br/>
- Use the QR algorithm, with the method of HouseHolders, which is a numerically stable algorithm rather than Givens ones.<br/>
- Is really efficient on small sized and symmetric matrix, else vectors suffers of round problems (we didn't find the solution to improve it) <br/><br/>


```R
#By Default
nantes_kmeans <- eigenfunction(M)
```