# R_data_analysis_project

## Global information
<br/>

### Authors :
<br/>
- **Aniss Ilyes BENTEBIB** _(1st year of Master degree [computer science] related to Natural Language Processing)_.<br/>
- **Camille-Amaury JUGE** _(1st year of Master degree [computer science] related to Data Science)_.<br/>
- This package is under the license : _CC-BY-4.0_, free to use except for commercial.<br/><br/>

### Objectives :
<br/>
- see the file nammed _"Exercise_Lise_Bellanger_Univ_Nantes.pdf"_.<br/>
- This project is related to the 1st year master Degree _"Analyse de de données"_ held by **Lise BELLANGER-HUSI** in 2019.<br/>
- see the final report of this project nammed _"JUGE_BENTEBIB_ADD_Project.pdf"_.<br/><br/>

### Language, Tools used :
<br/>
- R<br/>
- Anaconda 3 Environment<br/>
- Trello and Git<br/>
- Jupyter Notebook <br/><br/>

## Documentation
<br/>
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
<br/>
### Principal component analysis 
<br/>
<br/>
The pca function have three main categories :
- The base objects (matrix and calculation to prepare our pca)<br/>
- The pca's axes, factors, projection of variables and individuals and more ...<br/>
- The pca's interpreted datas (inertia, cumulative percentages, correlation)<br/><br/>

In the future, the pca's object could be used in other functions to draw multiple plots ... (coming soon)
<br/>

```R
#By Default
nantes_pca <- pca(M, bias=FALSE, scale=TRUE)
```

<br/>
This will return a list of pca components including :
<br/>

```R
# Base components
nantes_pca$base #the main matrix
nantes_pca$G #the gravity center of your matrix
nantes_pca$M_centered #the centered matrix
nantes_pca$S #the covariance matrix (if bias=TRUE S main diagonal is the biaised variance 1/n)
# Depending on the type of PCA you aim to do, the following objects won't appear
# --> "scale=TRUE" will enable thos following objects to appear
nantes_pca$1_div_sd #the matrix which contains on its diagonal the (1/standard deviation) of the column 
nantes_pca$M_scale #the centered and reduced matrix
nantes_pca$R #the correlation matrix
```