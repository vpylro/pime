## How to install PIME package

To install PIME first install the devtools package. 
Then load the library(devtools) and run install_github using the following commands.
```{r}
install.packages("devtools")
library(devtools)
install_github("microEcology/pime")
```
PIME uses a Phyloseq object as input. A description of the phyloseq object and a tutorial on how to create this file in R using OTU tables in many different formats is detailed into the Phyloseq website https://joey711.github.io/phyloseq/ 

##Step-by-step example.
The first step in PIME is to define whether the microbial community presents a high relative abundance of taxa with low prevalence, which is defined as noise in PIME analysis. This is defined by Random Forests analysis. In this example we run PIME using the restroom dataset (https://doi.org/10.1007%2Fs10482-017-0976-6) against the metadata variable called Environment (a variable with two categories: men’s and women’s restroom). 

The prediction using random forests trees on full dataset. Results in Out of Bag erro rate
```{r}
library(pime)
data("restroom")
pime.oob.error(restroom, "Environment")
```
If the OOB error rate <=0.1, the dataset present large differences, and pime might not be necessary. 
If the OOB error rate is greater, next functions will find the best prevalence interval for the dataset.
This function takes two parameters: The phyloseq object (restroom) and the predictor variable (Environment).

## Split the dataset by predictor variable
Two parameters are required to run this function: The phyloseq object (restroom) and the predictor variable (Environment).
```{r}
per_variable_obj= pime.split.by.variable(restroom, "Environment")
per_variable_obj
```

## Calculate the highest possible prevalence intervals
This function calculates prevalence for different intervals by increments of 5. 
The input file is the output from the pime.split.by.variable (per_variable_obj)
```{r}
prevalences=pime.prevalence(per_variable_obj)
prevalences
```

## Calculate the best prevalence interval for the dataset 
This function will return a table with PERMANOVA  and random forests results for each prevalence interval. The number of taxa and the number of remaining sequences for each prevalence interval are also computed. 
The best prevalence interval value provides the clearest separation of communities while still including a majority of the taxa in the analysis. If true differences are present.
It will be represented by the first interval in which the OOB error rate is zero or close to zero.
The input file is the list of prevalences generated by the pime.prevalence (prevalences) and the predictor variable (Environment).
```{r}
pime.best.prevalence(prevalences, "Environment")
```
## Within this dataset the best prevalence interval was 60%
To obtain the phyloseq object at this cutoff use the following command.

```{r}
prevalence.60 = prevalences$`60`
prevalence.60
```
## Estimating prediction error
To estimate error in prediction, we will use pime.error.prediction() to randomly assign treatments to samples and run random forests classification on each prevalence interval. The function returns in a boxplot and a table with results of each classification error. For the purposes of this example we are running only 10 randomizations for saving time but we recommend at least 100 randomizations to obtain reliable results.
```{r}
#randomized=pime.error.prediction(restroom, "Environment", bootstrap = 10, parallel = TRUE, max.prev = 95)
#randomized$Plot
#randomized$`Table results'
```
It is also possible to estimate the variation of OOB error with each prevalence interval fintering. This is done by running the random forests classification for n times, determined by the user. This function will return a boxplot figure and a table for each classificaiton error.
```{r}
replicated.oob.error= pime.oob.replicate(prevalences, "Environment", bootstrap = 10, parallel = TRUE)
