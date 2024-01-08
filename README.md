### 1. Data structures and study design
Innthis project the main objection was to integrate two data sets to have a progression view on different stages in non small cell lung cancer (NSCLC).
  

### **3.Batch effect removal** 

Since we integrated two data sets, one of the most concerning consideration should be the batch effect. We might see non-biological

#### 1. LOESS

locally weighed scatter plot smothing

#### 2. Quintile Normalization

We can simply use `Matrix::RowMean()` to add a new column to our count matrix, which represent the avarage expression of each Feature within the samples. After having the mean of each feature, assing each gene count to the mean we have.
