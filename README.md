# scWECTA
scWECTA is a soft-weighted ensemble model for single-cell cell-type annotation using five different informative gene sets and integrating five classifiers.<br /> 

## Usage
The workflow of scWECTA mainly contains three steps.<br />  

* Data preprocessing<br />
* Model training<br/>
* Model predicting<br/> 

### Data preprocessing
Suppose the working directory is path, and user should put all three folders in our project into it. Besides, it needs to include the other two folders, data and data_processed, which are used to store the original scRNAseq data and the results of data processing respectively.<br />  

The data preprocessing is mainly done with preprocess.R, which taking scRNA-seq raw count matrix (named train.csv and test.csv respectively) and the label matrix of train data (named train_label.csv) as input. The command line under linux is as follows,
```
Rscript path/R/preprocess.R path/data path/data_processed [T or F] 1000
```
The first parameter in this command is the data path of all input data, and user need to make sure train.csv, train_label.csv and test.csv are all in this directory. The second parameter is the output path of processed data. The third parameter denotes whether you want to process train and test dataset respectly (F) or you want to process train and test dataset together (T) with common genes. The final parameter is the number of highly variablegenes user want to choose.<br />

After data preprocessing, results related to training data will be put into path/data_processed/train and result relevant to the test data will be stored into path/data_processed/test. If data_processed does not exist, the code will automatically create the directory.<br />  

### Model training and model predicting
Make sure that the gene marker list, all_cell_markers.txt (downloaded from CellMarker) is in the path/marker. Besides that, user need to create a new output folder in path, named output, to store our final cell type annotation result. The result is a csv table, named "anno.csv", which contains a unique column ('label') of cell type annotation result.<br />  

Model training and model predicting contain in the scWECTA.py, and user could run it in linux with python3. The command line is as follows,
```
python3 path/python/scWECTA.py --train path/data_processed/train/ --test path/data_processed/test/ --marker path/marker/ --output path/output/ -s Human -t Pancreas --thred 0.5 --sim spearman --cal sigmoid
``` 

### Packages used in python3 and R (4.1.0)
Before using scWECTA, you should make sure that all packages listed below are installed in python and R.<br />    

#### python
* numpy 
* pandas  
* scipy 
* scikit-learn

#### R
* Seurat  
* scran
* M3Drop  
* SingleCellExperiment  
* xfun
