# scWECTA
scWECTA is a soft-weighted ensemble model for single-cell cell-type annotation using five different informative gene sets and integrating 
five classifiers.

## Usage
The workflow of scWECTA mainly contain three steps.<br /> 

(1) Data preprocessing<br />
(2) Model training<br/>
(3) Model predicting<br/>


### Data preprocessing
Suppose your working directory is named as path, and it needs to include two folders, data and output, which are used to store the original scRNAseq data and the results of data processing respectively.<br />
  
The data preprocessing is mainly done with preprocess.R, which taking scRNA-seq raw count matrix (train and test data) and the label matrix of train data as input. The command line under linux is as follows,
```
Rscript preprocess.R "path/data" "path/output" [T or F]
```
The first parameter in this command is the data path of all input data, and user need to make sure train.csv, train_label.csv and test.csv are all in this directory. The second parameter is the output path of processed data. The third parameter denotes whether you want to process train and test dataset respectly (T) or you want to process train and test dataset together (F) with common genes.<br />

After data preprocessing, results related to training data will be put into path/output/train and result related to the test data will be put into path/output/test.

### Model training and model predicting
Remember to put the marker list, all_cell_markers.txt in the data folder of our github project, into path/marker. Besides that, create a new output folder, named output_final, to store our final cell type annotation result. 

These two steps contain in the scWECTA.py, and user could run it in linux with python 3.7. The command line is as follows,
```
python scWECTA.py -train "path/train" -test "path/test/" -marker "path/marker/" -o "path/output_final/" -s "Human" -t "Pancreas" -thred 0.5
```
