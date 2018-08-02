# Ultrafast clustering of single-cell flow cytometry data using FlowGrid
    Authors: Xiaoxin Ye and Joshua W. K. Ho
    Contact: j.ho@victorchang.edu.au
    Copyright © 2018, Victor Chang Cardiac Research Institute
## Input data format
Our FlowGrid algorithm could be applied into many format data set but the sample code only accept csv format. In the csv file, the first row is feature name and each columns is seperated by "`,`". If you have true label file , you could use --l filename to input label file for testing the ARI of FlowGrid result

## Install
Before using the package, we need to install the dependent package sklearn and numpy.
``` python
sudo pip install -r requirements.txt
```
or
``` python
sudo pip install sklearn numpy
```
## Usage
A summary of the argument of sample code is included in the table below.
 

|argument | usage| Require |
| :----: | :----: | :----: |
|--f | the input file name| True|
|--n | number of bins | True |
|--eps | maximun distance between two bins| True |
|--t | threshold for high density bin| False(default:40) |
|--o | the output file name| False(default: out.csv) |
|--l | the true label file name| False|

## Sample
After installing the dependent packages, you could try to use the sample code to run FlowGrid on sample data.
``` bash
python sample_code.py --f sample_data.csv --n 4  --eps 1.1 --l sample_label.csv
```
The predicted label is saved at out.csv and the sample result is as follow.
``` 
The number of cells is: 23377
The number of dimensions is: 4
runing time: 0.027
ARI:0.9816
``` 
