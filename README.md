# Ultrafast clustering of single-cell flow cytometry data using FlowGrid
    Authors: Xiaoxin Ye and Joshua W. K. Ho
    Contact: j.ho@victorchang.edu.au
    Copyright © 2018, Victor Chang Cardiac Research Institute
## Input data format
Our FlowGrid algorithm could be applied into many format data set but the sample code only accept csv format. In the csv file, the first row is feature name and each columns is seperated by "`,`". If you have true label file , you could use --l filename to input label file for testing the ARI of FlowGrid result.

## Install
Before using the package, we need to install the dependent package sklearn and numpy.
``` Bash
pip install -r requirements.txt --user
```
or
``` Bash
pip install sklearn numpy scipy --user
```
## Usage
A summary of the argument of sample code is included in the table below.
 

|Argument | Usage| Required? |
| :----: | :----: | :----: |
|--f | the input file name| required |
|--n | number of bins | required |
|--eps | maximun distance between two bins| required |
|--t | threshold for high density bin| optional (default:40) |
|--o | the output file name| optional (default: out.csv) |
|--l | the true label file name| optional |

## Sample
After installing all the dependent packages, you could try to use the sample code to run FlowGrid on the sample data.
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
# FlowGrid for Scanpy usage
    Authors: Xiunan Fang and Joshua W. K. Ho
    Contact: xiunan@hkuh.k
    Copyright © 2019, Joshua W. K. Ho's Lab
## NOTE!
FlowGrid implementation in Scanpy can be found in github PAGE: https://github.com/holab-hku/FlowGrid

