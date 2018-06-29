# Ultrafast clustering of single-cell flow cytometry data using FlowGrid
    Authors: Xiaoxin Ye and Joshua W. K. Ho
    Contact: j.ho@victorchang.edu.au
    Copyright © 2018, Victor Chang Cardiac Research Institute
## Input data format
Our FlowGrid algorithm only csv format and each columns is seperated by `,`.

## Install
Before using the package, we need to install the dependent package sklearn and numpy.
``` python
pip install -r requirements.txt
```
## Usage
A summary of the argument of sample code is included in the table below.
 

|argument | usage|
| :----: | :----: |
|--f | the input file name|
|--n | number of bins |
|--t | threshold for high density bin|
|--d | maximun distance between two bins|
|--o | the output file name|

## Sample
After installing the dependent packages, you could try to use the sample code to run FlowGrid on sample data.
``` bash
python sample_code.py --f sample.csv --n 4 --t 40 --d 1.1 --o out.csv
```
The result is saved at out.csv.