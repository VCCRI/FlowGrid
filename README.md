# Ultrafast clustering of single-cell flow cytometry data using FlowGrid
    Authors: Xiaoxin Ye and Joshua W. K. Ho
    Contact: j.ho@victorchang.edu.au
    Copyright © 2018, Victor Chang Cardiac Research Institute
## Input data format
Our FlowGrid algorithm only csv format and each columns is seperated by `,`.
## Usages
Before using the package, we need to install the dependent package sklearn and numpy.
``` python
pip install -r requirements.txt
```
After installing the dependent packages, you could try to use the sample code to run FlowGrid on sample data.
``` bash
python sample_code.py --f sample.csv --n 4 --t 40 --d 1.1
```