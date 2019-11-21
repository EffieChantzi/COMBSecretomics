# COMBSecretomics

COMBSecretomics is a pragmatic theoretical and computational framework for secretome-based higher-order combination analysis. It offers the unique opportunity to study how multidrug treatments (i.e., more than 2 drugs) affect the protein release patterns of inactive/unprovoked and active/provoked cell populations. This molecular based evaluation approach of combination treatments opens for systematic studies associated with the secretome of complex diseases where there are still unmet diagnostic anf therapeutic needs. 

COMBSecretomics is compatible with any type of microtiter well format and measurement technology that gives values proportional to the corresponding protein concentrations of interest. There are two requirements for the performed experiments though:

1. Cell populations associated with healthy and disease cells must be studied in parallel on the same microtiter plate.
2. The experimental design must be based on a set of intra-plate (technical) replicate measurements. 


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
```
MATLAB R2018a or a later version is recommended
```
### Download 
```
$ git clone  https://github.com/EffieChantzi/COMBSecretomics.git
```

#### Directories
- Code: directory containing the source code 
- Case Study: directory containing the data of the case study

## Case Study
These instructions will help you test COMBSecretomics by using the data from our case study.

### Instructions

1. Open MATLAB

2. Navigate to the Code directory (see above)

3. Run the main.m script:
```
 main
```

4. Select interactively the raw data file 
```
3510.csv
``` 
from the Case Study directory (see above).

5. Enter cut-off threshold for blank filtering (%):
```
15
```

6. Enter cut-off threshold for coefficient of variation (%):
```
25
```

7. Enter number of resamplings:
```
500
```

8. Perform exhaustive subset search for K-means clustering?
```
2
```

9. Select mode for COMBSecretomics:
```
2
```

Running COMBSecretomics with the inputs from steps 5-9 above will reproduce the case study results for the Response Analysis.
All these inputs can be modified accordingly based on your preference. Feel free to explore!

#### Testing
All generated results will be saved in the subdirectory COMBSecretomics under the directory Case Study.

## Authors

* **Effie Chantzi**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
