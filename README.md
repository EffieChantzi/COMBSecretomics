# COMBSecretomics

COMBSecretomics is a computational framework for secretome-related higher-order drug combination analysis. It offers the unique opportunity to study how multidrug treatments (i.e., more than 2 drugs) affect the protein release patterns of inactive/unprovoked and active/provoked cell populations. This molecular based evaluation approach of combination treatments opens for systematic studies associated with the secretome of complex diseases where there are still unmet diagnostic and therapeutic needs. 

COMBSecretomics is compatible with any type of microtiter well format and measurement technology (typically antibody-based multiplex assays or mass spectrometry) that gives values proportional to the corresponding protein concentrations of interest. There are two requirements for the performed experiments:

1. Cell populations associated with healthy and disease cells must be studied in parallel on the same microtiter plate to avoid batch (inter-plate) variability.
2. The experimental design must be based on a set of intra-plate (technical) replicate measurements in order to perform the associated non-parametric resampling-based statistical analyses. 


## Getting Started

This guide explains how to employ COMBSecretomics using the data of our example case study 
or any other study provided that the raw protein release data are on the required format 
(see sections "Example raw data file" and "User-defined inputs" in the Supplementary 
Information for more details).

### Prerequisites for running the source code of COMBSecretomics
```
MATLAB: version R2018a or later, Machine Learning Toolbox, Bioinformatics Toolbox.
OS: Windows, Linux, Mac OS X.
```

### If you do not have access to MATLAB
```
COMBSecretomics is also provided as a command line tool that can be deployed as a 
standalone executable on Windows machines.
```

### Download 
```
$ git clone  https://github.com/EffieChantzi/COMBSecretomics.git
```

### Provided Directories
- Windows_Executable: command line version of COMBSecretomics for Windows when there is not access to MATLAB
- source_code: source code of COMBSecretomics 
- case_study: test data of our example case study

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
