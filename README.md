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
- case_study: test data of our example case study
- source_code: source code of COMBSecretomics 

### Getting started with the source code:

- Navigate to 'source_code', open the file 'main.m' with MATLAB and run it (or simply type 'main' in the MATLAB 
  command window).
  
### Getting started with the Windows executable:

- Navigate to 'Windows_Executable', right click on the the executable file COMBSecretomics.exe 
  and choose the option 'Run as Administrator'. Click 'Yes' on the User Account Control.
  
### Sequence of user-defined inputs

1. The user is prompted to select interactively the specification file, which is 
essentially the data file with the raw protein release measurements. It should be in 
a particular format in order to be valid (see section "Example raw data file" in the
Supplementary Information). It does not have to be saved in a particular directory for
COMBSecretomics to run.







Running COMBSecretomics with the inputs from steps 5-9 above will reproduce the case study results for the Response Analysis.
All these inputs can be modified accordingly based on your preference. Feel free to explore!

#### Testing
All generated results will be saved in the subdirectory COMBSecretomics under the directory Case Study.

## Authors

* **Effie Chantzi**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
