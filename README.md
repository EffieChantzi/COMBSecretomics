# COMBSecretomics

COMBSecretomics is a computational framework for secretome-related higher-order drug combination analysis. It offers the unique opportunity to study how multidrug treatments (i.e., more than 2 drugs) affect the protein release patterns of inactive/unprovoked and active/provoked cell populations. This molecular based evaluation approach of combination treatments opens for systematic studies associated with the secretome of complex diseases where there are still unmet diagnostic and therapeutic needs. 

COMBSecretomics is compatible with any type of microtiter well format and measurement technology (typically antibody-based multiplex assays or mass spectrometry) that gives values proportional to the corresponding protein concentrations of interest. There are two requirements for the performed experiments:

1. Healthy and disease associated cells must be studied in parallel on the same microtiter plate to avoid batch (inter-plate) variability.
2. The experimental design must be based on a set of intra-plate (technical) replicate measurements in order to perform the associated non-parametric resampling-based statistical analyses. 


## Getting Started

This guide explains how to employ COMBSecretomics using the data of our example case study 
or any other study provided that the raw protein release data are on the required format 
(see sections "Example raw data file" and "User-defined inputs" in the Supplementary 
Information of our article for more details).

### Prerequisites for running the source code of COMBSecretomics
```
MATLAB: version R2018a or later, Machine Learning Toolbox, Bioinformatics Toolbox.
OS: Windows, Linux, Mac OS X.
```

### If you do not have access to MATLAB
```
COMBSecretomics is also provided as a command line tool that can be deployed as a 
standalone executable on Windows machines. In this case, download and installl the
Windows version of MATLAB Runtime for R2019b from the following link on the MathWorks website:  
```
https://www.mathworks.com/products/compiler/mcr/index.html


### Download 
You can download the COMBsecretomics project on your local machine by executing (terminal/Git bash):
```
$ git clone  https://github.com/EffieChantzi/COMBSecretomics.git
```
You can also use the 

### Provided Directories
- Windows_Executable: command line version of COMBSecretomics for Windows when there is not access to MATLAB
- case_study: test data of our example case study
- source_code: source code of COMBSecretomics 

### Getting started with the source code:

- Open MATLAB, navigate to 'source_code' and run the file 'main.m' (or simply type 'main' in the MATLAB 
  command window).
  
### Getting started with the Windows executable:

- Navigate to 'Windows_Executable', right click on the the executable file COMBSecretomics.exe 
  and choose the option 'Run as Administrator'. Click 'Yes' on the User Account Control.
  
### Sequence of user-defined inputs

```
1. The user is prompted to select interactively the CSV specification file, which is 
essentially the data file with the raw protein release measurements. It should be in 
a particular format in order to be valid (see section "Example raw data file" in the
Supplementary Information). It can be selected from the directory of your choice. For 
the current case study, select the file '3510.csv' from the directory 'case_study',
which is the file with the raw protein release measurements for our example case study.
'3510' denotes the barcode of the microtiter plate used for the corresponding in vitro 
experiment. You can use your own specification file if you want to analyze a new 
experiment provided that the filename is in the form <barcode.csv> and the contents of
the csv file agree with the instructions under the section "Example raw data file" in 
the Supplementary Information.
```

```
2. The user is asked to select the cut-off threshold (in percent) for the blank filtering
(see section "Blank filtering" in the Supplementary Information for more details). In our 
case the value 15 was used. Other values can be set accordingly. A low value ensures that
protein release measurements with low levels of noise are kept for further analysis.
```
```
3. The user is asked to select the cut-off threshold (in percent) for the coefficient of 
variation of the protein release measurements (see section "Coefficient of variation" in
the Supplementary Information for more details). In our case study the value 25 was used.
Other values can be set accordingly. A low value ensures that proteins with low levels of
technical variability are kept for further analysis.
```

```
4. The user is asked to enter the number of resampling based validation datasets to be 
created (see section "Resampling statistics" in the main article for more details). In 
our study the value 500 was used. Other values can be set accordingly. The higher this 
number, the more datasets are going to be created, which is advisable especially if 
several intra-plate replicate measurements are employed. 
```

```
5. The user is asked to select if exhaustive subset search should be performed or not for 
the visualization of the results from the top-down hierarchical K-Means clustering. Enter
1 for yes and 2 for no. In our case study, option 1 was used for the main article (Fig. 5) 
and option 2 for the supplementary information (Supplementary Fig. S18). It is highly 
recommended to use option 1 for big exhaustive combination panels (i.e., when more than 4 
single drugs are used to desing an exhaustive experiment).
```

```
6. The user is asked to select analysis mode regarding the type of cells to be analyzed. 
Enter 1 for the analysis of protein release measurements for unstimulated cells and 2 
for stimulated cells. The user is advised to run both modes sequentially (i.e., one after
the other).
```

### Results

The generated results can be found inside the folder 'case_study' under the automatically 
created directory with the name 'Results'. All results from the quality control and intra-plate
averaging are saved as colored figures in EPS format in addition to a CSV file, which contains a
table with all raw protein release measurements after quality control and before intra-plate 
averaging. Under 'Results', two subdirectories are created with the name 'Stimulated' and 
'Unstimulated' containing all results (in EPS format) from the analyses using the protein release
measurements for stimulated and unstimulated cells, respectively. 


## Authors

* **Effie Chantzi** (efthymia.chantzi@medsci.uu.se)

## License

This project is licensed under the GNU GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details
