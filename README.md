# FEM
FEM is a Fast and Efficient short read Mapper. Currently, FEM can return all mapping locations of NGS single-end short reads with up to 7 errors.

## News
### 03/13/2018
We pushed the source code of FEM purely in C to replace the old hybrid one. 

### 03/10/2018
The source code of FEM is avaliable. We removed the dependence on Intel TBB library so that FEM can be easily used and tested. A parallel sorting without dependency on third party libraries will be added to FEM soon.

## Usage
### Indexing
```
FEM index <window size> <step size> <reference.fa> 
```

### Mapping
```
Usage:   FEM align [options] 

Options:
         -e        INT    error threshold 
         -t        INT    number of threads 
         -f        STR    seeding algorithm: "g" for group seeding and "vl" for variable-length seeding 
         -a               use one additional q-gram for filter

Input/output: 
         --ref     STR    Input reference file
         --read    STR    Input read file
         -o        STR    Output SAM file 
```


## Parameters
To reduce mapping time, we recommend to use the smallest step size as long as the index can fit into the memory. In next version, FEM will choose the step size according to given memory adaptively. 

## Citing FEM
If you use FEM, please cite [Fast and efficient short read mapping based on a succinct hash index](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2094-5).

## Contacts
Haowen Zhang <br />
Email: hwzhang@gatech.edu

We welcome any bug report and suggestion. Please start an issue in the repo.
