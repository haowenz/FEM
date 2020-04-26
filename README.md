# FEM
FEM is a Fast and Efficient short read Mapper. Currently, FEM can return all mapping locations of NGS single-end short reads with up to 7 errors.

## Usage
### Indexing
```
FEM index <window_size> <step_size> <reference> <output>
```

### Mapping
```
Usage:  FEM map [options]

Options:
        -e       INT  error threshold
        -t       INT  number of threads
        -f       STR  seeding algorithm: "g" for group seeding and "v" for variable-length seeding
        -a       INT  # additional q-grams (only for test)

Input/output:
        --ref    STR  Input reference file
        --index  STR  Input index file
        --read1  STR  Input read1 file
        -o       STR  Output SAM file
```

## Parameters
Note that there is upper bound on the step size. Given a read of length _l_, window size _k_ and error threashold _e_, the max step size is _l/(e+2) − k + 1_. More details on this can be found in the paper.

To reduce mapping time, we recommend to use the smallest step size as long as the index can fit into the memory. In next version, FEM will choose the step size according to given memory adaptively. 

Since group seeding is usually sensitive enough and more efficient than variable-length seeding, we removed the implementation of variable-length seeding in the latest version. But you can find it in v0.1.

## Citing FEM
If you use FEM, please cite:

Zhang, H., Chan, Y., Fan, K. et al. Fast and efficient short read mapping based on a succinct hash index. BMC Bioinformatics 19, 92 (2018). https://doi.org/10.1186/s12859-018-2094-5
