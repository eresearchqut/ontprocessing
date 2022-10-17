# ONTProcessing workflow

Roberto Barrero, 17/10/2022

Craig Windell, 17/10/2022

## Usage:
Run the command
```
nextflow run eresearchqut/ontprocessing {optional arguments}...
```

## Optional arguments:
-  -resume                           Resume a failed run
-  --outdir                          Path to save the output file
                                    'results'
-  --samplesheet '[path/to/file]'    Path to the csv file that contains the list of
                                    samples to be analysed by this pipeline.
                          Default:  'index.csv'

  Contents of samplesheet csv:
```
sampleid,sample_files,reference
SAMPLE01,/user/folder/*.fastq.gz,/path/to/reference.fasta
---
*sample_files can refer to a folder with a number of*
*files that will be merged in the pipeline*
```

- --flye_read_error: adjust parameters for given read error rate (as fraction e.g. 0.03)

Default:  0.03

- --flye_ont_mode: Select from nano-raw, nano-corr, nano-hq

Default:  'nano-hq'

- --nanoq_code_start: Start codon position in the reference sequence

Default:  1

- --nanoq_read_length: Length cut off for read size

Default:  9000

- --nanoq_num_ref: Number of references used in the alignment

Default:  1

- --nanoq_qual_threshhold: Base quality score cut off

Default:  5

- --nanoq_jump:Increase this to make larger read intervals

Default:  10
