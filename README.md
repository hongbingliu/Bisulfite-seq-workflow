# Bisulfite-seq-workflow

### Processing Bisulfite-seq data from Bismark

#### Contain :

main script : BS-seq-processing.py

#### Environment : 

Linux, Python 3+, pandas, numpy

#### Useages :

A bisulfite-seq data workflow.

  -i : genes info (split by tab)
  
  -o : output name
  
  -f : gene annotation file (eg : gtf)
  
  -s : bismark file ( the output file from step `bismark_methylation_extractor` : "XXX_context_test_dataset.fastq_bismark.txt" )
  
  -d : divide genes into <int> (default : 1)
