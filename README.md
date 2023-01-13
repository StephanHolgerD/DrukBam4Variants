# variantfinder
### `variantfinder` is a comandline tool to search for variant supporting alignment in bam files. The results can be added as an INFO field in the vcf file or exported as a tsv. If needed or for debugging, variant supporting alignments can be exported to bam files. 



# Runtime

* Depending on the seq. depth analyzing 1000 variants ~ 5 seconds 
* writing a 2 bam files one with supporting and one with non-supporting aligments per variant ~ 2000 bam files --> 60 seconds total runtime. 
* writing bam files is at the moment not multiprocessed, depending on the writing speed time of the used hard drive more cpu could speed up the process


# Installing

## requirements

* pysam
* numpy
* tqdm

install via pip:
```
pip install variantfinder

```



# Usage


  ```
  usage: variantfinder vcf [-h] -b BAM -v VCF -o VCFOUT -f FASTA [--cpu CPU] [--bamsout] [--bamsdir BAMSDIR]

options:
  -h, --help            show this help message and exit

required arguments:
  -b BAM, --bam BAM     Pos. sorted and indexed bam file
  -v VCF, --vcf VCF     vcf file with variants of interest
  -o VCFOUT, --vcfout VCFOUT
                        vcf output
  -f FASTA, --fasta FASTA
                        fasta file of reference

optional arguments:
  --cpu CPU             number of cpu's to run in paralell
  --bamsout             create bams with variant supporting alignments and non supporting aligments
  --bamsdir BAMSDIR     folder for bams

```

  ```
 usage: variantfinder tsv [-h] -b BAM -v VCF -o TSVOUT -f FASTA [--cpu CPU] [--bamsout] [--bamsdir BAMSDIR]

options:
  -h, --help            show this help message and exit

required arguments:
  -b BAM, --bam BAM     Pos. sorted and indexed bam file
  -v VCF, --vcf VCF     vcf file with variants of interest
  -o TSVOUT, --tsvout TSVOUT
                        vcf output
  -f FASTA, --fasta FASTA
                        fasta file of reference

optional arguments:
  --cpu CPU             number of cpu's to run in paralell
  --bamsout             create bams with variant supporting alignments and non supporting aligments
  --bamsdir BAMSDIR     folder for bams


```


