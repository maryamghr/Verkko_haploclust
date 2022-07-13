# Haploclust
Snakemake pipeline for overlap graph-based haplotype phasing using Strand-seq and CCS reads.

Collaborators: Maryam Ghareghani, David Porubsky, and Tobias Marschall. 

## Setup

1. **Setup**

	Clone the git repository and checkout to the verkko branch.
	```
	git clone https://github.com/maryamghr/haploclust
    cd haploclust/pipeline
	git checkout verkko
	```

	Then install BubbleGun from `https://github.com/fawaz-dabbaghieh/bubble_gun`.

2. **Add your single-cell Strand-seq data**

	Create the subdirectory `../../ss_fastq/{sample_name}/`. Your Strand-seq FASTQ files of this sample go into this folder. The file names should be `{lib_name}_1.fastq.gz` and `{lib_name}_2.fastq.gz` for the first and second mates of reads per Strand-seq library.


3. **Add your CCS reads or the overlap graph**

	For CCS reads input: create the folder `../../ccs_reads` containing CCS reads with the name `{sample_name}.fastq.gz`.<br/>
        For overlap graph input: create the folder `../../{sample}/{assembly}` containing the overlap graph named `asm.bp.r_utg.gfa`.


5. **Adapt the config file**

6. **Create and activate conda environment**
    from `haploclust.yml` file.

7. **Run snakemake**

### NOTE

The Haploclust package is currently under development and contains unpublished work. Any usage for publishing is strictly prohibited without permission.
