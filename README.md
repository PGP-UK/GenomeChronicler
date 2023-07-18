```ascii
 #####                                         #####
#     # ###### #    #  ####  #    # ######    #     # #    # #####   ####  #    # #  ####  #      ###### #####
#       #      ##   # #    # ##  ## #         #       #    # #    # #    # ##   # # #    # #      #      #    #
#  #### #####  # #  # #    # # ## # #####     #       ###### #    # #    # # #  # # #      #      #####  #    #
#     # #      #  # # #    # #    # #         #       #    # #####  #    # #  # # # #      #      #      #####
#     # #      #   ## #    # #    # #         #     # #    # #   #  #    # #   ## # #    # #      #      #   #
 #####  ###### #    #  ####  #    # ######     #####  #    # #    #  ####  #    # #  ####  ###### ###### #    #
```

# GenomeChronicler

This is the repository for GenomeChronicler, the Personal Genome Project United Kingdom (PGP-UK) genomic report generation scripts.

## Dependencies

* [Docker](https://docs.docker.com/) or [Singularity](https://sylabs.io/guides/3.1/user-guide)

## Input Files

* BAM file OR VEP file
* VEP generated summary HTML file (optional)

## Running GenomeChronicler

### Obtaining test data

Getting some test data (NA12878 from ENA, pre-mapped to GRCh38, and the respective reference) and store this in your current folder:

```bash
# CRAM file for NA12878
curl -L -o na12878wxs.cram ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram

# gVCF for NA12878
# TODO...

# HG38 Human Reference
curl -L -o hg38.fa ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
```

### With Singularity

1. Downloading pre-packaged GenomeChronicler from Github Packages

```bash
singularity pull docker://ghcr.io/pgp-uk/genomechronicler:latest
```

2. Converting data to BAM format

```bash
singularity exec GenomeChronicler_latest.sif samtools view -@ 8 -T hg38.fa -b -o NA12878wxs.bam na12878wxs.cram
```

Running GenomeChronicler on the data

```bash
singularity run GenomeChronicler_latest.sif --bamFile=NA12878wxs.bam --GATKthreads 8
```

> Note: you can use `singularity shell` to open an interactive linux to run commands directly

### With Docker

1. Downloading pre-packaged GenomeChronicler from Github Packages

```bash
docker pull ghcr.io/pgp-uk/genomechronicler:latest
```

2. Converting data to BAM format

```bash
docker run -v $PWD:/d/ genomechronicler samtools view -@ 8 -T /d/hg38.fa -b -o /d/NA12878wxs.bam /d/na12878wxs.cram
```

Running GenomeChronicler on the data

```bash
# With a GVCF
docker run -v $PWD:/d/ genomechronicler genomechronicler --vcfFile=/GenomeChronicler/NA12878wxs.g.vcf --resultsDir /GenomeChronicler/out2 --GATKthreads 8

# With a BAM
docker run -v $PWD:/d/ genomechronicler genomechronicler --bamFile=/d/NA12878wxs.bam --resultsDir /d/out --GATKthreads 8
```

> Note: you can use `docker run --rm -it -v $PWD:/d/ genomechronicler bash` to open an interactive linux to run commands directly.

## Command Line Options

|       Option      | Requirement | Description                                                                                                                                                                                                                                                                                                                                                                                                                      |
|:-----------------:|:------------:|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| --bamFile         | REQUIRED (if no gVCF)     | The path to a BAM file that has been preprocessed through markDuplicates and VariantQualityScoreRecalibration. This can be obtained by running the first step of the Sarek nextflow pipeline, or through other means that do respect the general principles of the GATK Variation Calling Best Practices workflow. Note that no variation calling is needed to run GenomeChronicler.
| --vcfFile         | REQUIRED (if no BAM)   | The path to a gVCF file produced by GATK or the GRCh38 reference genome. This can be obtained by running all steps from the Sarek nextflow pipeline, or through other means that do respect the general principles of the GATK Variation Calling Best Practices workflow. This avoids the need to run the GATK with GC.                                             |
| --vepFile         | OPTIONAL     | For the summary tables to appear in the report, a VEP summary HTML file must be provided. This will likely be generated if the data is from whole genome sequencing and variants were called (e.g. by running all the germline calling steps of the Sarek nextflow pipeline or other GATK Best Practices based workflow). If this isn't provided, summary tables and plots will automatically be excluded from the final report. |
| --resultsDir      | OPTIONAL     | For setting the absolute path of the results folder to be produced when running GenomeChronicler.                                                                                                                                                                                                                                                                                                                                |
| --customTemplate  | OPTIONAL     | For customising the output report, set this variable to the path of a custom LaTeX file to act as a template for the report. The default templates bundled with this software can also be found in the project github page.                                                                                                                                                                                                      |
| --GATKthreads     | OPTIONAL     | Number of threads to use for the GATK genotyping steps of this processing pipeline.                                                                                                                                                                                                                                                                                                                                              |

## Development

### Working with a Docker Image

There is a separate docker image that is designed for development:

1. Clone the repository:

2. Build the development docker image

```bash
docker build -f Dockerfile-dev -t gc_dev .
```

3. Start docker environment and map local directory inside.

```bash
docker run --rm -it -v $PWD:/GenomeChronicler gc_dev bash
```

4. Now any changes you make to the codebase in your host machine, will be immediately reflected in the docker image.

### Release a new Docker build

1. Login using a token with access to Github Packages. See [here](https://docs.github.com/en/packages/working-with-a-github-packages-registry/working-with-the-container-registry) for more info.

```bash
export CR_PAT=YOUR_TOKEN

echo $CR_PAT | docker login ghcr.io -u USERNAME --password-stdin
```

2. Build docker image and upload to Github Packages

```bash
docker build -t ghcr.io/pgp-uk/genomechronicler:latest .
docker build -t ghcr.io/pgp-uk/genomechronicler:v0.0.3 .

docker push ghcr.io/pgp-uk/genomechronicler:latest
docker push ghcr.io/pgp-uk/genomechronicler:v0.0.3
```
