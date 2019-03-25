# long-reads

To test run the pipeline use the command:
```
nextflow run main.nf --reads s3://giab/data/AshkenazimTrio/HG002_NA24385_son/CORNELL_Oxford_Nanopore/giab.hg002.2D.fastq --fasta s3://deepvariant-data/genomes/hg19/hg19.fa --model ont-learningRate1e-4-epoch1499 --min_support 3
```
(Don't use the data in the `testdata` folder as `sniffles` will fail due to too few reads being detected in the BAM file)

The models for the pipeline can be found here: `s3://lifebit-featured-datasets/pipelines/long-reads/trainedModels/`
