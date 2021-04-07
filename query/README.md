This script is for querying of the ENTEx AS catalog.
Please enter the path to the input file, followed by at most three of the following types of keywords:
      Assay (e.g. H3K27ac)
      Individual (e.g. ENC-001)
      Tissue (e.g. spleen)
The output will count the number of AS SNPs in the catalog from the query.

Example:
```
bash query_snp.sh hetSNVs_default_AS.tsv H3K27ac ENC-001 Spleen
```

To print the help information:
```
bash query_snp.sh 
```
