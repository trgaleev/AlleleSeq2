This script is for querying the ENTEx AS catalog.

Please enter the path to the input file, followed by at most three of the following types of keywords:

 * Assay (e.g. H3K27ac)
 * Individual (e.g. ENC-001)
 * Tissue (e.g. spleen)
      
The output is all AS SNPs in the catalog from the query.

Example:
```
bash query_snp.sh hetSNVs_default_AS.tsv H3K27ac ENC-001 spleen
```

To write the output to file:
```
bash query_snp.sh hetSNVs_default_AS.tsv H3K27ac ENC-001 spleen > output.tsv
```

To count the number of SNPs:
```
bash query_snp.sh hetSNVs_default_AS.tsv H3K27ac ENC-001 spleen | wc -l
```

To print the help information:
```
bash query_snp.sh 
```
