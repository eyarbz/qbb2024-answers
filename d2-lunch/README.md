# day 2 lunch

## Answer 1

- `cut -f 7 hg38-gene-metadata-feature.tsv | sort | uniq -c`

- There are 19618 protein_coding. Ribozyme looks interesting because catalytic RNA is cool.

## Answer 2

- `cut  -f 1  hg38-gene-metadata-go.tsv | uniq -c | sort -n`

- ENSG00000168036 has the most, with 273

- `grep ENSG00000168036 hg38-gene-metadata-go.tsv| - sort -k3 -f  | > hg38-gene-metadata-go-subset.tsv`

- Cell fate decisions

## Answer 1

- `grep  "IG_._gene" genes.gtf | cut -f 1 | uniq -c`

- Most are on chromosome 14, 91 genes.

- `grep "IG.*_pseudogene" genes.gtf | cut -f 1 | uniq -c`

- Less chromosomes have pseudogenes compared to normal Ig. Both distributions have chr14 and chr2 as their highest

## Answer 2

- Because it will pull out lines that have a substring of pseudogene in the tag key value pair, not just in the gene_type. This would be better:

- `grep "gene_type.*pseudogene.*tag" genes.gtf`

## Answer 3

- `cut -f1,4,5,14 gene-tabs.gtf > genes.bed`


