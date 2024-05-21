[TOC]


# Iterate over each sample ID in 'metadata_dxal.tsv' file, starting from the second line
time for i in `tail -n+2 ./metadata_dxal.tsv | cut -f 1`; do
  # Use vsearch to merge forward and reverse reads of each sample
  vsearch --fastq_mergepairs seq_dxal_16s/${i}.R1.fq --reverse seq_dxal_16s/${i}.R2.fq \
  --fastqout temp/${i}.merged.fq --relabel ${i}.
done &

# Concatenate all merged fastq files into a single file
cat temp/*.fq > temp/all.fq

## Step 3: Find and remove PCR duplicates
time usearch11 \
  -search_pcr2 temp/all.fq \
  -fwdprimer AGRGTTTGATCMTGGCTCAG \
  -revprimer GWATTACCGCGGCKGCTG  \
  -minamp 180                      \
  -maxamp 480                      \
  -strand both                     \
  -maxdiffs 2                      \
  -tabbedout ./temp/qc.txt   \
  -fastqout  ./temp/qc.fastq \
  -log ./temp/qc.log

## Step 4: Dereplicate sequences
time vsearch --derep_fulllength temp/qc.fastq \
  --output temp/uniques.fa --relabel Uni --minuniquesize 10 --sizeout

## Step 5: Run UNOISE3 to denoise sequences and generate ZOTUs
time usearch11 -unoise3 temp/uniques.fa \
  -zotus temp/zotus.fa

## Step 6: Build an OTU table
time usearch11 -otutab temp/qc.fastq -otus result/raw/otus.fa \
  -otutabout result/raw/otutab.txt -threads 8

## Step 7a: Calculate distance matrix for OTUs
time usearch11 -calc_distmx result/raw/otus.fa -tabbedout TREE/distmx.txt \
  -maxdist 0.2 -termdist 0.3
  
## Step 7b: Calculate distance matrix in Phylip format
time usearch11 -calc_distmx result/raw/otus.fa -tabbedout TREE/dist.tree \
  -format phylip_lower_triangular

## Step 7c: Cluster OTUs and generate a phylogenetic tree
time usearch11 -cluster_aggd TREE/distmx.txt -treeout TREE/otu.tree \
  -clusterout TREE/clusters.txt -id 0.80 -linkage min

## Step 8: Annotate OTUs with taxonomy using SINTAX
time vsearch --sintax result/raw/otus.fa --db ${db}/usearch/silva_16s_v123.fa \
  --tabbedout result/raw/otus.sintax --sintax_cutoff 0.8
head result/raw/otus.sintax

# Copy OTU files to the final results directory
cp result/raw/otu* result/

## Step 9: Simplify and format the SINTAX taxonomy file
cut -f 1,4 result/otus.sintax \
  |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
  > result/taxonomy2.txt
head -n3 result/taxonomy2.txt

## Step 10: Parse and format the taxonomy file into a table
awk 'BEGIN{OFS=FS="\t"}{
  delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
  split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
  print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];
}' result/taxonomy2.txt > temp/otus.tax

# Format the taxonomy table and add header
sed 's/;/\t/g;s/.__//g;' temp/otus.tax | cut -f 1-8 | \
  sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
  > result/taxonomy.txt

# Display the first three lines of the formatted taxonomy table
head -n3 result/taxonomy.txt





