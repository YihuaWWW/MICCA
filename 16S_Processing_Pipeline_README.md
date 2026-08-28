# 16S rRNA Amplicon Data Processing Pipeline

This document provides a practical workflow for processing publicly available 16S rRNA amplicon sequencing data from SRA/ENA downloads to OTU tables, taxonomic assignments, and QIIME 2 diversity analyses.

The workflow is organized into the following stages:

1. Download raw sequencing data
2. Convert SRA/ENA files to FASTQ
3. Organize sequencing files
4. Merge paired-end reads
5. Remove primers
6. Assess sequencing quality
7. Quality filtering
8. Dereplication and OTU/ASV processing
9. Chimera removal
10. Generate an OTU table
11. Assign taxonomy
12. Prepare a BIOM table
13. Import data into QIIME 2
14. Build a phylogenetic tree
15. Perform alpha- and beta-diversity analyses
16. Assess sequencing depth and alpha-diversity differences

> **Important:** Replace all paths, filenames, sample IDs, database locations, and sequencing-specific parameters with values appropriate for your own dataset. The commands below reflect the original workflow and should be validated against your software versions and experimental design before use.

---

## 1. Required Software and Databases

The original workflow uses:

- SRA Toolkit (`prefetch`, `fasterq-dump`)
- `gzip`
- VSEARCH
- USEARCH
- QIIME 1
- QIIME 2
- BIOM
- Dendroscope
- RDP/SILVA reference databases

Recommended environment setup should be performed according to the software versions available on your server.

### Example QIIME 1 environment

```bash
conda activate qiime1
```

### Example QIIME 2 environment

```bash
conda activate qiime2-2023.7
```

---

# Part I. Download and Prepare Raw Data

## 2. Download SRA/ENA Data

### 2.1 Prepare accession lists

Create one text file for each dataset or cohort, for example:

```text
No2.txt
No17.txt
No19.txt
No24.txt
No28.txt
No33.txt
No47.txt
No71.txt
No74.txt
No75.txt
```

Each file should contain the corresponding accession IDs required by `prefetch`.

### 2.2 Download the data

Run:

```bash
prefetch --option-file No2.txt
prefetch --option-file No17.txt
prefetch --option-file No19.txt
prefetch --option-file No24.txt
prefetch --option-file No28.txt
prefetch --option-file No33.txt
prefetch --option-file No47.txt
prefetch --option-file No71.txt
prefetch --option-file No74.txt
prefetch --option-file No75.txt
```

For a large project, it is usually more convenient to maintain a single accession list and process it systematically.

---

# Part II. Convert SRA/ENA Files to FASTQ

## 3. Convert SRA Files to FASTQ

### 3.1 Convert paired-end SRA files

For directories containing `SRR` accessions:

```bash
for id in SRR*; do
    nohup fasterq-dump \
        -O ./ \
        --split-files \
        -e 2 \
        "./${id}" \
        --include-technical &
done
```

For `ERR` accessions:

```bash
for id in ERR*; do
    nohup fasterq-dump \
        -O ./ \
        --split-files \
        -e 2 \
        "./${id}" \
        --include-technical &
done
```

### 3.2 Compress FASTQ files

```bash
for file in *.fastq; do
    nohup gzip "$file" &
done
```

After compression, paired-end files should typically appear as:

```text
Sample1_1.fastq.gz
Sample1_2.fastq.gz
```

or, before compression:

```text
Sample1_1.fastq
Sample1_2.fastq
```

---

# Part III. Organize the Project

## 4. Recommended Directory Structure

A simple project structure is:

```text
project/
├── raw/
├── seq/
├── temp/
├── result/
├── metadata/
└── db/
```

Where:

- `raw/` — downloaded SRA/ENA files
- `seq/` — paired-end FASTQ files
- `temp/` — intermediate files
- `result/` — final analysis outputs
- `metadata/` — sample metadata
- `db/` — reference databases

Create the main directories with:

```bash
mkdir -p seq temp result metadata db
```

### 4.1 Copy sequencing data

Replace the example paths with your own:

```bash
cp -r /path/to/rawdata/* /path/to/project/
```

If several cohorts are stored separately:

```bash
cp -r /path/to/rawdata/cohort1/* /path/to/project/
cp -r /path/to/rawdata/cohort2/* /path/to/project/
```

### 4.2 Move paired-end reads into the `seq/` directory

```bash
mv /path/to/project/*_1.fastq /path/to/project/seq/
mv /path/to/project/*_2.fastq /path/to/project/seq/
```

If the files are compressed:

```bash
mv /path/to/project/*_1.fastq.gz /path/to/project/seq/
mv /path/to/project/*_2.fastq.gz /path/to/project/seq/
```

---

# Part IV. Merge Paired-End Reads

## 5. Merge Forward and Reverse Reads

Prepare a metadata file such as:

```text
metadata_all.txt
```

The first column should contain sample IDs corresponding to the FASTQ filenames.

For example:

```text
SampleID    Group    Cohort
Sample01    Control  Cohort1
Sample02    Disease  Cohort1
```

### 5.1 Merge paired-end reads with VSEARCH

```bash
for sample in $(tail -n +2 metadata_all.txt | cut -f 1); do
    vsearch \
        --fastq_mergepairs "seq/${sample}_1.fastq" \
        --reverse "seq/${sample}_2.fastq" \
        --fastqout "temp/${sample}.merged.fq"
done
```

If your FASTQ files are located directly in the working directory:

```bash
for sample in $(tail -n +2 metadata_all.txt | cut -f 1); do
    vsearch \
        --fastq_mergepairs "${sample}_1.fastq" \
        --reverse "${sample}_2.fastq" \
        --fastqout "temp/${sample}.merged.fq"
done
```

> If your files are compressed, check the VSEARCH version and supported input format before running the command.

---

# Part V. Combine Samples

## 6. Concatenate Merged Reads

After all samples have been merged:

```bash
cat temp/*.merged.fq > temp/all.fq
```

Check the resulting file:

```bash
ls -lh temp/all.fq
```

The original workflow also used several intermediate concatenation commands. In most cases, one consistently named merged FASTQ file is easier to track and reduces the risk of accidentally combining the same reads multiple times.

---

# Part VI. Remove Primers

## 7. Remove Primer Sequences

If the forward primer length is 19 bp and the reverse primer length is 20 bp, the original workflow used:

```bash
usearch \
    -fastx_truncate temp/all.fq \
    -stripleft 19 \
    -stripright 20 \
    -fastqout temp/stripped.fastq
```

### Important

The values `19` and `20` are **dataset-specific**. They should correspond to the actual primer lengths in your sequencing protocol.

Do not use these values automatically for a different primer set.

---

# Part VII. Assess Sequence Quality

## 8. Basic FASTQ Quality Summary

```bash
usearch \
    -fastx_info temp/stripped.fastq \
    -output temp/allreads_info.txt
```

This provides a summary of the sequence dataset after primer removal.

---

## 9. Estimate Expected Errors

Before selecting a maximum expected error threshold (`maxee`), inspect the expected-error distribution:

```bash
usearch \
    -fastq_eestats2 temp/stripped.fastq \
    -fastq_ascii 33 \
    -fastq_qmin 0 \
    -fastq_qmax 47 \
    -output temp/eestats_all.txt \
    -ee_cutoffs 0.6,0.8,1.0 \
    -length_cutoffs 350,390,10
```

The expected-error threshold should be selected based on the quality profile of the dataset.

---

# Part VIII. Quality Filtering

## 10. Filter Low-Quality Reads

The original workflow used `maxee = 1.0`:

```bash
usearch \
    -fastq_filter temp/stripped.fastq \
    -fastq_ascii 33 \
    -fastq_qmin 0 \
    -fastq_qmax 47 \
    -fastq_maxee 1.0 \
    -fastaout temp/filtered.fa
```

An equivalent VSEARCH-based approach used:

```bash
vsearch \
    --fastq_filter temp/stripped.fq \
    --fastq_maxee 1 \
    --fastq_maxlen 553 \
    --fastq_minlen 200 \
    --fastaout temp/filtered.fa \
    --fastq_qmax 42
```

A simpler VSEARCH version is:

```bash
vsearch \
    --fastq_filter temp/stripped.fq \
    --fastq_maxee 1 \
    --fastaout temp/filtered.fa \
    --fastq_qmax 42
```

### Recommended practice

Use **one filtering strategy consistently** throughout the project. Do not run several alternative filtering commands sequentially unless you intentionally want to compare their effects.

---

# Part IX. Dereplication

## 11. Dereplicate Full-Length Sequences

### USEARCH

```bash
usearch \
    --fastx_uniques temp/filtered.fa \
    --fastaout temp/uniques.fa \
    -sizeout \
    -relabel Uniq
```

### VSEARCH

```bash
vsearch \
    --derep_fulllength temp/filtered.fa \
    --output temp/uniques.fa \
    --sizeout \
    --relabel Uniq
```

The original workflow retained unique sequences with a minimum abundance of 6 during downstream processing.

---

# Part X. OTU/ASV Processing

## 12. Denoising / ASV Generation

The original workflow used the USEARCH UNOISE3 algorithm:

```bash
usearch \
    -unoise3 temp/uniques.fa \
    -zotus temp/zotus.fa
```

This produces denoised representative sequences.

> The original pipeline subsequently used the filename `otus.fa` in the chimera-removal step. Make sure that the output filename is consistent. For example, rename the UNOISE output explicitly if you intend to use it as `otus.fa`:

```bash
cp temp/zotus.fa temp/otus.fa
```

---

## 13. Alternative: 97% OTU Clustering

The original workflow also included a conventional OTU clustering approach:

```bash
usearch \
    -cluster_otus temp/uniques.fa \
    -otus result/otus.fa \
    -uparseout result/otu_results.txt \
    -relabel OTU
```

> **Choose the strategy that matches your study design.** UNOISE3 generates denoised sequence variants, whereas 97% clustering produces conventional OTUs. These approaches should not be mixed without a clear methodological reason.

---

# Part XI. Chimera Removal

## 14. Reference-Based Chimera Detection

Using the RDP reference database:

```bash
usearch \
    -uchime2_ref temp/otus.fa \
    -db db/rdp_16s_v16_sp.fa \
    -chimeras temp/otus_chimeras.fa \
    -strand plus \
    -mode balanced
```

Identify non-chimeric sequences:

```bash
cat temp/otus.fa temp/otus_chimeras.fa |
    grep '>' |
    sort |
    uniq -u |
    sed 's/>//' \
    > temp/non_chimeras.id
```

Extract the non-chimeric representative sequences:

```bash
usearch \
    -fastx_getseqs temp/otus.fa \
    -labels temp/non_chimeras.id \
    -fastaout result/otus.fa
```

---

# Part XII. Generate the OTU Table

## 15. Map Reads Back to Representative Sequences

Using VSEARCH:

```bash
vsearch \
    --usearch_global temp/filtered.fa \
    --db result/otus.fa \
    --otutabout result/otutab.txt \
    --id 0.97 \
    --threads 4
```

Alternatively, using USEARCH:

```bash
usearch \
    -otutab temp/filtered.fa \
    -otus result/otus.fa \
    -otutabout result/otutable.txt \
    -mapout result/map.txt
```

### Output

The main abundance table is:

```text
result/otutab.txt
```

or:

```text
result/otutable.txt
```

Use one filename consistently in downstream analyses.

---

# Part XIII. Taxonomic Assignment

## 16. Taxonomy Assignment with RDP

Assign taxonomy using the RDP database:

```bash
usearch \
    -sintax result/otus.fa \
    -db db/rdp_16s_v16_sp.fa \
    -strand both \
    -tabbedout temp/sintax.txt \
    -sintax_cutoff 0.6
```

The original workflow also tested:

```bash
usearch \
    -sintax result/otus.fa \
    -db db/rdp_gold.fa \
    -strand both \
    -tabbedout temp/sintax_rdp214.txt \
    -sintax_cutoff 0.6
```

---

## 17. Taxonomy Assignment with SILVA

```bash
usearch \
    -sintax result/otus.fa \
    -db db/silva_species_assignment_v138.1.fa \
    -strand both \
    -tabbedout temp/sintax_silva.txt \
    -sintax_cutoff 0.6
```

> Use the reference database that is appropriate for your study and document its version in the final analysis report.

---

# Part XIV. Format Taxonomy Tables

## 18. Create a Two-Column Taxonomy Table

The original workflow converted the SINTAX output into a two-column table:

```bash
cut -f 1,4 temp/sintax.txt |
sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
> result/taxonomy2.txt
```

The output contains:

```text
OTUID    Taxonomy
```

---

## 19. Create an Eight-Column Taxonomy Table

The desired format is:

```text
OTUID
Kingdom
Phylum
Class
Order
Family
Genus
Species
```

The original workflow used the following command:

```bash
awk 'BEGIN{OFS=FS="\t"}{
    delete a;
    a["k"]="Unassigned";
    a["p"]="Unassigned";
    a["c"]="Unassigned";
    a["o"]="Unassigned";
    a["f"]="Unassigned";
    a["g"]="Unassigned";
    a["s"]="Unassigned";
    split($2,x,";");
    for(i in x){
        split(x[i],b,"__");
        a[b[1]]=b[2];
    }
    print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];
}' result/taxonomy2.txt > temp/otus.tax
```

Then format the table:

```bash
sed 's/;/\t/g;s/.__//g;' temp/otus.tax |
cut -f 1-8 |
sed '1s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
> result/taxonomy.txt
```

The final taxonomy table is:

```text
result/taxonomy.txt
```

---

# Part XV. Prepare a BIOM Table

## 20. Combine the OTU Table and Taxonomy Table

The original workflow used Excel to combine:

```text
otutable.txt
taxonomy.txt
```

If the OTU identifiers in the two tables do not match exactly, use the **intersection of OTUs**.

The final table should contain:

- OTU abundance information
- Corresponding taxonomy
- Consistent OTU identifiers

Save the final table as:

```text
otutable_tax.txt
```

with tab delimiters.

> For reproducible analyses, a scripted merge using R, Python, or command-line tools is preferable to manual Excel operations.

---

## 21. Convert the OTU Table to BIOM

Activate the QIIME 1 environment:

```bash
conda activate qiime1
```

Then:

```bash
biom convert \
    -i otutable_tax.txt \
    -o otutable_tax.biom \
    --to-hdf5 \
    --table-type="OTU table" \
    --process-obs-metadata taxonomy
```

---

# Part XVI. QIIME 1 Taxonomic Summaries

## 22. Prepare Sample Metadata

Prepare:

```text
sample_metadata.tsv
```

The metadata file should contain sample IDs matching the IDs in the OTU table.

---

## 23. Summarize Counts at Different Taxonomic Levels

```bash
summarize_taxa.py \
    -i otutable_tax.biom \
    -o ./Sumtax_counts/ \
    -a
```

---

## 24. Summarize Relative Abundance

```bash
summarize_taxa_through_plots.py \
    -o Sumtax_reabun \
    -i otutable_tax.biom \
    -m sample_metadata.tsv
```

---

# Part XVII. QIIME 2 Analysis

## 25. Activate QIIME 2

```bash
conda activate qiime2-2023.7
```

---

## 26. Convert OTU Table to BIOM

If `otutable.txt` is already a valid BIOM-compatible OTU table:

```bash
biom convert \
    -i otutable.txt \
    -o otutable.biom \
    --to-hdf5 \
    --table-type="OTU table"
```

---

## 27. Import the Feature Table

```bash
qiime tools import \
    --input-path otutable.biom \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path otutable.qza
```

---

## 28. Summarize the Feature Table

Prepare:

```text
sample_metadata.tsv
```

Then:

```bash
qiime feature-table summarize \
    --i-table otutable.qza \
    --o-visualization otutable.qzv \
    --m-sample-metadata-file sample_metadata.tsv
```

View the result:

```bash
qiime tools view otutable.qzv
```

This visualization allows you to inspect:

- Number of samples
- Number of features/OTUs
- Sequencing depth
- Sample-specific feature counts

---

# Part XVIII. Import Representative Sequences

## 29. Prepare Representative Sequences

Use the final representative sequence file:

```text
otus.fa
```

If archaeal sequences are required for rooting the phylogenetic tree, the original workflow combined the representative sequences with an archaeal reference file:

```text
Archaea.fa
```

Create a combined FASTA file as appropriate for your study.

---

## 30. Import Representative Sequences into QIIME 2

```bash
qiime tools import \
    --input-path otus.fa \
    --output-path otus.qza \
    --type 'FeatureData[Sequence]'
```

---

## 31. Optional: Summarize Representative Sequences

```bash
qiime feature-table tabulate-seqs \
    --i-data otus.qza \
    --o-visualization otus.qzv
```

View the result:

```bash
qiime tools view otus.qzv
```

---

# Part XIX. Construct a Phylogenetic Tree

## 32. Multiple Sequence Alignment

```bash
qiime alignment mafft \
    --i-sequences otus.qza \
    --o-alignment aligned-otus.qza
```

---

## 33. Mask Highly Variable Regions

```bash
qiime alignment mask \
    --i-alignment aligned-otus.qza \
    --o-masked-alignment masked-aligned-otus.qza
```

---

## 34. Build an Unrooted Tree

```bash
qiime phylogeny fasttree \
    --i-alignment masked-aligned-otus.qza \
    --o-tree unrooted-tree.qza
```

---

## 35. Export the Tree

```bash
qiime tools export \
    --input-path unrooted-tree.qza \
    --output-path exported-tree
```

The exported Newick tree can be found in:

```text
exported-tree/
```

Typically:

```text
tree.nwk
```

---

# Part XX. Root the Phylogenetic Tree

## 36. Root the Tree Using an Archaeal Outgroup

The original workflow used Dendroscope to:

1. Open the Newick tree.
2. Locate the archaeal sequences.
3. Use the archaeal sequences as the outgroup.
4. Reroot the tree.
5. Remove the archaeal sequences after rooting.
6. Export the rooted tree in Newick format.

Save the resulting tree as:

```text
rooted-tree.tre
```

Upload it to the QIIME 2 working directory.

---

## 37. Import the Rooted Tree into QIIME 2

```bash
qiime tools import \
    --input-path rooted-tree.tre \
    --output-path rooted-tree.qza \
    --type 'Phylogeny[Rooted]'
```

### Alternative: Midpoint Rooting

The original workflow also contained:

```bash
qiime phylogeny midpoint-root \
    --i-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza
```

However, if an appropriate biological outgroup is available, the outgroup-rooted tree should be used consistently with the study design.

---

# Part XXI. Alpha and Beta Diversity

## 38. Select a Sampling Depth

Before running core diversity metrics, inspect the feature-table summary and choose an appropriate sampling depth.

The original workflow used:

```text
8279
```

This is **dataset-specific** and should not automatically be reused for another dataset.

---

## 39. Calculate Core Diversity Metrics

```bash
qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree.qza \
    --i-table otutable.qza \
    --p-sampling-depth 8279 \
    --m-metadata-file sample_metadata.tsv \
    --output-dir core-metrics-results
```

The output directory contains multiple alpha- and beta-diversity results.

---

## 40. Main Diversity Outputs

### Alpha diversity

Common outputs include:

```text
evenness_vector.qza
faith_pd_vector.qza
observed_features_vector.qza
shannon_vector.qza
```

These represent:

- **Evenness** — distributional evenness of features
- **Faith's PD** — phylogenetic diversity
- **Observed features** — number of observed features
- **Shannon diversity** — richness and evenness

### Beta diversity

Common outputs include:

```text
bray_curtis_distance_matrix.qza
jaccard_distance_matrix.qza
unweighted_unifrac_distance_matrix.qza
weighted_unifrac_distance_matrix.qza
```

These represent different measures of between-sample community dissimilarity.

---

# Part XXII. Export QIIME 2 Results

## 41. Export Diversity Results

Create a separate directory for each exported result:

```bash
cd core-metrics-results
```

For example:

```bash
qiime tools export \
    --input-path bray_curtis_distance_matrix.qza \
    --output-path exported_bray_curtis_distance_matrix

qiime tools export \
    --input-path bray_curtis_pcoa_results.qza \
    --output-path exported_bray_curtis_pcoa_results

qiime tools export \
    --input-path evenness_vector.qza \
    --output-path exported_evenness_vector

qiime tools export \
    --input-path faith_pd_vector.qza \
    --output-path exported_faith_pd_vector

qiime tools export \
    --input-path jaccard_distance_matrix.qza \
    --output-path exported_jaccard_distance_matrix

qiime tools export \
    --input-path jaccard_pcoa_results.qza \
    --output-path exported_jaccard_pcoa_results

qiime tools export \
    --input-path observed_features_vector.qza \
    --output-path exported_observed_features_vector

qiime tools export \
    --input-path shannon_vector.qza \
    --output-path exported_shannon_vector

qiime tools export \
    --input-path unweighted_unifrac_distance_matrix.qza \
    --output-path exported_unweighted_unifrac_distance_matrix

qiime tools export \
    --input-path unweighted_unifrac_pcoa_results.qza \
    --output-path exported_unweighted_unifrac_pcoa_results

qiime tools export \
    --input-path weighted_unifrac_distance_matrix.qza \
    --output-path exported_weighted_unifrac_distance_matrix

qiime tools export \
    --input-path weighted_unifrac_pcoa_results.qza \
    --output-path exported_weighted_unifrac_pcoa_results
```

---

# Part XXIII. Rarefaction Analysis

## 42. Generate Alpha-Rarefaction Curves

Return to the directory containing `otutable.qza`:

```bash
cd ..
```

Run:

```bash
qiime diversity alpha-rarefaction \
    --i-table otutable.qza \
    --i-phylogeny rooted-tree.qza \
    --p-max-depth 8279 \
    --p-min-depth 10 \
    --p-iterations 828 \
    --m-metadata-file sample_metadata.tsv \
    --o-visualization alpha-rarefaction.qzv \
    --output-dir alpha-rarefaction
```

View the result:

```bash
qiime tools view alpha-rarefaction.qzv
```

### Parameter note

The values:

```text
max-depth = 8279
min-depth = 10
iterations = 828
```

come from the original workflow and should be adapted to the size and sequencing-depth distribution of the new dataset.

---

# Part XXIV. Test Alpha-Diversity Differences

## 43. Faith's Phylogenetic Diversity

```bash
qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
    --m-metadata-file sample_metadata.tsv \
    --o-visualization core-metrics-results/faith-pd-group-significance.qzv
```

View the result:

```bash
qiime tools view \
    core-metrics-results/faith-pd-group-significance.qzv
```

---

## 44. Evenness

```bash
qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results/evenness_vector.qza \
    --m-metadata-file sample_metadata.tsv \
    --o-visualization core-metrics-results/evenness-group-significance.qzv
```

View the result:

```bash
qiime tools view \
    core-metrics-results/evenness-group-significance.qzv
```

---

# Part XXV. Recommended Final Output Structure

At the end of the workflow, organize the important outputs as follows:

```text
project/
├── metadata/
│   ├── metadata_all.txt
│   └── sample_metadata.tsv
│
├── seq/
│   ├── Sample01_1.fastq
│   ├── Sample01_2.fastq
│   └── ...
│
├── temp/
│   ├── all.fq
│   ├── stripped.fastq
│   ├── filtered.fa
│   ├── uniques.fa
│   ├── otus.fa
│   ├── sintax.txt
│   └── ...
│
├── result/
│   ├── otutab.txt
│   ├── otutable.txt
│   ├── otutable_tax.txt
│   ├── taxonomy.txt
│   └── otus.fa
│
├── db/
│   ├── rdp_16s_v16_sp.fa
│   ├── rdp_gold.fa
│   └── silva_species_assignment_v138.1.fa
│
└── qiime2/
    ├── otutable.qza
    ├── otutable.qzv
    ├── otus.qza
    ├── rooted-tree.qza
    ├── core-metrics-results/
    └── alpha-rarefaction/
```

---

# Quick Workflow Summary

The complete workflow can be summarized as:

```text
SRA/ENA accession IDs
        │
        ▼
    prefetch
        │
        ▼
   fasterq-dump
        │
        ▼
     FASTQ
        │
        ▼
 Merge paired-end reads
        │
        ▼
 Remove primers
        │
        ▼
 Quality assessment
        │
        ▼
 Quality filtering
        │
        ▼
   Dereplication
        │
        ▼
 Denoising / OTU clustering
        │
        ▼
 Chimera removal
        │
        ▼
 Representative sequences
        │
        ▼
    OTU table
        │
        ├──────────────► Taxonomic assignment
        │                       │
        │                       ▼
        │                  Taxonomy table
        │
        ▼
      BIOM
        │
        ▼
      QIIME 2
        │
        ├──────────────► Feature-table summary
        │
        ├──────────────► Phylogenetic tree
        │
        ├──────────────► Alpha diversity
        │
        ├──────────────► Beta diversity
        │
        └──────────────► Rarefaction analysis
```

---

# Important Reproducibility Notes

## 1. Keep sample IDs consistent

The sample ID should be identical across:

- FASTQ filenames
- metadata files
- OTU table
- taxonomy table
- QIIME 2 artifacts

For example:

```text
Sample01_1.fastq
Sample01_2.fastq
```

should correspond to:

```text
Sample01
```

in the metadata.

## 2. Record all dataset-specific parameters

At minimum, record:

- Primer sequences and lengths
- Read length
- Maximum expected error (`maxee`)
- Minimum and maximum sequence length
- Dereplication threshold
- OTU clustering threshold, if applicable
- Chimera-removal database and version
- Taxonomic reference database and version
- QIIME 2 version
- USEARCH version
- VSEARCH version
- Sampling depth

## 3. Avoid mixing alternative processing strategies

The original notes contain several alternative commands, including:

- USEARCH vs. VSEARCH
- UNOISE3 denoising vs. 97% OTU clustering
- Different FASTQ filtering commands
- Different taxonomy databases

These should be treated as **alternative methods**, not as consecutive steps.

For a reproducible analysis, select one strategy and apply it consistently to all samples.

## 4. Validate parameters before processing a new dataset

Values such as:

```text
--fastq_maxee 1
--fastq_minlen 200
--fastq_maxlen 553
--id 0.97
--p-sampling-depth 8279
--p-max-depth 8279
```

were used in the original workflow. They are not universal defaults and should be reconsidered for each dataset.

## 5. Preserve intermediate files

Keep intermediate files until the complete analysis has been successfully validated. This makes it possible to identify problems at individual processing stages without repeating the entire pipeline.
