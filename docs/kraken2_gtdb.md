# Preparation of a Kraken2 database with GTDB data

## 1. Set Up the Conda Environment

Create a new conda environment with the required tools:

```bash
conda create -n meta_index -c conda-forge -c bioconda mash biopython
```

## 2. Prepare GTDB Files for Kraken2 (Using rrwick's Metagenomics-Index-Correction Tools)

### 2.1 Identify the Genome Files

Locate the folder containing genomes in a recent release of the GTDB database.
In release 220, the genomes are located in `skani/database`. Look for folders named GCF and GCA, and create a variable named `SEQ`: 
```bash
SEQ='/absolute/path/to/fastani/database'
```

### 2.2 Locate the Taxonomy File

Locate a file named `gtdb_taxonomy.tsv` which is a concatenation of the archea and bacteria taxonomy .tsv files. In release 220, it is located in the `taxonomy` folder. Create another variable called `TAX`:
```bash
TAX='/absolute/path/to/taxonomy/gtdb_taxonomy.tsv'
```

### 2.3 Generate Kraken2-Compatible Taxonomy Files

To build a Kraken2 database, NCBI-style taxonomy files are generated using the GTDB definitions:

```bash
# Git clone rrwick tool
git clone https://github.com/rrwick/Metagenomics-Index-Correction.git

# Generate taxonomy files
cd Metagenomics-Index-Correction/

conda activate meta_index

./tax_from_gtdb.py \
  --gtdb $TAX \
  --assemblies $SEQ \
  --nodes nodes.dmp \
  --names names.dmp \
  --kraken_dir kraken_genomes
```
To run on an HPC using Slurm:

```bash
sbatch \
  -D $PWD \
  --output $PWD/kraken2_gtdb-%j.out \
  --export=ALL \
  -c 2 \
  -p your_partition \
  --account=your_account \
  -t 300 \
  --wrap="./tax_from_gtdb.py \
    --gtdb $TAX \
    --assemblies $SEQ \
    --nodes nodes.dmp \
    --names names.dmp \
    --kraken_dir kraken_genomes"
```

This will generate three outputs:
- `kraken_genomes`
- `names.dmp`
- `nodes.dmp`

## 3. Build the Custom Kraken2 Database
### 3.1 Set Up the Database Directory
Create an empty folder where you want your kraken db to be located. Specify it in the variable DBNAME.

```bash
DBNAME='/absolute/path/to/kraken2_gtdb'
```

You don't need to download any taxonomy because we will use the GTDB one that was created in the previous steps. Create the directory and move taxonomy files:

```bash
mkdir -p $DBNAME/taxonomy
mv names.dmp nodes.dmp $DBNAME/taxonomy
```

### 3.2 Add Sequences to the Kraken2 Library

Adding all sequences may take more than a day. While in the `Metagenomics-Index-Correction` directory:

```bash
conda activate kraken2
for f in kraken_genomes/*.fa; do 
  kraken2-build --add-to-library $f --db $DBNAME
done
```
### 3.3 Final Build Step

Building the database may also take over 24 hours:

```bash
conda activate kraken2
sbatch \
  -D $PWD \
  --output $PWD/kraken2_gtdb-%j.out \
  --export=ALL \
  -c 40 \
  -p your_partition \
  --account=your_account \
  -t 1800 \
  --wrap="kraken2-build --build --threads 40 --db $DBNAME"
```
Expected output:

```text
Creating sequence ID to taxonomy ID map (step 1)...
Sequence ID to taxonomy ID map already present, skipping map creation.
Estimating required capacity (step 2)...
Estimated hash table requirement: 415553606216 bytes
Capacity estimation complete. [1h11m28.949s]
Building database files (step 3)...
Taxonomy parsed and converted.
CHT created with 17 bits reserved for taxid.
Completed processing of 11987232 sequences, 273877474788 bp
Writing data to disk...  complete.
Database files completed. [26h42m33.513s]
Database construction complete. [Total: 28h1m3.773s]
```
 
## 4. Optional: Build a Confidence-Level Specific Bracken Database

Kraken2 databases default to a confidence threshold of 0. To generate Bracken databases for higher thresholds (e.g., 0.5), follow the steps below. For more context, see [this issue](https://github.com/jenniferlu717/Bracken/issues/154) and [Bracken's documentation](https://github.com/jenniferlu717/Bracken#step-1a-search-all-library-input-sequences-against-the-database).

### 4.1 Concatenate Library Files
In your Kraken2 DB folder:

```bash
sbatch \
  -D $PWD \
  --output $PWD/kraken2_gtdb-%j.out \
  --export=ALL \
  -c 2 \
  -p your_partition \
  --account=your_account \
  -t 100 \
  --wrap="find -L library/ -name "*.fna" | xargs cat > input.fasta"
```

### 4.2 Generate custom confidence-level db

On a node with 1TB of memory:
```bash
DBNAME='/absolute/path/to/kraken2/gtdb/release220'
conda activate kraken2
sbatch \
  -D $PWD \
  --output $PWD/kraken2_gtdb-%j.out \
  --export=ALL \
  -c 40 \
  -p your_partition \
  --account=your_account \
  -t 100 \
  --wrap="kraken2 \
    --db=$DBNAME \
    --threads=40 \
    --confidence=0.5 \
    input.fasta  > database.confidence0.5.kraken"
```

Expected output:
```shell
Loading database information... done.
11986902 sequences (273864.53 Mbp) processed in 3097.882s (232.2 Kseq/m, 5304.23 Mbp/m).
  11971677 sequences classified (99.87%)
  15225 sequences unclassified (0.13%)
```

### 4.3 Generate Read Classifications

Ensure the file `seqid2taxid.map` exists. If needed, recreate it using the `make_seqid2taxid_map.pl` script from Kraken2.

On a HPC with Slurm  with 512GB of memory:

```bash
conda activate kraken2
DBNAME='/absolute/path/to/kraken2/gtdb/release220'
THREADS=20
KMER_LEN=35
READ_LEN=150
sbatch \
  -D $PWD \
  --output $PWD/kraken2_gtdb-%j.out \
  --export=ALL \
  -c 20 \
  -t 250 \
  -p your_partition \
  --account=your_account \
  --wrap="kmer2read_distr \
    --seqid2taxid ${DBNAME}/seqid2taxid.map \
    --taxonomy ${DBNAME}/taxonomy \
    --kraken database.confidence0.5.kraken \
    --output database${READ_LEN}mers.confidence0.5.kraken \
    -k ${KMER_LEN} \
    -l ${READ_LEN} \
    -t ${THREADS}"
```
Expected output:

```shell
>>STEP 0: PARSING COMMAND LINE ARGUMENTS
        Taxonomy nodes file: /path/to/kraken2_gtdb/release220/taxonomy/nodes.dmp
        Seqid file:          /path/to/kraken2_gtdb/release220/seqid2taxid.map
        Num Threads:         20
        Kmer Length:         35
        Read Length:         150
>>STEP 1: READING SEQID2TAXID MAP
        11986902 total sequences read
>>STEP 2: READING NODES.DMP FILE
        113318 total nodes read
>>STEP 3: CONVERTING KMER MAPPINGS INTO READ CLASSIFICATIONS:
        150mers, with a database built using 35mers
        11986923 sequences converted (finished: GCF_009769735.1_NZ_WSRO01000002.1|kraken:taxid|98704))9))
Time Elaped: 212 minutes, 6 seconds, 0.00000 microseconds        
```

### 4.4 Generate the kmer distribution file

Note that the input file for this step contains the string *mers*.

On a HPC with Slurm  with 512GB of memory:

```bash
conda activate kraken2
READ_LEN=150
sbatch \
  -D $PWD \
  --output $PWD/kraken2_gtdb-%j.out \
  --export=ALL \
  -c 20 \
  -t 100 \
  -p your_partition \
  --account=your_account \
  --wrap="generate_kmer_distribution.py \
    -i database${READ_LEN}mers.confidence0.5.kraken \
    -o database${READ_LEN}mers.confidence0.5.kmer_distrib"
```
Expected output:

```shell
PROGRAM START TIME: 10-22-2023 14:09:44
...1 total genomes read from kraken output file
...creating kmer counts file -- lists the number of kmers of each classification per genome
...creating kmer distribution file -- lists genomes and kmer counts contributing to each genome
PROGRAM END TIME: 10-22-2023 14:10:19
```
### 4.5 Rename the kmer distribution file

Rename the distribution file to the expected Bracken format:
```bash
mv database150mers.confidence0.5.kmer_distrib database150mers.kmer_distrib
```


## 5. Optional Cleanup

Kraken2 provides a way to delete intermediate files, but as you may need to generate several confidence-level specific databases, it is recommended to keep the intermediate files.  

```bash
kraken2-build --clean --db $DBNAME
```

Expected output:

```bash
kraken2-build --clean --db $DBNAME
Database disk usage: 778G
After cleaning, database uses 444G
```

