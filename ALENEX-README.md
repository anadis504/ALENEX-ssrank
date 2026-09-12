# ALENEX Artifact Evaluation

This artifact contains the implementation and experimental data for the ALENEX paper `New space-time tradeoffs for subset rank and k-mer lookup`.

The artifact evaluation focuses on three reproducibility checks:

*  Construction of the SBWT variants and their resulting sizes.
*  Search performance on the variants using positive query sets.
*  s-rank query performance on the variants.


The experimental data archive contains three preconstructed plain-matrix SBWT indexes. These indexes are themselves variants and should also be included when evaluating search and s-rank functionality.

The experiments do not require the original genomic datasets; the preconstructed indexes and query data needed for evaluation are provided in the experimental-data archive.

## Requirements

The easiest way to run the artifact is using Docker.

The experiments require:

- Docker
- Sufficient disk space for the provided indexes and generated variants
- Sufficient memory for the selected dataset and experiment

The absolute runtime depends on the hardware used for evaluation. The resulting index sizes and query results should be reproducible independently of hardware.

## Obtain the experimental data

The experimental data are provided separately from the GitHub repository:

[dropbox-archive](https://tinyurl.com/ALENEX-2027-AE-Diseth-Puglisi)

Download and extract the archive.

The extracted directory contains the three plain-matrix SBWT indexes, positive-query datasets, and s-rank query data described below. 

* The three provided SBWT indexes are:

| Dataset | File | Size (bytes) |
|---|---|---:|
| *E. coli* | `Ecoli_31.sbwt` | 225,872,080 |
| 661k *Salmonella* | `661kSalmonella_31.sbwt` | 419,670,120 |
| Human | `Human_31.sbwt` | 3,990,246,312 |

The variants described in the paper can be constructed from these indexes using the `build-variant` command documented in the main repository README.


 * The positive query files contain 50,000 positive query reads of length 1,000 bp sampled from the corresponding datasets:

```text
E. coli:	        Ecoli_50k_pos_queries_len1k.fasta
661k Salmonella:	661kSalmonella_50k_pos_queries_len1k.fasta
Human:          	Human_50k_pos_queries_len1k.fasta
```

* The archive also contains the data required for the s-rank experiments.
These files contain the position and character queries used for the
s-rank experiments:

```text
test_indices.bin
test_characters.bin
```

After extraction, choose a directory for the experimental data, referred to below as:

```/path/to/data```

## Build the Docker image

Clone the repository together with its submodules:

```bash
git clone --recurse-submodules https://github.com/anadis504/ALENEX-ssrank.git
cd ALENEX-ssrank
```


Build the Docker image:

```bash
docker build -t sbwt .
```

Verify the installation:

```bash
docker run --rm sbwt --help
```

The relevant commands are:
```bash
sbwt build
sbwt build-variant
sbwt search
sbwt s-rank-queries
```

The Docker image uses `sbwt` as its entrypoint. Input and output data are supplied by mounting the experimental-data directory into the container.

For example:
```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    --help
```

## Reproducibility check 1: Constructing variants

The first reproducibility check is the construction of the SBWT variants from the supplied plain-matrix indexes.

The archive contains one plain-matrix index for each dataset:

```text
Ecoli_31.sbwt
661kSalmonella_31.sbwt
Human_31.sbwt
```

For example, to construct the correction-sets variant for E. coli:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    build-variant \
    -i /data/Ecoli_31.sbwt \
    -o /data/Ecoli_correction-sets.sbwt \
    --variant correction-sets
```

The same command can be used for the other datasets by changing the input and output filenames.

The available variant names are:
```text
plain-matrix
mef-matrix
plain-split
mef-split
ef-split
pred8-WT-split
pino-WT-split
pred8-split-packed
pred8-split-w-packed
pred8-split-transposed
pino-split-transposed
correction-sets
blocked-correction-sets
fixed-block-correction-setsA-smaller
fixed-block-correction-setsA
fixed-block-correction-setsB
fixed-block-correction-setsC
blocked8-split
blocked9-split
```


The complete list is also available with:

```bash
docker run --rm sbwt build --help
```

## Result to record

The construction command reports the size of the resulting variant.

Record the reported size for each constructed variant. These sizes are the datapoints to compare with the corresponding results in the paper.

The three supplied plain-matrix indexes are already constructed and should be treated as the plain-matrix baseline variants.

## Reproducibility check 2: Positive-query search

The second reproducibility check measures the query performance of the different SBWT variants.

The archive contains 50,000 positive query reads of length 1,000 bp for each dataset:

```text
Ecoli_50k_pos_queries_len1k.fasta
661kSalmonella_50k_pos_queries_len1k.fasta
Human_50k_pos_queries_len1k.fasta
```

For example, to search the E. coli plain-matrix index:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    search \
    -i /data/Ecoli_31.sbwt \
    -q /data/Ecoli_50k_pos_queries_len1k.fasta \
    -o /data/Ecoli_plain-matrix_search.txt
```

To search a newly constructed variant:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    search \
    -i /data/Ecoli_correction-sets.sbwt \
    -q /data/Ecoli_50k_pos_queries_len1k.fasta \
    -o /data/Ecoli_correction-sets_search.txt
```

Repeat for the desired variants and datasets.

The provided plain-matrix indexes should be included as baseline variants.

### Result to record

The primary result of this experiment is the query performance, reported by the program in microseconds per query (us/query).

The program reports two timing measurements. For example:
```text
us/query: 0.263166 (excluding I/O etc)
us/query end-to-end: 0.316309
```

For reproducing the performance results reported in the paper, record the value following us/query: and not the us/query end-to-end: value.

The us/query end-to-end value includes additional I/O and other overhead and is not the primary performance datapoint.

The output file produced with -o is not the primary evaluation result, but can be used to verify that the search completed successfully.

Record the us/query value for each selected variant and compare it with the corresponding performance datapoint in the paper.

Evaluators may run all variants or a representative subset if evaluating every variant is too time-consuming.

## Reproducibility check 3: s-rank queries

The third reproducibility check measures the performance of s-rank queries on the supplied plain-matrix indexes and on the constructed variants.

The archive contains:

```text
test_indices.bin
test_characters.bin
```

These files contain the position and character queries used for the s-rank experiments and can be used with all three datasets and their corresponding variants.

For example, to run the s-rank queries on the E. coli plain-matrix index:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    s-rank-queries \
    -i /data/Ecoli_31.sbwt \
    -p /data/test_indices.bin \
    -c /data/test_characters.bin \
    -o /data/Ecoli_plain-matrix_s-rank.txt
```

To run the same queries on a constructed variant:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    s-rank-queries \
    -i /data/Ecoli_correction-sets.sbwt \
    -p /data/test_indices.bin \
    -c /data/test_characters.bin \
    -o /data/Ecoli_correction-sets_s-rank.txt
```

Repeat for the desired variants and datasets.

The supplied plain-matrix indexes should be included as baseline variants.

### Result to record

The primary result of this experiment is the s-rank query performance reported by the s-rank-queries command, expressed as time per query.

The reported query time excludes input/output time and therefore measures the performance of the s-rank operation itself.

Record the reported time per query for each variant and compare it with the corresponding datapoint in the paper.

The output file produced with -o can be used as a sanity check that the queries completed successfully, but the primary result is the reported query performance.

Evaluators may run all variants or a representative subset if evaluating every variant is too time-consuming.
