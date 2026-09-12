# ALENEX Artifact Evaluation Data

This [archive](link_to_server) contains the experimental data required to reproduce the experiments in the paper *New space-time tradeoffs for subset rank and k-mer lookup*.

The three datasets used in the experiments are described in detail in the paper.

## Files

### Preconstructed SBWT indexes

These are plain-matrix SBWT indexes used as inputs for constructing the variants evaluated in the paper.

| Dataset | File | Size (bytes) |
|---|---|---:|
| *E. coli* | `Ecoli_31.sbwt` | 225,872,080 |
| 661k *Salmonella* | `661kSalmonella_31.sbwt` | 419,670,120 |
| Human | `Human_31.sbwt` | 3,990,246,312 |

The variants described in the paper can be constructed from these indexes using the `build-variant` command documented in the main repository README.

### Positive query reads

The following files contain 50,000 positive query reads of length 1,000 bp sampled from the corresponding datasets:

```text
Ecoli_50k_pos_queries_len1k.fasta
661kSalmonella_50k_pos_queries_len1k.fasta
Human_50k_pos_queries_len1k.fasta
```

### s-rank query data

### s-rank query data

The following files contain the position and character queries used for the
s-rank experiments:

```text
test_indices.bin
test_characters.bin
```

These query files can be used with all three supplied SBWT indexes and with
the SBWT variants constructed from them.

For example:

```bash
./bin/sbwt s-rank-queries \
    -i /data/Ecoli_31.sbwt \
    -p /data/test_indices.bin \
    -c /data/test_characters.bin \
    -o /data/s-rank-results.txt
```

The same command can be used with Ecoli_31.sbwt, 661kSalmonella_31.sbwt,
Human_31.sbwt, or any of the corresponding variants constructed from these
indexes.


### Reproducibility

The source code and Docker environment are available in the accompanying GitHub repository:

https://github.com/anadis504/ALENEX-ssrank

The .sbwt files are provided as preconstructed plain-matrix indexes so that the SBWT variants evaluated in the paper can be constructed directly without requiring the original datasets.