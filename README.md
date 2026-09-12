# SBWT

This is the code for the paper "New space-time tradeoffs for subset rank and k-mer lookup". The repository includes implementations of the various SBWT variants described in the paper. The data structures answer k-mer membership queries on the input data. Note that contrary to many other k-mer membership data structures, our code is not aware of DNA reverse complements. That is, it considers a k-mer and its reverse complement as separate k-mers.

This construction algorithm is based on the lightning-fast [k-mer counter KMC](https://github.com/refresh-bio/KMC). We call the KMC binaries directly from our code. The construction is very disk-heavy, so it is recommended to run construction code off a fast SSD drive.

# Compiling

Download the repository with its submodules

```bash
git clone --recurse-submodules https://github.com/anadis504/ALENEX-ssrank.git

cd ALENEX-ssrank

```

## Building with Docker

A Dockerfile is provided for building and running SBWT in a self-contained Ubuntu 22.04 environment. The Docker image builds the sbwt executable automatically.

Building the Docker image

From the repository root, run:

```bash
docker build -t sbwt .
```


After the image has been built, verify that the executable is available:

```bash
docker run --rm sbwt --help
```

The container exposes sbwt as its entrypoint, so Docker arguments are passed directly to the SBWT executable.

### Building an SBWT index with Docker

The repository contains example input data in example_data/. For example:
```bash
docker run --rm \
    sbwt \
    build \
    -i /SBWT/example_data/coli3.fna \
    -o /tmp/index.sbwt \
    -k 30
```

For data outside the repository, mount the directory containing the input data into the container. For example, if `/path/to/data` contains `input.fna`:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    build \
    -i /data/input.fna \
    -o /data/index.sbwt \
    -k 30
```

The directory `/data` is a directory inside the container mapped to `/path/to/data` on the host. Therefore, the resulting `index.sbwt` will be available directly in `/path/to/data` after the container exits.

Other SBWT commands can be run in the same way. For example, to search an existing index:

```bash
docker run --rm \
    -v "/path/to/data:/data" \
    sbwt \
    search \
    -i /data/index.sbwt \
    -q /data/queries.fastq \
    -o /data/out.txt
```

To see the available commands and options:

```bash
docker run --rm sbwt --help
```

To see the options for a particular command:

```bash
docker run --rm sbwt build --help
```

or:

```bash
docker run --rm sbwt search --help
```

## Manual compilation

Alternatively, SBWT can be compiled directly on Ubuntu.

```bash
apt-get update
apt-get install -y g++ gcc cmake git python3-dev g++-8 libz-dev

git clone --recurse-submodules \
    https://github.com/anadis504/ALENEX-ssrank.git

cd ALENEX-ssrank
mkdir build
cd build

cmake .. -DCMAKE_CXX_COMPILER=g++-8 -DMAX_KMER_LENGTH=32
make -j8
```

On MacOS `cmake` with the following flags:

```bash
cmake .. -DCMAKE_CXX_COMPILER=g++-12 -DMAX_KMER_LENGTH=32 -DCMAKE_EXE_LINKER_FLAGS=-Wl,-ld_classic
```


Change the parameter `-DMAX_KMER_LENGTH=32` to increase the maximum allowed k-mer length, up to 255. Larger values lead to slower construction and higher disk usage during construction.

**Troubleshooting**: If you run into problems involving the `<filesystem>` header, you probably need to update your compiler. The compiler `g++-8` should be sufficient. Install a new compiler and direct CMake to use it with the `-DCMAKE_CXX_COMPILER` option. For example, to set the compiler to `g++-8`, run CMake with the option `-DCMAKE_CXX_COMPILER=g++-8`.

Note: the Elias-Fano variants make use of the `_pext_u64` instruction in the BMI2 instruction set. Older CPUs might not support this instruction. In that case, we fall back to a simple software implementation, which will ruin the performance of the Elias-Fano variants (those whose variant name starts with "mef").

## Index construction

Below is the command to build the SBWT for input data `example_data/coli3.fna` provided in this repository, with k = 30. The index is written to the file `index.sbwt`.

```
./build/bin/sbwt build -i example_data/coli3.fna -o index.sbwt -k 30
```

This builds the default variant, which is the plain matrix SBWT. Other variants can be specified with the `--variant` option.
The list of all command line options and parameters is below:

```
Construct an SBWT variant.
Usage:
  build [OPTION...]

  -i, --in-file arg             The input sequences as a FASTA or FASTQ 
                                file, possibly gzipped. If the file 
                                extension is .txt, the file is interpreted 
                                as a list of input files, one file on each 
                                line. All input files must be in the same 
                                format.
  -o, --out-file arg            Output file for the constructed index.
  -k, --kmer-length arg         The k-mer length.
  -p, --precalc-length arg      Precalculate SBWT intervals of strings of 
                                this length. Speeds up query, but takes 
                                4^(p+2) bytes of memory. (default: 8)
      --variant arg             The SBWT variant to build. Available 
                                variants:
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
                                blocked9-split (default: plain-matrix)
      --add-reverse-complements
                                Also add the reverse complement of every 
                                k-mer to the index. Warning: this creates a 
                                temporary reverse-complemented duplicate of 
                                each input file before construction. Make 
                                sure that the directory at --temp-dir can 
                                handle this amount of data. If the input is 
                                gzipped, the duplicate will also be 
                                compressed, which might take a while.
      --no-streaming-support    Save space by not building the streaming 
                                query support bit vector. This leads to 
                                slower queries if the query reads are positive (consequtive k-mers are present in the SBWT).
  -t, --n-threads arg           Number of parallel threads. (default: 1)
  -a, --min-abundance arg       Discard all k-mers occurring fewer than 
                                this many times. By default we keep all 
                                k-mers. Note that we consider a k-mer 
                                distinct from its reverse complement. 
                                (default: 1)
  -b, --max-abundance arg       Discard all k-mers occurring more than this 
                                many times. (default: 1000000000)
  -m, --ram-gigas arg           RAM budget in gigabytes (not strictly 
                                enforced). Must be at least 2. (default: 2)
  -d, --temp-dir arg            Location for temporary files. (default: .)
  -v, --verbose                 Print more verbose output.
  -h, --help                    Print usage

Usage example: build -i example_data/coli3.fna -o index.sbwt -k 30
```



To build, e.g., the corrections-sets variant from the already constructed plain-matrix variant (Ecoli_31.sbwt) run:


``` diff
./build/bin/sbwt build-variant \
    -i /data/Ecoli_31.sbwt \
    -o /data/correction-sets.sbwt \
    --variant correction-sets
```

The input file is one of the supplied plain-matrix SBWT indexes corresponding to a dataset used in the experimental evaluation.

The large input indexes required to reproduce the experiments are distributed separately. See the instructions in the ALENEX Artifact Evaluation README.md for information about obtaining the experimental data.

# Running queries

To query for existence of all k-mers in an index for all sequences in a fastq-file, run the following command:

```
./build/bin/sbwt search -i ./correction-sets.sbwt -q example_data/queries.fastq -o out.txt
```

This prints for each query of length n in the input a line containing n-k+1 space-separated integers, which are the ranks of the columns representing the k-mer in the index. If the k-mer is not found, -1 is printed. If the index was built with streaming support (which is the default), the faster streaming query algorithm is automatically used. The full options are:

```
Query all k-mers of all input reads.
Usage:
  search [OPTION...]

  -o, --out-file arg    Output filename.
  -i, --index-file arg  Index input file.
  -q, --query-file arg  The query in FASTA or FASTQ format, possibly
			gzipped. Multi-line FASTQ is not supported. If the
			file extension is .txt, this is interpreted as a
			list of query files, one per line. In this case,
			--out-file is also interpreted as a list of output
			files in the same manner, one line for each input
			file.
  -z, --gzip-output     Writes output in gzipped form. This can shrink the
			output files by an order of magnitude.
  -h, --help            Print usage
```



# Acknowledgements

The command-line parsing and the gzip support are implemented using [cxxopts](https://github.com/jarro2783/cxxopts) and [zstr](https://github.com/mateidavid/zstr) respectively.
