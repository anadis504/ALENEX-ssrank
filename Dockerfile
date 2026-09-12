FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc-10 \
    g++-10 \
    cmake \
    git \
    python3-dev \
    zlib1g-dev \
    libbz2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /SBWT

COPY . .

RUN mkdir -p build && \
    cd build && \
    cmake .. \
        -DCMAKE_CXX_COMPILER=g++-10 \
        -DMAX_KMER_LENGTH=32

RUN cmake --build build -j8

RUN ln -s /SBWT/build/bin/sbwt /usr/local/bin/sbwt

CMD ["sbwt"]
