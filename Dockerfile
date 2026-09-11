FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        tzdata \
        g++ \
        gcc \
        cmake \
        git \
        python3-dev \
        g++-10 \
        libz-dev \
        libbz2-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /SBWT

COPY . .

RUN mkdir -p build && \
    cd build && \
    cmake .. \
        -DCMAKE_CXX_COMPILER=g++-10 \
        -DMAX_KMER_LENGTH=32 \
        -DBUILD_TESTS=1

RUN cmake --build build -j8

CMD ["/SBWT/build/bin/sbwt"]
