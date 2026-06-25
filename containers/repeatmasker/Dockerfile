# RepeatMasker with the HMMER (nhmmer) search engine against Dfam.
#
# The Dfam library (FamDB .h5) is NOT baked in - it is large and versioned
# independently, so it is bind-mounted at runtime over $RM/Libraries/famdb
# (see config/repeatmasker.config container_options). The minimal curated Dfam
# bundled in the RepeatMasker tarball is kept only so `configure` and the build
# smoke test succeed.

FROM debian:trixie-slim AS builder
ARG REPEATMASKER_VERSION

ENV RNA=/rna
WORKDIR $RNA

RUN apt update && apt install -y \
    curl \
    gcc \
    make \
    perl \
    python3 \
    python3-pip \
    python3-h5py \
    hmmer \
    trf \
    && rm -rf /var/lib/apt/lists/*

# RepeatMasker is a Perl/Python tool - "building" is unpacking + configure.
RUN curl -OL https://www.repeatmasker.org/RepeatMasker/RepeatMasker-${REPEATMASKER_VERSION}.tar.gz && \
    tar -xzf RepeatMasker-${REPEATMASKER_VERSION}.tar.gz -C $RNA && \
    rm RepeatMasker-${REPEATMASKER_VERSION}.tar.gz && \
    cd $RNA/RepeatMasker && \
    perl ./configure \
      -default_search_engine hmmer \
      -hmmer_dir /usr/bin \
      -trf_prgm /usr/bin/trf \
      -libdir $RNA/RepeatMasker/Libraries

FROM debian:trixie-slim
RUN apt update && apt install -y \
    perl \
    python3 \
    python3-h5py \
    hmmer \
    trf \
    && rm -rf /var/lib/apt/lists/*

ENV RNA=/rna
WORKDIR $RNA

COPY --from=builder $RNA/RepeatMasker $RNA/RepeatMasker
ENV PATH="$RNA/RepeatMasker:$PATH"

# Smoke test: tool resolves (mirrors the infernal container's `cmscan -h`).
# A full -engine hmmer run is exercised in the pipeline, not at build time, to
# avoid depending on a complete Dfam being present in the image.
RUN RepeatMasker -v

ENTRYPOINT ["/bin/bash"]
