FROM bioconductor/bioconductor_docker@sha256:67ab741591c5cd068f8993ca5f0371794c300da9cc5e0b3b2c1c90c9058e4c6f
WORKDIR /work
COPY release/ /work/release/
COPY dist/TALLBMDeconv_0.1.1.tar.gz /work/
COPY examples/PUBLIC_EXAMPLE.R /work/
RUN mkdir -p /work/tmp && TMPDIR=/work/tmp Rscript --vanilla /work/release/restore_locked.R
RUN TMPDIR=/work/tmp Rscript --vanilla /work/release/run_supported_example.R /work/example_check
CMD ["R"]
RUN bash /work/release/check_supported.sh
