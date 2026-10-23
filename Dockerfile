# Release build requires the tested wheel as a separate, explicit BuildKit context:
# docker build --build-context wheel=./dist -t viralscan:3.0.0rc1 .
# The main context contains only this file, .dockerignore and the committed lock.
FROM condaforge/miniforge3@sha256:1bfb5a539841983db847561c797f0ee6deca7c1a652b32d6eb457ada5d14069a

LABEL maintainer="emma.vonk@hotmail.nl" \
      version="3.0.0rc1" \
      description="ViralScan — viral load quantification from single-cell RNA-seq"

COPY environment.explicit.txt /tmp/environment.explicit.txt
RUN mamba create --yes --name viralscan --file /tmp/environment.explicit.txt \
    && mamba clean --all --force-pkgs-dirs --yes

SHELL ["conda", "run", "-n", "viralscan", "/bin/bash", "-c"]
RUN --mount=type=bind,from=wheel,source=.,target=/tmp/dist,ro \
    python -m pip install --no-deps /tmp/dist/*.whl \
    && python -m pip check \
    && viralscan doctor --profile full

WORKDIR /data
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "viralscan", "viralscan"]
CMD ["--help"]
