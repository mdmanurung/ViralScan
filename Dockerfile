# Release build requires the tested wheel as a separate, explicit BuildKit context:
# docker build --build-context wheel=./dist -t viralscan:3.0.0rc1 .
# The main context contains only this file, .dockerignore and the committed lock.
# 26.9.4-0 multi-arch index (Ubuntu 24.04, built 2026-10-24). The previous pin was a 2024-06 Ubuntu 20.04
# image; grype reports that distro as end of life. Refresh: docker buildx imagetools inspect condaforge/miniforge3:<tag>.
FROM condaforge/miniforge3@sha256:63d72eadf0ab73ee63bad8f5cfdfe1754530415b4aeba99e466e9a6f77ceec95

LABEL maintainer="emma.vonk@hotmail.nl" \
      version="3.0.0rc1" \
      description="ViralScan — viral load quantification from single-cell RNA-seq"

# Security updates for the base layers (PLAN REL-11): the Ubuntu libbrotli1 deb and the
# Miniforge base-env Python, neither of which the locked viralscan env replaces.
RUN apt-get update \
    && apt-get install -y --only-upgrade libbrotli1 \
    && rm -rf /var/lib/apt/lists/*

COPY environment.explicit.txt /tmp/environment.explicit.txt
RUN conda install --yes --name base 'python>=3.14.10,<3.15' \
    && mamba create --yes --name viralscan --file /tmp/environment.explicit.txt \
    && mamba clean --all --force-pkgs-dirs --yes

SHELL ["conda", "run", "-n", "viralscan", "/bin/bash", "-c"]
RUN --mount=type=bind,from=wheel,source=.,target=/tmp/dist,ro \
    python -m pip install --no-deps /tmp/dist/*.whl \
    && python -m pip check \
    && viralscan doctor --profile full

WORKDIR /data
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "viralscan", "viralscan"]
CMD ["--help"]
