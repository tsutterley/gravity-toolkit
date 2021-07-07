FROM python:3.9-slim-buster

LABEL Tyler Sutterley "tsutterl@uw.edu"

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/Los_Angeles"

RUN useradd --create-home --shell /bin/bash gravity

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        libgeos-dev \
        libhdf5-dev \
        libnetcdf-dev \
        libproj-dev \
        libxml2-dev \
        libxslt1-dev \
        openssl \
        proj-data \
        proj-bin && \
    apt-get clean

WORKDIR /home/gravity

RUN pip3 install --no-cache-dir --no-binary=cartopy \
        cartopy \
        future \
        h5py \
        lxml \
        matplotlib \
        netCDF4 \
        numpy \
        python-dateutil \
        pyyaml \
        scipy && \
    pip3 install --no-cache-dir --no-deps git+https://github.com/tsutterley/read-GRACE-geocenter.git && \
    pip3 install --no-cache-dir --no-deps git+https://github.com/tsutterley/geoid-toolkit.git

COPY . .

RUN --mount=source=.git,target=.git,type=bind \
    pip install --no-cache-dir --no-deps .

USER gravity

CMD ["bash"]
