# syntax=docker/dockerfile:1.6
FROM alpine:3.21

RUN apk update && apk upgrade
RUN apk add --no-cache \
  gcc \
  g++ \
  cmake \
  ninja \
  armadillo-dev \
  onetbb-dev \
  boost-dev \
  bzip2 \
  tar \
  openblas-dev \
  git

WORKDIR /src
ARG DLIB_VERSION=19.24
ADD --checksum=sha256:28fdd1490c4d0bb73bd65dad64782dd55c23ea00647f5654d2227b7d30b784c4 \
  https://dlib.net/files/dlib-${DLIB_VERSION}.tar.bz2 /src/
RUN tar axf dlib-${DLIB_VERSION}.tar.bz2
WORKDIR /src/dlib-${DLIB_VERSION}
RUN cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
RUN cmake --build build
RUN cmake --install build
RUN rm -rf /src/dlib*

RUN mkdir /draco
WORKDIR /draco
COPY CMakeLists.txt .
COPY src ./src
COPY test ./test

ARG BUILD_TYPE=Release
ARG EXTRA_CMAKE_ARGS="-DCMAKE_BUILD_TYPE=Release -DLINK_TIME_OPTIMIZATIONS=ON -DNATIVE_BUILD=ON -DARMA_NO_WRAPPER=ON"

RUN cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=${BUILD_TYPE} ${EXTRA_CMAKE_ARGS}
RUN cmake --build build
RUN cmake --install build

WORKDIR /data
RUN rm -rf /draco

ENTRYPOINT ["draco"]
