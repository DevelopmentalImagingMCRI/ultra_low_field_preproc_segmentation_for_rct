FROM ubuntu:18.04
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update

RUN apt-get install -y curl build-essential cmake-curses-gui python3-dev git virtualenv libpng-dev wget libopenblas-dev liblapack-dev gfortran

RUN wget -O- http://neuro.debian.net/lists/bionic.us-ca.full | tee /etc/apt/sources.list.d/neurodebian.sources.list && apt-key adv --recv-keys --keyserver hkps://keyserver.ubuntu.com 0xA5D32F012649A5A9

RUN apt-get update && apt-get install -y fsl-core

RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null

RUN echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ bionic main' | tee /etc/apt/sources.list.d/kitware.list >/dev/null

RUN apt-get update && apt-get install -y cmake

RUN mkdir -p /code/
WORKDIR /code/
RUN mkdir c3d && curl -L https://sourceforge.net/projects/c3d/files/c3d/1.0.0/c3d-1.0.0-Linux-x86_64.tar.gz | tar xz -C c3d --strip-components=1

ARG REPO=ITK_NiftyMIC
ARG COMMIT=26415ac2e6197de7b07ffcb0c3f740aa937ba7e9
RUN mkdir "${REPO}"  && curl -L "https://codeload.github.com/gift-surg/"${REPO}"/tar.gz/${COMMIT}" | tar -xz -C "${REPO}" --strip-components=1

ARG REPO=NSol
ARG COMMIT=d47b4840b37035281000f35768aa8b0c960ea583
RUN mkdir "${REPO}"  && curl -L "https://codeload.github.com/richardbeare/"${REPO}"/tar.gz/${COMMIT}" | tar -xz -C "${REPO}" --strip-components=1

ARG REPO=NiftyMIC
ARG COMMIT=fb341053a6504d08042fbba48dd204f4cb9af284
RUN mkdir "${REPO}"  && curl -L "https://codeload.github.com/richardbeare/"${REPO}"/tar.gz/${COMMIT}" | tar -xz -C "${REPO}" --strip-components=1
 
ARG REPO=PySiTK
ARG COMMIT=6bb9419b9e5ec7c20aab443172b86da69a268b9a
RUN mkdir "${REPO}"  && curl -L "https://codeload.github.com/gift-surg/"${REPO}"/tar.gz/${COMMIT}" | tar -xz -C "${REPO}" --strip-components=1

ARG REPO=niftyreg
ARG COMMIT=6db8b16c17884a9b8859c980c4f1c408f62bd9ca
RUN mkdir "${REPO}"  && curl -L "https://codeload.github.com/KCL-BMEIS/"${REPO}"/tar.gz/${COMMIT}" | tar -xz -C "${REPO}" --strip-components=1

RUN virtualenv --python=python3.6 /code/venv
RUN (. /code/venv/bin/activate ; pip install SimpleITK==1.2.4)


ENV NIFTYREG_INSTALL=/code/niftyreg_install
ENV NIFTYMIC_ITK_DIR=/code/ITK_NiftMIC-build

RUN . /code/venv/bin/activate && mkdir /code/ITK_NiftyMIC_build && cd /code/ITK_NiftyMIC_build &&\
  cmake   -D CMAKE_BUILD_TYPE=Release   -D BUILD_TESTING=OFF   -D BUILD_EXAMPLES=OFF \
  -D BUILD_SHARED_LIBS=ON   -D ITK_WRAP_PYTHON=ON   -D ITK_LEGACY_SILENT=ON \
  -D ITK_WRAP_float=ON   -D ITK_WRAP_double=ON   -D ITK_WRAP_signed_char=ON \
  -D ITK_WRAP_signed_long=ON   -D ITK_WRAP_signed_short=ON \
  -D ITK_WRAP_unsigned_char=ON   -D ITK_WRAP_unsigned_long=ON \
  -D ITK_WRAP_unsigned_short=ON   -D ITK_WRAP_vector_float=ON \
  -D ITK_WRAP_vector_double=ON   -D ITK_WRAP_covariant_vector_double=ON \
  -D Module_ITKReview=ON   -D Module_SmoothingRecursiveYvvGaussianFilter=ON \
  -D Module_BridgeNumPy=ON   ../ITK_NiftyMIC/ && \
  make -j8 && \
  cp Wrapping/Generators/Python/WrapITK.pth /code/venv/lib/python3.6/site-packages/

RUN .  /code/venv/bin/activate && mkdir /code/niftyreg_build && mkdir /code/niftyreg_install && \
    cd  /code/niftyreg_build && \
    cmake -D CMAKE_INSTALL_PREFIX=/code/niftyreg_install ../niftyreg && \
    make -j4 && \
    make install 

RUN cd /code/PySiTK/ && (. /code/venv/bin/activate; pip install .)

RUN cd /code/NSoL && (. /code/venv/bin/activate ; pip install -e .)

RUN cd /code/NiftyMIC && (. /code/venv/bin/activate ; pip install -e .)
  
COPY entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]

