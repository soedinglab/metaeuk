FROM alpine:latest AS metaeuk-builder

RUN apk add --no-cache gcc g++ cmake musl-dev git ninja 

WORKDIR /opt/metaeuk
ADD . .

RUN git submodule init && git submodule update

WORKDIR build_sse
RUN cmake -G Ninja -DHAVE_SSE=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN ninja && ninja install

WORKDIR ../build_avx
RUN cmake -G Ninja -DHAVE_AVX=1 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..
RUN ninja && ninja install

