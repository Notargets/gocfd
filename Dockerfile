FROM alpine

RUN apk update
RUN apk add --no-cache git make musl-dev go
RUN apk add libx11-dev
RUN apk add libxi-dev
RUN apk add libxinerama-dev
RUN apk add libxrandr-dev
RUN apk add libxcursor-dev

# Configure Go
ENV GOROOT /usr/lib/go
ENV GOPATH /go
ENV PATH /go/bin:$PATH

RUN mkdir -p ${GOPATH}/src ${GOPATH}/bin

WORKDIR $GOPATH

CMD ["make test"]
