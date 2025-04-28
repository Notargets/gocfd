gocfd:
	go fmt ./...
	CGO_LDFLAGS="-lopenblas -llapacke -lgfortran -lm -lpthread" CGO_CFLAGS="-march=native -mavx -mavx2 -mfma" go install -ldflags "-X main.GitVersion=$(git describe --tags --always --dirty)" ./...
	@printf "run this -> %s\n" "$$GOPATH/bin/$@"

test:
	go test -cover ./...

tidy:
	go mod tidy

builder:
	docker build -f Dockerfile_Ubuntu -t gcr.io/gocfd-275017/builder .

push:
	docker push gcr.io/gocfd-275017/builder

generate:
	go generate ./...

bench:
	 go test github.com/notargets/gocfd/model_problems/Euler2D/benchmarks/... -bench=GetFlowFunction
	 #go test github.com/notargets/gocfd/model_problems/Euler2D/benchmarks/... -bench=EulerSolve

libsInstall:
	# Make sure OpenBLAS (with Fortran ABI) is on your system:
	sudo apt update
	sudo apt install libopenblas-dev liblapacke-dev gfortran



	# Install the Netlib BLAS wrapper:
	CGO_LDFLAGS="-lopenblas -lgfortran" go install gonum.org/v1/netlib/blas/netlib
	# Install the Netlib LAPACK wrapper:
	CGO_LDFLAGS="-lopenblas -lgfortran" go install gonum.org/v1/netlib/lapack/netlib

