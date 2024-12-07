gocfd:
	go fmt ./...  && go install ./...
	@printf "run this -> %s \n" "$$"\GOPATH/bin/$@

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
