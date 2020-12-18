gocfd:
	go fmt ./...  && go install ./...
	@printf "run this -> %s \n" "$$"\GOPATH/bin/$@

test:
	go test -cover ./...

builder:
	docker build -f Dockerfile_Ubuntu -t gcr.io/gocfd-275017/builder .

push:
	docker push gcr.io/gocfd-275017/builder

generate:
	go generate ./...

bench:
	go test -run BenchmarkEuler_Solve -bench . ./model_problems/Euler2D/ -v
