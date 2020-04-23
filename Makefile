gocfd:
	go fmt ./...  && go install ./...
	@printf "run this -> %s \n" "$$"\GOPATH/bin/$@

test:
	go test -cover ./...

builder:
	docker build -f Dockerfile_Ubuntu -t gcr.io/gocfd-275017/builder .

push:
	docker push gcr.io/gocfd-275017/builder

