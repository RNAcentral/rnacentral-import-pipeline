image=rnacentral/rnacentral-import-pipeline
tag=latest
sif=$(tag).sif
docker=$(image):$(tag)

env: requirements.txt requirements-dev.txt

requirements.txt: requirements.in
	pip-compile -o requirements.txt requirements.in

requirements-dev.txt: requirements-dev.in
	pip-compile -o requirements-dev.txt requirements-dev.in

rust:
	cargo build --release
	cp target/release/json2fasta bin
	cp target/release/split-ena bin
	cp target/release/kv bin
	cp target/release/expand-urs bin
	cp target/release/precompute bin
	cp target/release/search-export bin
	cp target/release/ftp-export bin

docker: Dockerfile requirements.txt .dockerignore
	docker build -t "$(docker)" .

shell: docker
	docker run -v `pwd`:/rna/import-pipeline -i -t "$(docker)"

publish: docker
	docker push "$(docker)"

.PHONY: docker publish clean env rust
