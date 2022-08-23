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
	mv target/release/json2fasta bin
	mv target/release/split-ena bin
	mv target/release/expand-urs bin
	mv target/release/precompute bin
	mv target/release/search-export bin
	mv target/release/ftp-export bin
	mv target/release/json2dfasta bin
	mv target/release/expression-parse bin

docker: Dockerfile requirements.txt .dockerignore
	docker build -t "$(docker)" .

shell: docker
	docker run -v `pwd`:/rna/import-pipeline -i -t "$(docker)"

publish: docker
	docker push "$(docker)"

.PHONY: docker publish clean env rust
