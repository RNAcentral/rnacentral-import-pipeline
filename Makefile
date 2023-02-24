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
	mv -f target/release/json2fasta bin
	mv -f target/release/split-ena bin
	mv -f target/release/expand-urs bin
	mv -f target/release/precompute bin
	mv -f target/release/search-export bin
	mv -f target/release/ftp-export bin
	mv -f target/release/json2dfasta bin
	mv -f target/release/expression-parse bin
	mv -f target/release/bed-expander bin

clean:
	rm bin/json2fasta
	rm bin/split-ena
	rm bin/expand-urs
	rm bin/precompute
	rm bin/search-export
	rm bin/ftp-export
	rm bin/json2dfasta
	rm bin/expression-parse
	cargo clean

docker: Dockerfile requirements.txt .dockerignore
	docker build -t "$(docker)" .

shell: docker
	docker run -v `pwd`:/rna/import-pipeline -i -t "$(docker)"

publish: docker
	docker push "$(docker)"

.PHONY: docker publish clean env rust
