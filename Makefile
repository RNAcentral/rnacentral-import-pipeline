image=rnacentral/rnacentral-import-pipeline
tag=latest
sif=$(tag).sif
docker=$(image):$(tag)

rust:
	cargo build --release
	mv -f target/release/json2fasta bin
	mv -f target/release/split-ena bin
	mv -f target/release/expand-urs bin
	mv -f target/release/precompute bin
	mv -f target/release/search-export bin
	mv -f target/release/ftp-export bin
	mv -f target/release/json2dfasta bin
	mv -f target/release/bed-expander bin

clean:
	rm bin/json2fasta
	rm bin/split-ena
	rm bin/expand-urs
	rm bin/precompute
	rm bin/search-export
	rm bin/ftp-export
	rm bin/json2dfasta
	rm bin/bed-expander
	cargo clean

docker: Dockerfile
	docker build -t "$(docker)" .

shell: docker
	docker run -v `pwd`:/rna/import-pipeline -i -t "$(docker)"

publish: docker
	docker push "$(docker)"

.PHONY: docker publish clean env rust
