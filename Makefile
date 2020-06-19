image=rnacentral/rnacentral-import-pipeline
tag=latest
sif=$(tag).sif
docker=$(image):$(tag)

docker: Dockerfile requirements.txt .dockerignore
	docker build --build-arg CACHE_DATE="$(date)" -t "$(docker)" .

shell:
	docker run -v `pwd`:/rna/import-pipeline -i -t "$(docker)" bash

publish: docker
	docker push "$(docker)"

.PHONY: docker publish clean
