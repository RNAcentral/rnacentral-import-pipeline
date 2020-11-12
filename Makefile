image=rnacentral/rnacentral-import-pipeline
tag=latest
sif=$(tag).sif
docker=$(image):$(tag)

env: requirements.txt requirements-dev.txt

requirements.txt: requirements.in
	pip-compile -o requirements.txt requirements.in

requirements-dev.txt: requirements-dev.in
	pip-compile -o requirements-dev.txt requirements-dev.in

docker: Dockerfile requirements.txt .dockerignore
	docker build -t "$(docker)" .

shell: docker
	docker run -v `pwd`:/rna/import-pipeline -i -t "$(docker)"

publish: docker
	docker push "$(docker)"

.PHONY: docker publish clean env
