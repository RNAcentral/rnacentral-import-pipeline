image=rnacentral/rnacentral-import-pipeline
tag=latest
sif=$(tag).sif
docker=$(image):$(tag)

singularity/$(sif): singularity/ssh-config publish
	ssh -F $< default "sudo singularity build $(sif) docker://$(docker)"
	scp -F $< "default:$(sif)" $@

singularity/Vagrantfile:
	cd singularity && vagrant init singularityware/singularity-2.4 --box-version 2.4

singularity/ssh-config: singularity/Vagrantfile
	cd singularity && vagrant up && vagrant ssh-config > ssh-config

docker: Dockerfile requirements.txt .dockerignore 
	docker build -t "$(docker)" .

shell:
	docker run -i -t "$(docker)" bash

publish: docker
	docker push "$(docker)"

.PHONY: docker publish
