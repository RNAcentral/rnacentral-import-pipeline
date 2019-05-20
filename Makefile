image=rnacentral/rnacentral-import-pipeline
tag=latest
sif=$(tag).sif
docker=$(image):$(tag)

singularity/Vagrantfile:
	cd singularity && vagrant init singularityware/singularity-2.4 --box-version 2.4

singularity/ssh-config: singularity/Vagrantfile
	cd singularity && vagrant up && vagrant ssh-config > ssh-config

docker: Dockerfile requirements.txt .dockerignore 
	docker build -t "$(docker)" .
	docker push "$(docker)"

singularity/$(sif): singularity/ssh-config Dockerfile requirements.txt
	ssh -F singularity/ssh-config default "sudo singularity build $(sif) docker://$(docker)"
	scp -F singularity/ssh-config "default:$(sif)" $@

.PHONY: docker
