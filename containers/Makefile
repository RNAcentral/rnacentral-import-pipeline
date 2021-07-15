cpat=rnacentral/cpat:latest

cpat: cpat/Dockerfile
	docker build -t "$(cpat)" -f $^ .
	docker push $(cpat)

.PHONY: cpat
