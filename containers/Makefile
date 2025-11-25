# Container image names
cpat=rnacentral/cpat:latest

# Tool versions
INFERNAL_VERSION=1.1.2
SAMTOOLS_VERSION=1.18

# Docker repository
DOCKER_REPO=rnacentral

# Tool container image names
INFERNAL_IMAGE=$(DOCKER_REPO)/infernal:$(INFERNAL_VERSION)
INFERNAL_LATEST=$(DOCKER_REPO)/infernal:latest
SAMTOOLS_IMAGE=$(DOCKER_REPO)/samtools:$(SAMTOOLS_VERSION)
SAMTOOLS_LATEST=$(DOCKER_REPO)/samtools:latest

.PHONY: all cpat infernal samtools push-infernal push-samtools push-all clean test

all: infernal samtools

# Build CPAT container (existing target)
cpat: cpat/Dockerfile
	docker build -t "$(cpat)" -f $^ .
	docker push $(cpat)

# Build Infernal container
infernal:
	@echo "Building Infernal $(INFERNAL_VERSION) container..."
	docker build -t $(INFERNAL_IMAGE) -f infernal/Dockerfile .
	docker tag $(INFERNAL_IMAGE) $(INFERNAL_LATEST)

# Build Samtools container
samtools:
	@echo "Building Samtools $(SAMTOOLS_VERSION) container..."
	docker build -t $(SAMTOOLS_IMAGE) -f samtools/Dockerfile .
	docker tag $(SAMTOOLS_IMAGE) $(SAMTOOLS_LATEST)

# Push Infernal to Docker Hub
push-infernal: infernal
	@echo "Pushing Infernal images to Docker Hub..."
	docker push $(INFERNAL_IMAGE)
	docker push $(INFERNAL_LATEST)

# Push Samtools to Docker Hub
push-samtools: samtools
	@echo "Pushing Samtools images to Docker Hub..."
	docker push $(SAMTOOLS_IMAGE)
	docker push $(SAMTOOLS_LATEST)

# Push all tool containers
push-all: push-infernal push-samtools

# Test tool containers locally
test: infernal samtools
	@echo "Testing Infernal container..."
	docker run --rm $(INFERNAL_IMAGE) cmscan -h > /dev/null
	docker run --rm $(INFERNAL_IMAGE) cmbuild -h > /dev/null
	@echo "✓ Infernal tests passed"
	@echo ""
	@echo "Testing Samtools container..."
	docker run --rm $(SAMTOOLS_IMAGE) samtools --version
	docker run --rm $(SAMTOOLS_IMAGE) tabix --version
	@echo "✓ Samtools tests passed"

# Clean up local tool container images
clean:
	@echo "Removing tool container images..."
	docker rmi $(INFERNAL_IMAGE) $(INFERNAL_LATEST) || true
	docker rmi $(SAMTOOLS_IMAGE) $(SAMTOOLS_LATEST) || true
