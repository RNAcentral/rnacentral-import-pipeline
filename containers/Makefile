# Container image names
cpat=rnacentral/cpat:latest

# Tool versions
INFERNAL_VERSION=1.1.2
SAMTOOLS_VERSION=1.22.1
CPAT_VERSION=3.0.5

# Docker repository
DOCKER_REPO=rnacentral

# Tool container image names
INFERNAL_IMAGE=$(DOCKER_REPO)/infernal:$(INFERNAL_VERSION)
INFERNAL_LATEST=$(DOCKER_REPO)/infernal:latest
SAMTOOLS_IMAGE=$(DOCKER_REPO)/samtools:$(SAMTOOLS_VERSION)
SAMTOOLS_LATEST=$(DOCKER_REPO)/samtools:latest
RUST_UTILS_IMAGE=$(DOCKER_REPO)/rust-utils:latest

.PHONY: all cpat infernal samtools rust-utils push-infernal push-samtools push-rust-utils push-all clean test

all: infernal samtools rust-utils

# Build CPAT container (existing target)
cpat: cpat/Dockerfile
	docker buildx build -t "$(cpat)" -f $^ --platform linux/amd64 .
	docker push $(cpat)

# Build Infernal container
infernal:
	@echo "Building Infernal $(INFERNAL_VERSION) container..."
	docker buildx build -t $(INFERNAL_IMAGE) -f infernal/Dockerfile --platform linux/amd64 .
	docker tag $(INFERNAL_IMAGE) $(INFERNAL_LATEST)

# Build Samtools container
samtools:
	@echo "Building Samtools $(SAMTOOLS_VERSION) container..."
	docker buildx build -t $(SAMTOOLS_IMAGE) -f samtools/Dockerfile --platform linux/amd64 .
	docker tag $(SAMTOOLS_IMAGE) $(SAMTOOLS_LATEST)

# Build Rust utilities container
rust-utils:
	@echo "Building Rust utilities container..."
	docker buildx build -t $(RUST_UTILS_IMAGE) -f rust-utils/Dockerfile --platform linux/amd64 ../
	@echo "✓ Rust utilities built successfully"

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

# Push Rust utilities to Docker Hub
push-rust-utils: rust-utils
	@echo "Pushing Rust utilities image to Docker Hub..."
	docker push $(RUST_UTILS_IMAGE)

# Push all tool containers
push-all: push-infernal push-samtools push-rust-utils

# Test tool containers locally
test: infernal samtools rust-utils
	@echo "Testing Infernal container..."
	docker run --rm $(INFERNAL_IMAGE) cmscan -h
	docker run --rm $(INFERNAL_IMAGE) cmbuild -h
	@echo "✓ Infernal tests passed"
	@echo ""
	@echo "Testing Samtools container..."
	docker run --rm $(SAMTOOLS_IMAGE) samtools --version
	docker run --rm $(SAMTOOLS_IMAGE) tabix --version
	@echo "✓ Samtools tests passed"
	@echo ""
	@echo "Testing Rust utilities container..."
	docker run --rm $(RUST_UTILS_IMAGE) json2fasta --help
	docker run --rm $(RUST_UTILS_IMAGE) precompute --help
	docker run --rm $(RUST_UTILS_IMAGE) search-export --help
	@echo "✓ Rust utilities tests passed"

# Clean up local tool container images
clean:
	@echo "Removing tool container images..."
	docker rmi $(INFERNAL_IMAGE) $(INFERNAL_LATEST) || true
	docker rmi $(SAMTOOLS_IMAGE) $(SAMTOOLS_LATEST) || true
	docker rmi $(RUST_UTILS_IMAGE) || true
