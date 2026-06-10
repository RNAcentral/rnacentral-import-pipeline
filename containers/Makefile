# Tool versions - read from VERSION files
INFERNAL_VERSION := $(shell cat infernal/VERSION)
SAMTOOLS_VERSION := $(shell cat samtools/VERSION)
CPAT_VERSION := $(shell cat cpat/VERSION)
EMBOSS_VERSION := $(shell cat emboss/VERSION)

# Rust utilities version - auto mode uses git SHA
RUST_VERSION_RAW := $(shell cat rust-utils/VERSION)
ifeq ($(RUST_VERSION_RAW),auto)
    RUST_GIT_SHA := $(shell git rev-parse --short=8 HEAD 2>/dev/null || echo "unknown")
    ifneq ($(shell git diff-index --quiet HEAD -- 2>/dev/null; echo $$?),0)
        RUST_VERSION := $(RUST_GIT_SHA)-dirty
    else
        RUST_VERSION := $(RUST_GIT_SHA)
    endif
else
    RUST_VERSION := $(RUST_VERSION_RAW)
endif

# Docker repository
DOCKER_REPO=rnacentral

# Tool container image names
CPAT_LATEST=$(DOCKER_REPO)/cpat:latest
CPAT_IMAGE=$(DOCKER_REPO)/cpat:$(CPAT_VERSION)
EMBOSS_IMAGE=$(DOCKER_REPO)/emboss:$(EMBOSS_VERSION)
EMBOSS_LATEST=$(DOCKER_REPO)/emboss:latest
INFERNAL_IMAGE=$(DOCKER_REPO)/infernal:$(INFERNAL_VERSION)
INFERNAL_LATEST=$(DOCKER_REPO)/infernal:latest
SAMTOOLS_IMAGE=$(DOCKER_REPO)/samtools:$(SAMTOOLS_VERSION)
SAMTOOLS_LATEST=$(DOCKER_REPO)/samtools:latest
RUST_UTILS_IMAGE=$(DOCKER_REPO)/rust-utils:$(RUST_VERSION)
RUST_UTILS_LATEST=$(DOCKER_REPO)/rust-utils:latest

.PHONY: all cpat emboss infernal samtools rust-utils push-cpat push-emboss push-infernal push-samtools push-rust-utils push-all clean test validate-versions

all: infernal samtools rust-utils emboss

# Build CPAT container (existing target)
cpat: cpat/Dockerfile
	docker buildx build -t $(CPAT_IMAGE) -f $^ --platform linux/amd64 .
	docker tag $(CPAT_IMAGE) $(CPAT_LATEST)

# Push CPAT to Docker Hub
push-cpat: cpat
	@echo "Pushing CPAT images to Docker Hub..."
	docker push $(CPAT_IMAGE)
	docker push $(CPAT_LATEST)

# Build EMBOSS container with procps for Nextflow compatibility
emboss: emboss/Dockerfile
	@echo "Building EMBOSS $(EMBOSS_VERSION) container..."
	docker buildx build -t $(EMBOSS_IMAGE) -f $^ --platform linux/amd64 .
	docker tag $(EMBOSS_IMAGE) $(EMBOSS_LATEST)

# Push EMBOSS to Docker Hub
push-emboss: emboss
	@echo "Pushing EMBOSS images to Docker Hub..."
	docker push $(EMBOSS_IMAGE)
	docker push $(EMBOSS_LATEST)

# Build Infernal container
infernal:
	@echo "Building Infernal $(INFERNAL_VERSION) container..."
	docker buildx build -t $(INFERNAL_IMAGE) -f infernal/Dockerfile --platform linux/amd64  --build-arg INFERNAL_VERSION=$(INFERNAL_VERSION) .
	docker tag $(INFERNAL_IMAGE) $(INFERNAL_LATEST)

# Build Samtools container
samtools:
	@echo "Building Samtools $(SAMTOOLS_VERSION) container..."
	docker buildx build -t $(SAMTOOLS_IMAGE) -f samtools/Dockerfile --platform linux/amd64 --build-arg SAMTOOLS_VERSION=$(SAMTOOLS_VERSION) .
	docker tag $(SAMTOOLS_IMAGE) $(SAMTOOLS_LATEST)

# Build Rust utilities container
rust-utils:
	@echo "Building Rust utilities container (version: $(RUST_VERSION))..."
	docker buildx build \
		--build-arg RUST_VERSION=$(RUST_VERSION) \
		-t $(RUST_UTILS_IMAGE) \
		-f rust-utils/Dockerfile \
		--platform linux/amd64 ../
	docker tag $(RUST_UTILS_IMAGE) $(RUST_UTILS_LATEST)
	@echo "✓ Rust utilities $(RUST_VERSION) built successfully"

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
	docker push $(RUST_UTILS_LATEST)

# Push all tool containers
push-all: push-cpat push-emboss push-infernal push-samtools push-rust-utils

# Test tool containers locally
test: cpat emboss infernal samtools rust-utils
	@echo "Testing CPAT container..."
	docker run --rm $(CPAT_IMAGE) cpat -h
	@echo "✓ CPAT tests passed"
	@echo ""
	@echo "Testing EMBOSS container..."
	docker run --rm $(EMBOSS_IMAGE) tcode -help
	docker run --rm $(EMBOSS_IMAGE) ps --version
	@echo "✓ EMBOSS tests passed"
	@echo ""
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

# Validate VERSION files
	validate-versions:
		@echo "Validating VERSION files..."
		@for tool in infernal samtools rust-utils cpat emboss; do \
			if [ ! -f $$tool/VERSION ]; then \
				echo "ERROR: Missing VERSION file: $$tool/VERSION"; \
			exit 1; \
		fi; \
		echo "✓ $$tool/VERSION: $$(cat $$tool/VERSION)"; \
	done

# Clean up local tool container images
	clean:
		@echo "Removing tool container images..."
		docker rmi $(CPAT_IMAGE) $(CPAT_LATEST) || true
		docker rmi $(EMBOSS_IMAGE) $(EMBOSS_LATEST) || true
		docker rmi $(INFERNAL_IMAGE) $(INFERNAL_LATEST) || true
		docker rmi $(SAMTOOLS_IMAGE) $(SAMTOOLS_LATEST) || true
		docker rmi $(RUST_UTILS_IMAGE) $(RUST_UTILS_LATEST) || true
