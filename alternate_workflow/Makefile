.PHONY: help build build-nocache run resume validate dry-run test unit-test unit-test-fast diagnose-h5ad shell clean logs

STEP ?= cell_type_annotation

IMAGE_NAME := t_cell_workflow
IMAGE_TAG := latest
CONFIG ?= config_standard.json
CONFIG_DIR ?= $(PWD)
DATA_DIR ?= $(PWD)
OUTPUT_DIR ?= $(PWD)
LOGS_DIR ?= $(OUTPUT_DIR)/logs
DATA_DIR := $(shell echo $(DATA_DIR))
OUTPUT_DIR := $(shell echo $(OUTPUT_DIR))
LOGS_DIR := $(shell echo $(LOGS_DIR))
CPUS ?= 4
MEM ?= 90g
DOCKER_RM ?= --rm

ifeq ($(filter /%,$(CONFIG)),)
  _CONFIG_HOST := $(CONFIG_DIR)/$(CONFIG)
else
  _CONFIG_HOST := $(CONFIG)
endif
_CONFIG_NAME := $(notdir $(_CONFIG_HOST))

DOCKER_RUN := docker run $(DOCKER_RM) \
	--platform linux/amd64 \
	$(if $(MEM),--memory=$(MEM) --memory-swap=-1) --cpus=$(CPUS) --shm-size=8g \
	--cap-add=SYS_PTRACE \
	--entrypoint bash \
	-e SF_LOG_DIR=/logs \
	-e HOST_UID=$(shell id -u) \
	-e HOST_GID=$(shell id -g) \
	-e HDF5_USE_FILE_LOCKING=FALSE \
	-v $(DATA_DIR):/raw \
	-v $(OUTPUT_DIR):/results \
	-v $(LOGS_DIR):/logs \
	-v $(_CONFIG_HOST):/app/$(_CONFIG_NAME):ro

help:
	@echo ""
	@echo "T cell Workflow - Docker Operations"
	@echo "==========================================="
	@echo "make build Build Docker image"
	@echo "make build-nocache Rebuild image without cache"
	@echo ""
	@echo "Pipeline Commands"
	@echo "-----------------"
	@echo "make run Run full workflow (default: 4 CPUs)"
	@echo "make validate Validate configuration"
	@echo "make dry-run Simulate pipeline without execution"
	@echo "make test Alias for dry-run"
	@echo "make unit-test Rebuild the image, then run the R unit tests inside it"
	@echo "make unit-test-fast Run tests against host source (no rebuild; for iteration only)"
	@echo ""
	@echo "Utilities"
	@echo "---------"
	@echo "make shell Open interactive shell in container"
	@echo "make diagnose-h5ad Diagnose .h5ad reading (BPCells + metadata backends)"
	@echo "make logs Show recent docker logs"
	@echo "make clean Remove stopped containers"
	@echo ""
	@echo "Configuration Variables"
	@echo "-----------------------"
	@echo "CPUS=N Set number of CPUs (default: 4)"
	@echo "MEM=Ng Set container memory limit (default: 90g; pass MEM= empty for unlimited)"
	@echo "DOCKER_RM= Disable auto-remove for post-mortem debugging (default: --rm)"
	@echo "CONFIG=file.json Set config file name (default: config_standard.json)"
	@echo "CONFIG_DIR=/path Directory containing the config file (default: current directory)"
	@echo "DATA_DIR=/path Host directory mounted as /raw (default: current directory)"
	@echo "OUTPUT_DIR=/path Host directory mounted as /results (default: current directory)"
	@echo "LOGS_DIR=/path Host directory mounted as /logs (default: OUTPUT_DIR/logs)"
	@echo ""
	@echo "Examples:"
	@echo " make run CONFIG=my_config.json"
	@echo " make run CONFIG=my_config.json CONFIG_DIR=/mnt/configs"
	@echo " make run CONFIG=/mnt/citeseq/data/my_config.json DATA_DIR=/mnt/citeseq/data"
	@echo " make run CONFIG=/mnt/citeseq/data/my_config.json DATA_DIR=/mnt/citeseq/data OUTPUT_DIR=/mnt/citeseq/results"
	@echo " make run CONFIG=/mnt/citeseq/data/my_config.json DATA_DIR=/mnt/citeseq/data OUTPUT_DIR=/mnt/results LOGS_DIR=/mnt/logs"
	@echo " make run CONFIG=/mnt/citeseq/data/my_config.json DATA_DIR=/mnt/citeseq/data OUTPUT_DIR=/mnt/citeseq CPUS=8"
	@echo ""

define ensure_mounts
	@mkdir -p "$(OUTPUT_DIR)" "$(LOGS_DIR)" "$(DATA_DIR)"
	@test -f "$(_CONFIG_HOST)" || { \
	  echo "ERROR: config file not found: $(_CONFIG_HOST)"; \
	  echo " pass CONFIG=<file> and optionally CONFIG_DIR=<dir>, or CONFIG=/abs/path.json"; \
	  exit 1; }
endef

build:
	@echo "Building Docker image: $(IMAGE_NAME):$(IMAGE_TAG)"
	docker build --platform linux/amd64 -t $(IMAGE_NAME):$(IMAGE_TAG) .
	@echo "[OK] Build complete!"

build-nocache:
	@echo "Rebuilding image without cache..."
	docker build --platform linux/amd64 --no-cache -t $(IMAGE_NAME):$(IMAGE_TAG) .

run:
	@echo "Running workflow with config: $(_CONFIG_HOST)"
	$(ensure_mounts)
	$(DOCKER_RUN) \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	/app/entrypoint.sh run $(_CONFIG_NAME)

resume:
	@echo "Resuming workflow from step: $(STEP)"
	$(ensure_mounts)
	$(DOCKER_RUN) \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	/app/entrypoint.sh resume $(_CONFIG_NAME) $(STEP)

validate:
	@echo "Validating configuration: $(_CONFIG_HOST)"
	$(ensure_mounts)
	$(DOCKER_RUN) \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	/app/entrypoint.sh validate $(_CONFIG_NAME)

dry-run:
	@echo "Performing pipeline dry-run..."
	$(ensure_mounts)
	$(DOCKER_RUN) \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	/app/entrypoint.sh dry-run $(_CONFIG_NAME)

test: dry-run

unit-test: build
	@echo "Running unit tests inside $(IMAGE_NAME):$(IMAGE_TAG)..."
	docker run --rm \
	--platform linux/amd64 \
	--entrypoint bash \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	-c "cd /app && Rscript tests/testthat.R"

unit-test-fast:
	@echo "Running unit tests against HOST source (no rebuild - not what ships)..."
	docker run --rm \
	--platform linux/amd64 \
	--entrypoint bash \
	-v "$(PWD)/tests:/app/tests:ro" \
	-v "$(PWD)/script:/app/script:ro" \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	-c "cd /app && Rscript tests/testthat.R"

diagnose-h5ad:
	@echo "Diagnosing h5ad handling..."
	docker run --rm \
	--platform linux/amd64 \
	--entrypoint bash \
	-v "$(PWD)/..:/raw:ro" \
	-v "$(PWD)/diagnose_h5ad.R:/app/diagnose_h5ad.R:ro" \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	-c "cd /app && Rscript diagnose_h5ad.R"

shell:
	@echo "Opening container shell (image /app, not host)..."
	@mkdir -p "$(OUTPUT_DIR)" "$(LOGS_DIR)" "$(DATA_DIR)"
	docker run --rm -it \
	--platform linux/amd64 \
	--entrypoint bash \
	-e HOST_UID=$(shell id -u) \
	-e HOST_GID=$(shell id -g) \
	-v $(DATA_DIR):/raw \
	-v $(OUTPUT_DIR):/results \
	-v $(LOGS_DIR):/logs \
	$(IMAGE_NAME):$(IMAGE_TAG) \
	/app/entrypoint.sh shell

logs:
	@echo "Showing recent container logs..."
	docker logs $$(docker ps -lq)

clean:
	@echo "Cleaning Docker containers..."
	docker container prune -f
	@echo "[OK] Cleanup complete!"
