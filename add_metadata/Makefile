.PHONY: docker run genaiz genaiz-test help

IMAGE_NAME := cite-seq-add-metadata
IMAGE_TAG  := latest

# amd64 pinned: the jammy binary repo only serves x86_64; arm64 source builds fail in this base image.
PLATFORM   := linux/amd64

help:
	@echo "Targets:"
	@echo "  docker      Build Docker image"
	@echo "  run         Run locally with example config"
	@echo "  genaiz      Build GenAIz smart function"
	@echo "  genaiz-test Run GenAIz test"

docker:
	docker build --platform $(PLATFORM) -t $(IMAGE_NAME):$(IMAGE_TAG) .

run:
	Rscript main.R --config configs/config.example.json

genaiz:
	genaiz sf build

genaiz-test: run/in/configs/config.json
	genaiz sf test

run/in/configs/config.json:
	mkdir -p run/in/configs && cp configs/config.example.json run/in/configs/config.json
