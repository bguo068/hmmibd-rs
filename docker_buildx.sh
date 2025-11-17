#! /usr/bin/env bash
docker buildx build \
  --platform linux/amd64,linux/arm64 \
  -t bguo068/hmmibd-rs:v0.1.5 \
  -t bguo068/hmmibd-rs:latest \
  --push .
