#! /usr/bin/env bash
docker buildx build \
  -t bguo068/hmmibd-rs:v0.1.7 \
  -t bguo068/hmmibd-rs:latest \
  --platform linux/amd64,linux/arm64 \
  --push .
