#!/bin/bash

docker build -t masters-gene-expression . && \
  docker tag masters-gene-expression ctlab.registry.jetbrains.space/bioinfdocker/masters-gene-expression:latest && \
  docker push ctlab.registry.jetbrains.space/bioinfdocker/masters-gene-expression:latest