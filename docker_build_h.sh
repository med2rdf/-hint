#!/bin/bash
DOCKER_DIR="docker_build_h"

cp -p tsv2rdf_hint.py $DOCKER_DIR
cp -p templ_hint.ttl $DOCKER_DIR
cp -p templ_hint.ttl.prefix $DOCKER_DIR
cp -p templ_hint.ttl.evi $DOCKER_DIR
sudo docker build -t med2rdf/hint $DOCKER_DIR
