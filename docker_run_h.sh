#!/bin/bash
sudo docker run -v `pwd`:/mnt med2rdf/hint -c /mnt/tsv2rdf_hint.json
