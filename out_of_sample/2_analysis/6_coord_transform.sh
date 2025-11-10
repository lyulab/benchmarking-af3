#!/bin/bash

for d in finished_outputs/*/; do
	python align_structures.py "$d"
done
