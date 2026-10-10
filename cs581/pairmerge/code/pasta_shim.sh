#!/bin/bash
# PASTA's [muscle] path: forwards "-in1 A -in2 B -out OUT -quiet -profile" to pasta_merger.py
exec /usr/bin/python3 "$(dirname "$(readlink -f "$0")")/pasta_merger.py" "$@"
