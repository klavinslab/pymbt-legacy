#!/bin/bash
# Runs a simple server for the docs - not stable / production worthy
source venv/bin/activate
python bin/build_docs.py
cd docs
make html && cd _build/html && nohup python -m SimpleHTTPServer 3089
