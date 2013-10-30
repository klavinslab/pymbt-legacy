#!/bin/bash
# Runs a simple server for the docs - not stable / production worthy
source venv/bin/activate
python bin/convert_ipynb_to_rst.py
cd docs
make html && cd _build/html && nohup python -m SimpleHTTPServer 3089 &
