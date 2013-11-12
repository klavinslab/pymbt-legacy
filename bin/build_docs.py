#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Uses nbconvert to convert all ipynbs in the ipython_examples dir until I
can figure out a way to use a proper sphinx extension."""

# Doesn"t even use IPython API (TODO!)
# TODO: catch conversion errors (right now they pass silently)
import os
import subprocess


def ipynb_to_rst(directory, filename):
    """Converts ipynb to rst."""
    os.chdir(directory)
    subprocess.Popen(["ipython", "nbconvert", "--to", "rst",
                      filename],
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE)

# The ipython_examples dir has to be in the same dir as this script
script_dir = os.path.dirname(os.path.abspath(__file__))
docs_dir = os.path.join(os.path.dirname(script_dir), "docs")
dirs_files = []
for root, subfolders, files in os.walk(docs_dir):
    for f in files:
        if f.endswith("ipynb"):
            ipynb_to_rst(root, f)
