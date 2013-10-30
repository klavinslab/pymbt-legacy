#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Uses nbconvert to convert all ipynbs in the ipython_examples dir until I
can figure out a way to use a proper sphinx extension.'''

# Doesn't even use IPython API (TODO!)
# TODO: catch conversion errors (right now they pass silently)
import os
import subprocess

# The ipython_examples dir has to be in the same dir as this script
script_dir = os.path.dirname(os.path.abspath(__file__))
docs_dir = os.path.join(os.path.dirname(script_dir), "docs")
ipynb_dir = os.path.join(docs_dir, "examples")
os.chdir(ipynb_dir)
ipynb_files = filter(lambda x: x.endswith(".ipynb"), os.listdir(ipynb_dir))
for ipynb in ipynb_files:
    converter = subprocess.Popen(['ipython', 'nbconvert', '--to', 'rst',
                                  ipynb],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
