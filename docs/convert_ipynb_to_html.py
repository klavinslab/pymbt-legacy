# -*- coding: utf-8 -*-
#!/usr/bin/env python
'''Uses nbconvert to convert all ipynbs in the ipython_examples dir until I
can figure out a way to use a proper sphinx extension.'''

# Doesn't even use IPython API (TODO!)
import os
import subprocess

# The ipython_examples dir has to be in the same dir as this script
basedir = os.path.dirname(os.path.abspath(__file__))
ipynb_dirname = "examples"
ipynb_dir = os.path.join(basedir, ipynb_dirname)
os.chdir(ipynb_dir)
ipynb_files = filter(lambda x: x.endswith(".ipynb"), os.listdir(ipynb_dir))
for ipynb in ipynb_files:
    converter = subprocess.Popen(['ipython', 'nbconvert', '--to', 'rst',
                                  ipynb],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
