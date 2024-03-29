try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    cython = True
    from Cython.Build import cythonize
except ImportError:
    cython = False

config = {
    'description': 'pymbt',
    'author': 'Nick Bolten',
    'url': 'https://github.com/klavinslab/pymbt',
    'download_url': 'https://github.com/klavinslab/pymbt.git',
    'author_email': 'nbolten _at_ gmail',
    'version': '0.1',
    'install_requires': ['nose', 'numpy', 'biopython'],
    'extras_require': {'plotting': ['matplotlib'],
                       'documentation': ['sphinx'],
                       'alignment': ['cython'],
                       'yeast_database': ['intermine', 'requests']},
    'packages': ['pymbt',
                 'pymbt.analysis',
                 'pymbt.analysis._sequence',
                 'pymbt.analysis._sequencing',
                 'pymbt.analysis._structure',
                 'pymbt.constants',
                 'pymbt.database',
                 'pymbt.design',
                 'pymbt.design._oligo_synthesis',
                 'pymbt.design._sequence_generation',
                 'pymbt.seqio',
                 'pymbt.reaction',
                 'pymbt.sequence'],
    'package_data': {'pymbt': ['pymbt/analysis/_sequencing/data/*']},
    'include_package_data': True,
    'scripts': [],
    'name': 'pymbt',
    'license': 'Apache 2.0'
}

if cython:
    import numpy
    setup(ext_modules=cythonize(['pymbt/analysis/_sequencing/calign.pyx']),
          test_suite='nose.collector',
          include_dirs=[numpy.get_include()],
          **config)
else:
    setup(test_suite='nose.collector',
          **config)
