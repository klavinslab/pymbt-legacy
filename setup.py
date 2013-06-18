try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'pymbt',
    'author': 'Nick Bolten',
    'url': 'https://github.com/klavinslab/pymbt',
    'download_url': 'https://github.com/klavinslab/pymbt.git',
    'author_email': 'nbolten _at_ gmail',
    'version': '0.1',
    'install_requires': ['nose'],
    'requires': ['numpy',
                 'biopython',
                 'matplotlib'],
    'packages': ['pymbt',
                 'pymbt.analysis',
                 'pymbt.design',
                 'pymbt.io',
                 'pymbt.reaction',
                 'pymbt.sequence'],
    'scripts': [],
    'name': 'pymbt',
    'license': 'GPLv3'
}

setup(test_suite='nose.collector',
      **config)
