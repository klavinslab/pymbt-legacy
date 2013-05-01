try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages

config = {
    'description': 'pymbt',
    'author': 'Nick Bolten',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'nbolten _at_ gmail',
    'version': '0.1',
    'install_requires': ['nose'],
    'requires': ['numpy',
                 'biopython',
                 'matplotlib'],
    'packages': ['pymbt',
                 'pymbt.oligo_synthesis',
                 'pymbt.sequence_analysis',
                 'pymbt.sequence_generation'],
    'scripts': [],
    'name': 'pymbt',
    'license': 'GPLv3'
}

setup(**config)
