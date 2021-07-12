# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import subprocess
import os
import sys
from sphinx.builders.html import StandaloneHTMLBuilder
import sphinx_rtd_theme
sys.path.insert(0, os.path.abspath('./rocs/python/'))
# sys.path.insert(0, os.path.abspath('./python/'))
# sys.path.insert(0, os.path.abspath('./src/'))


# -- Project information -----------------------------------------------------

project = 'ROCS'
copyright = '2021, Yinan Li'
author = 'Yinan Li'

# The full version, including alpha/beta/rc tags
release = '2.0'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# subprocess.call('pip install sphinx-prompt', shell=True)
# subprocess.call('pip install breathe', shell=True)
extensions = [
    'sphinx-prompt',
    'breathe',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    # 'sphinx.ext.imgmath',
    'sphinx.ext.todo',
    'sphinx.ext.napoleon',
    'sphinx.ext.graphviz'
]

# Breathe configurations
subprocess.call('make clean', shell=True)
subprocess.call('cd ./rocs/doc/ ; doxygen', shell=True)
breathe_projects = { "ROCS": "./rocs/doc/xml/" }
# subprocess.call('cd ./doxygen ; doxygen', shell=True)
# breathe_projects = { "ROCS": "./doxygen/xml/" }
breathe_default_project = "ROCS"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


StandaloneHTMLBuilder.supported_image_types = [
    'image/svg+xml',
    'image/gif',
    'image/png',
    'image/jpeg'
]
