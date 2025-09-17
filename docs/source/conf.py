# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# Mock imports
from unittest.mock import MagicMock

class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()

MOCK_MODULES = ['numpy', 'cython', 'pysam','shutil','stat', 'setuptools', 'sc3dg.bwa', 'sc3dg.bowtie2', 'sc3dg.minimap2', 'sc3dg.nanoplexer', 'sc3dg.bedtools']
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)

project = 'STARK'
copyright = '2024, Wulab'
author = 'Wulab'
release = '1.0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx_click.ext',
    'recommonmark',
    'nbsphinx',
    'sphinx_rtd_theme'
]

templates_path = ['_templates']
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
pygments_style = "sphinx"
master_doc = "index"
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
