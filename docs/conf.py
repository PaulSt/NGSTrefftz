from os.path import dirname, isdir, join
import os
import re
import subprocess

from setuptools_scm import get_version
_cmake_args = []


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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'NGSTrefftz'
copyright = '2022, Paul Stocker'
author = 'Paul Stocker'

# The full version, including alpha/beta/rc tags
release = get_version(root='..', relative_to=__file__).split('+')[0]

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [ "myst_nb" ]
extensions = ["sphinx.ext.autodoc","sphinx.ext.mathjax","sphinx.ext.todo","sphinx.ext.githubpages",
              "IPython.sphinxext.ipython_console_highlighting", "IPython.sphinxext.ipython_directive",
              # "jupyter_sphinx.execute",
              "jupyter_sphinx",
              "nbsphinx",
              # "m2r2",
              "myst_parser",
              "sphinxemoji.sphinxemoji",
              ]
# source_suffix = ['.rst', '.md']

# Run notebook configuration

# The template used when exporting from nbconvert
#   full  - Outputs the full HTML document [Default]
#   basic - Outputs a single div (with no additional resources)
run_notebook_export_template = 'basic'  # Default: 'full'

# Display the source links to the generated evaluated files
run_notebook_display_source_links = False  # Default: True

# Whether or not to evaluate the notebooks prior to embedding them
evaluate_notebooks = True  # Default: True

# START nbsphinx stuff
#increase timeout for cell execution, since some files take long to execute
nbsphinx_timeout = 100000

# If True, the build process is continued even if an exception occurs:
nbsphinx_allow_errors = False

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
.. raw:: html

    <style>
        .p-Widget {
            height: 400px;
        }
        .dg.main {
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li {
            list-style: none;
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li div.dg {
            margin-bottom: 0px;
        }

        .MathJax { font-size: 0.9em !important; }
    </style>

.. only:: html
    .. role:: raw-html(raw)
        :format: html
"""

nbsphinx_widgets_path = ""
# END nbsphinx stuff

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_static']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'paper', 'env', 'jupyter_execute']

autoclass_content = 'both'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_context = {
  'display_github': True,
  'github_user': 'PaulSt',
  'github_repo': 'NGSTrefftz',
  # 'github_version': 'main/docs/',
}
# html_theme_options = {
    # 'github_user': 'PaulSt',
    # 'github_repo': 'NGSTrefftz',
    # # 'github_banner':True,
    # # 'travis_button':True,
    # 'github_button':True,
    # 'fixed_sidebar':True,
    # 'font_size':10
    # # 'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    # # 'analytics_anonymize_ip': False,
# }
# html_sidebars = {
   # 'index': [
        # 'about.html',
        # 'localtoc.html',
        # # 'relations.html',
       # ],
   # '**': ['about.html',
          # 'localtoc.html',
       # ],
# }


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [ 'css/mytheme.css' ]
html_style = 'css/mytheme.css'


# html_js_files = ['webgui_jupyter_widgets.js', 'webgui.js', 'tentswebgui.js']
