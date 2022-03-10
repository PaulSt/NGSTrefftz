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
release = 'v0.0.4'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [ "myst_nb" ]
extensions = ["sphinx.ext.autodoc","sphinx.ext.mathjax","sphinx.ext.todo","sphinx.ext.githubpages",
              "IPython.sphinxext.ipython_console_highlighting", "IPython.sphinxext.ipython_directive",
              "jupyter_sphinx.execute",
              "RunNotebook",
              "nbsphinx",
              "m2r2"
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
evaluate_notebooks = False  # Default: True

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
    </style>

.. only:: html
    .. role:: raw-html(raw)
        :format: html

"""

# END nbsphinx stuff

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'paper','env']

# LaTeX
latex_engine = 'xelatex'
latex_elements = {
    'fontpkg': r'''
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}
''',
    'preamble': r'''
\usepackage[titles]{tocloft}
\cftsetpnumwidth {1.25cm}\cftsetrmarg{1.5cm}
\setlength{\cftchapnumwidth}{0.75cm}
\setlength{\cftsecindent}{\cftchapnumwidth}
\setlength{\cftsecnumwidth}{1.25cm}
''',
    'fncychap': r'\usepackage[Bjornstrup]{fncychap}',
    'printindex': r'\footnotesize\raggedright\printindex',
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'
html_theme_options = {
    'github_user': 'PaulSt',
    'github_repo': 'NGSTrefftz',
    # 'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    # 'analytics_anonymize_ip': False,
}
html_sidebars = {
   'index': ['localtoc.html'],
   '**': [],
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_js_files = ['webgui_jupyter_widgets.js', 'webgui.js', 'tentswebgui.js']
