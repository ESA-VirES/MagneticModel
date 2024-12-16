# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import os
import inspect
import shutil
#import sphinx

try:  # for Sphinx >= 1.7
    from sphinx.ext import apidoc
except ImportError:
    from sphinx import apidoc

__location__ = os.path.join(
    os.getcwd(), os.path.dirname(inspect.getfile(inspect.currentframe()))
)

output_dir = os.path.join(__location__, "api")
module_dir = os.path.join(__location__, "../eoxmagmod/eoxmagmod")


try:
    shutil.rmtree(output_dir)
except FileNotFoundError:
    pass

#try:
#    args = (
#        f"sphinx-apidoc --implicit-namespaces -f -o {output_dir} {module_dir}"
#    ).split(" ")
#    if tuple(sphinx.__version__.split(".")) >= ("1", "7"):
#        # This is a rudimentary parse_version to avoid external dependencies
#        args = args[1:]
#    apidoc.main(args)
#except Exception as e:
#    print("Running `sphinx-apidoc` failed!\n{}".format(e))


def setup(app):
    app.add_css_file("css/custom.css")


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'VirES for Swarm - Magnetic Models, Coordinates and Other Utilities'
copyright = '2022-2024, EOX IT Services GmbH'
author = 'Martin Pačes, Týna Doležalova'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.coverage",
    "sphinx.ext.doctest",
    "sphinx.ext.ifconfig",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    #"myst_parser",
    "sphinx_rtd_theme",
    "sphinx-jsonschema",
    "frigate.sphinx.ext",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "tests"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
#html_static_path = ["static"]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}
