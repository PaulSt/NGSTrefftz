.. NGSTrefftz documentation master file, created by
   sphinx-quickstart on Wed Mar  2 15:24:26 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

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



Welcome 
===========================
NGSTrefftz is an add-On to NGSolve for Trefftz methods.
Take a look at the notebooks! 
You can run them online `here`_. 

.. _here: https://mybinder.org/v2/gh/PaulSt/NGSTrefftz/HEAD?filepath=doc%2Fnotebooks%2Findex.ipynb


.. toctree::
   :caption: Introduction
   :maxdepth: 1

   readme
   contrib

.. toctree::
   :caption: Documentation

   docu



.. toctree::
   :caption: Notebooks
   :maxdepth: 1

   notebooks/laplace.ipynb
   notebooks/helmholtz.ipynb
   notebooks/twave.ipynb
   notebooks/qtwave.ipynb
   notebooks/twavetents.ipynb
   notebooks/embTrefftz.ipynb
   notebooks/embTrefftz-poi.ipynb
   notebooks/embTrefftz-wave.ipynb
   notebooks/embTrefftz-helm.ipynb
   notebooks/embTrefftz-adv.ipynb

