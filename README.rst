.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/casskit.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/casskit
    .. image:: https://readthedocs.org/projects/casskit/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://casskit.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/casskit/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/casskit
    .. image:: https://img.shields.io/pypi/v/casskit.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/casskit/
    .. image:: https://img.shields.io/conda/vn/conda-forge/casskit.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/casskit
    .. image:: https://pepy.tech/badge/casskit/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/casskit
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/casskit

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

=======
casskit
=======


    Add a short description here!


A longer description of your project goes here...


.. _pyscaffold-notes:

=========
Summary of related packages
=========


.. list-table:: Fetching and parsing TCGA
   :widths: 30 10 10 65
   :header-rows: 1

   * - Name
     - Language
     - Active?
     - Description
   * - `TCGAbiolinks`_
     - R
     - Yes
     - An R/Bioconductor package for integrative analysis with TCGA data
   * - `TCGAutils`_
     - R
     - Yes
     - ``MultiAssayExperiment`` for TCGA


.. list-table:: Biological data structures
   :widths: 30 10 10 65
   :header-rows: 1

   * - Name
     - Language
     - Active?
     - Description
   * - ExperimentHub
     - R
     - Yes
     - A repository of curated biological data
   * - pyGeno
     - Python
     - Yes
     - precision medicine and proteogenomics
   * - `MultiAssayExperiment`_
     - R
     - Yes
     - A Bioconductor package for the representation of multi-assay experiments
   * - `scverse`_
     - R, Python
     - Yes
     - ``AnnData``, ``muon`` with ``PyTorch`` for single-cell RNA-seq
   * - `sgkit`_
     - Python
     - Yes
     - ``xarray`` for VCFs with some statgen functionality (eg GWAS)
   * - `scikit-allel`_
     - Python
     - No
     - Succeeded by ``sgkit``
   * - `scikit-genome`_
     - Python
     - Yes
     - add-on to CNVkit
   * - `scikit-bio`_
     - Python
     - Yes
     - mainly microbial genomics
   * - `dalmatian`_
     - Python
     - Yes
     - a collection of high-level functions for interacting with Firecloud via Pandas dataframes
   * - `HTSeq 2.0`_
     - Python
     - Yes
     - HTSeq is a Python package for analysis of high-throughput sequencing data


.. list-table:: Annotations
   :widths: 30 10 10 65
   :header-rows: 1

   * - Name
     - Language
     - Active?
     - Description
   * - `pypath`_ / `OmniPath`_
     - Python, R
     - Yes
     - A Python module for molecular signaling prior knowledge processing
   * - `pyensembl`_
     - Python
     - Yes
     - annotation
   * - eDGAR
     - Python
     - Yes
     - a database of Disease-Gene Associations


.. list-table:: Other
   :widths: 30 10 10 65
   :header-rows: 1

   * - Name
     - Language
     - Active?
     - Description
   * - `PyBDA`_
     - Python
     - Yes
     - A Python package for the analysis of biological data
   * - PyBEL
     - Python
     - Yes
     - A Python module for biological expression language
   * - pycellbase
     - Python
     - Yes
     - mainly microbial genomics
   * - pygenometracks
     - Python
     - Yes
     - 


=========
Development roadmap
=========

see here


Note
====

This project has been set up using PyScaffold 4.3.1. For details and usage
information on PyScaffold see https://pyscaffold.org/.


.. Refs
.. =====
.. _MultiAssayExperiment: https://github.com/waldronlab/MultiAssayExperiment
.. _OmniPath: https://omnipathdb.org
.. _PyBDA: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3087-8
.. _pycellbase: https://pypi.org/project/pycellbase/
.. _pyensembl: https://raw.githubusercontent.com/openvax/pyensembl/0e750e50105c22666fcd43181183719876e15e6a/README.md
.. _pypath: https://github.com/saezlab/pypath
.. _scikit-allel: https://scikit-allel.readthedocs.io/en/stable/
.. _scikit-bio: http://scikit-bio.org
.. _scikit-genome: https://cnvkit.readthedocs.io/en/stable/skgenome.html
.. _scverse: https://scverse.org
.. _sgkit: https://pystatgen.github.io/sgkit/latest/
.. _TCGAutils: https://github.com/waldronlab/TCGAutils
.. _TCGAbiolinks: https://github.com/BioinformaticsFMRP/TCGAbiolinks
.. _dalmatian: https://github.com/getzlab/dalmatian
.. _HTSeq 2.0: https://htseq.readthedocs.io/en/master/index.html