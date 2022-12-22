=========
TO-DO
=========

- [ ] Flexible passing of keyword arguments to underlying ``dask`` functions
- [ ] Extend test suite with ``hypothesis`` [`hypothesis docs`_] and ``pytest-benchmark`` [`pytest-benchmark docs`_]
- [ ] Identifying mislabeled samples with `COSMO`_ (split out from nf implementation). Include known TCGA mislabeling (Are known samples eg `Subchallenge 1 Training Key`_ corrected?)
- [ ] TCGA CPTAC data
- [ ] Add native support for `SMOTE` and other oversampling methods
- [ ] Support custom Data structure `serializations in dask`
- [ ] Add tutorial using `PyDESeq2`
- [ ] *cis*-mapping tutorial
- [ ] Add stratified sampling function using `imblearn.FunctionSampler`


.. Refs
.. =====
.. _Subchallenge 1 Training Key: https://precision.fda.gov/challenges/4
.. _COSMO : https://github.com/bzhanglab/COSMO
.. _hypothesis docs: https://hypothesis.readthedocs.io/en/latest/index.html
.. _imblearn.FunctionSampler: https://imbalanced-learn.org/stable/references/generated/imblearn.FunctionSampler.html
.. _PyDESeq2: https://github.com/FedeGerva/pydeseq2
.. _pytest-benchmark docs: https://pytest-benchmark.readthedocs.io/en/latest/
.. _serializations in dask: https://distributed.dask.org/en/stable/serialization.html#dask-serialization-family
.. _SMOTE: https://github.com/analyticalmindsltd/smote_variants