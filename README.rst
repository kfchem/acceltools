acceltools
==========

**acceltools** is an extension toolkit for `ACCeL <https://github.com/kfchem/accel>`_,
designed to enhance molecular structure analysis and computational chemistry workflows.

It supports document export, ECD/UV/NMR analysis, reaction plotting, puckering
parameter evaluation, and structural scanning, all based on molecular data
managed by ACCeL.

Features
--------

- Export molecular data to TXT / DOCX formats
- Calculate and visualize ECD and UV spectra
- Assign and compare NMR chemical shifts
- Generate reaction energy diagrams and scatter plots
- Analyze Cremer–Pople puckering parameters
- Perform bond/angle/dihedral coordinate scanning
- Output structured results as pandas.DataFrame

Installation
------------

.. code-block:: bash

    pip install acceltools

Dependencies
------------

- ACCeL (https://github.com/kfchem/accel)
- Python >= 3.9
- numpy, pandas, matplotlib, scipy, python-docx

License
-------

MIT License

Repository
----------

GitHub: https://github.com/kfchem/acceltools
