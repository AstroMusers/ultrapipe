ULTRApipe
=========
**ULTRApipe** is a Python library for tools and pipelines created for the ULTRASAT mission. This repository will be actively updated as more resources become available. API (via readthedocs) will be available soon.

Currently, we offer tools for:

- Conversions from TESS magntiude to ULTRASAT magnitude

.. image:: /media/graphics/TESS_to_ULTRASAT_Mag.png
  :width: 400
  :alt: Histogram of the TESS/ULTRASAT magnitude of known or PC exoplanetary systems.

Histogram of the TESS/ULTRASAT magnitude of known or PC exoplanetary systems.

- Dust Maps

.. image:: /media/graphics/Dust_Map.png
  :width: 600
  :alt: 3D dust map using the Pan-STARRS 1, 2MASS, and Gaia mission generated using data by `Green et al. 2019 <http://argonaut.rc.fas.harvard.edu/>`. The above maps show the magntiude extinction in the UV band (~230-300 nm) for declinations above -30 degrees at varying distances.

3D dust map using the Pan-STARRS 1, 2MASS, and Gaia mission generated using data by `Green et al. 2019 <http://argonaut.rc.fas.harvard.edu/>`. The above maps show the magntiude extinction in the UV band (~230-300 nm) for declinations above -30 degrees at varying distances.

- CVZ Cost Functions for Exoplanet Systems

.. image:: /media/graphics/Cost_Function.png
  :width: 600
  :alt: A cost function map of known or PC TOI systems for different proposed 7x7 degree continous viewing zones. The weights and parameterization of the cost map can be found in the `utils.py` file.

A cost function map of known or PC TOI systems for different proposed 7x7 degree continous viewing zones. The weights and parameterization of the cost map can be found in the `utils.py` file.
- ULTRASAT Systematics and Throughput Data
Several useful data files are also hosted here, detailing the throughput of ULTRASAT over its passband, as well as the limiting and saturation magnitudes of different radial positions on the detector for different colored targets.
