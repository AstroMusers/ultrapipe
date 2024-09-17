ULTRApipe
=========
**ULTRApipe** is a Python library for the ULTRASAT mission to perform photometry and light curve analysis relevant to exoplanet research. Its documentation and API will be available via readthedocs.

Currently, we offer tools for:

+++++++++++++++++++++++++++++++++++++++++++++++++++++
Color-transformations to/from ULTRASAT passband in NUV (230-290 nm)
+++++++++++++++++++++++++++++++++++++++++++++++++++++

A plot of the TESS and ULTRASAT magnitudes of stars hosting confirmed exoplanets is below.

.. image:: /media/graphics/TESS_to_ULTRASAT_Mag.png
  :width: 400
  :alt: Histogram of the TESS/ULTRASAT magnitude of known or PC exoplanetary systems.

++++++++++
Predicting ULTRASAT magnitudes
++++++++++
.. image:: /media/graphics/Dust_Map.png
  :width: 800
  :alt: 3D dust map using the Pan-STARRS 1, 2MASS, and Gaia mission generated using data by `Green et al. 2019 <http://argonaut.rc.fas.harvard.edu/>`. The above maps show the magntiude extinction in the UV band (~230-300 nm) for declinations above -30 degrees at varying distances.

Computing the ULTRASAT magnitudes taking into account the color transformation and dust extinction in the ULTRASAT band. We use the 3D dust map derived from Pan-STARRS 1, 2MASS, and Gaia `Green et al. 2019 <http://argonaut.rc.fas.harvard.edu/>`_. The above maps show the magntiude extinction in the UV band for declinations above -30 degrees at varying distances.

++++++++++++++++++++++++++++++++++++++++
Estimate of exoplanet science potential for potential ULTRASAT pointings
++++++++++++++++++++++++++++++++++++++++

.. image:: /media/graphics/Cost_Function.png
  :width: 800
  :alt: A cost function map of known or PC TOI systems for different proposed 7x7 degree continous viewing zones. The weights and parameterization of the cost map can be found in the `utils.py` file.

An estimation of the exoplanet science potential of ULTRASAT using the known and candidate exoplanets. The weights and parameterization of the map can be found in the `utils.py` file.

+++++++++++++++++++++++++++++++++++++++++
ULTRASAT Systematics and Throughput
+++++++++++++++++++++++++++++++++++++++++

Several useful data files are also hosted here, detailing the throughput of ULTRASAT over its passband, as well as the limiting and saturation magnitudes of different radial positions on the detector for different colored targets.
