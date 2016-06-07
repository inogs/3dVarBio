Basic Instructions
=====================

In order to run on PICO machine @Cineca you need to load some modules with the following commands:

.. code::

   module load profile/advanced
   module load autoload netcdf/4.1.3--intel--cs-xe-2015--binary

The only version of VTune that works with this modules is vtune/15.1 (because of the intel compiler automatically loaded by netcdf...).
