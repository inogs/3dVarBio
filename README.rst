Basic Instructions
=====================

In order to run on PICO machine @Cineca you need to load some modules with the following commands:

.. code::

   module load profile/advanced
   module load autoload netcdf/4.1.3--intel--cs-xe-2015--binary

The only version of VTune that works with this modules is vtune/15.1 (because of the intel compiler automatically loaded by netcdf...).

NetCDF 4.3.2
===============

Added **x86_64.LINUX.netcdf4.3.2.inc** file.

After load the proper modules with the following command:

.. code::

   module load autoload netcdff/4.4.2--intel--cs-xe-2015--binary

with:

.. code::

   cp x86_64.LINUX.netcdf4.3.2.inc compiler.inc
   make


you will compile with the latest version of the netcdf library available on PICO.
