Installation
============


System requirements
-------------------

- Operating system:

  - Mac OSX (tested on >=10.11), or
  - Ubuntu (tested on 20.04 LTS), or
  - Manjaro (tested on 5.9.16).

- Compilers: gcc (clang), g++ (clang++)



Prerequisite packages and softwares
-----------------------------------

.. note::
   To simplify the editing of the configuration file *cflags*, we suggest install **Boost**, **Armadillo** and **HDF5** by using *Homebrew* (for Mac OSX) or *apt* (for Ubuntu).


- `Boost c++ libraries <http://www.boost.org>`_ >= 1.69.0: libraries for graph searching, ode solvers;
- `Armadillo <http://arma.sourceforge.net>`_ (tested on 9.900.2, 9.850.1, and 10.2.2): a c++ linear algebra library;
- `HDF5 <https://www.hdfgroup.org/downloads/hdf5/>`_ (1.10.0/1.12.0): a high-performance data management and storage suite. Install the package from source code may take a while.

  + For Mac OSX users, using `Homebrew <https://formulae.brew.sh/formula/hdf5>`_ for installation is stongly recommended.
  + For Ubuntu users, it is also suggested to install this `package in Focal <https://launchpad.net/ubuntu/focal/+source/hdf5>`_ by running

    .. prompt:: bash $

       sudo apt install libhdf5-dev

  + Manjaro users can also use this `package <https://archlinux.org/packages/community/x86_64/hdf5>`_

- MATLAB: control simulation and graphics. **Installation of MATLAB is no longer mandatory for saving control synthesis results. Data saved in HDF5 format can also be loaded in MATLAB for simulation.** Update the environment variable `DYLD_LIBRARY_PATH` in terminal if you have MATLAB installed and wish to save data in MATLAB data format:

  ``export DYLD_LIBRARY_PATH=“Your_MATLAB_install_path/bin/maci64:$DYLD_LIBRARY_PATH"``

  ``export DYLD_LIBRARY_PATH=“Your_MATLAB_install_path/sys/os/maci64:$DYLD_LIBRARY_PATH”``

 - Python 3: an alternative of MATLAB. Make sure the packages **scipy**, **numpy**, **matplotlib**, and **h5py** are installed. The following command can be used:

   .. prompt:: bash $

      pip3 install scipy numpy matplotlib h5py


Configuration of makefiles
--------------------------
In the configuration file *cflags* in the root path, the paths assigned to `INCS` and `LDFLAGS` might need to be changed accordingly.

- For Mac OSX users, modify the path for the external packages to your actual install path, e.g.,

  ``EXTDIR := Your_package_install_path``

  If the prerequisites are installed in different paths, make sure the paths are given correctly to `INCS` and `LDFLAGS`.

- For Linux users, rename the file *cflags-linux* to *cflags*.

  * On Ubuntu, if HDF5 is installed via `apt` and `mpi` is not available, HDF5 may be installed only for serial compuation. In this case, the following changes may apply to the *cflags*:

    + Add to `INCS` the following path:

      ``-I$(INCDIR)/hdf5/serial``

    + Change `-lhdf5_hl` and `-lhdf5` to `-lhdf5_serial_hl` and `-lhdf5_serial`, respectively.

  * If HDF5 is installed by compiling the source code from HDF5 official website, make sure that *zlib* and *szip* are installed as in their instruction (*release_docs/INSTALL*). Use

    .. prompt:: bash $

       ./configure --prefix=Your_install_path --enable-cxx --with-szip=Path_to_szip

  to configure *cmake*. After installation, add to `INCS` the following paths:

  ``-I$(Your_install_path)/include -I$(Path_to_szip)/include``

  and add to `LDFLAGS` the following paths:

  ``-I$(Your_install_path)/lib -I$(Path_to_szip)/lib``

- For all users, if you choose to save control synthesis results in MATLAB data format, make sure that MATLAB is installed and edit:

  ``MATDIR := Your_MATLAB_install_path``

  Otherwise, comment out or delete the line starting with `MATDIR`.
