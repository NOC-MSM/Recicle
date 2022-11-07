============================================================================
Fork of NEMO_ 3.6 stable branch with modifcations for FABM_ and shelf systems.
============================================================================

This master branch includes modifications for shelf systems and FABM_ with ERSEM_ enabled.

The pure NEMO_-code without shelf and FABM modifications is in the ``feat/NEMO-svn`` branch.
It should remain untouched in order to enable conflict-free syncing with the central NEMO repository.
The current state is based on svn revision ``@6232``.

Shelf-enabled code without the FABM coupler is in ``feat/shelf-enabled``.

Compilation using ``makenemo``
==============================

To compile on a typical PML workstation using the ``makenemo`` tool provided with NEMO::

  module load mpi #required on fedora
  export XIOS_HOME #set to the basepath of XIOS-1 (compile before, 2 doesn't work!)
  # typical global configuration:
  ./makenemo -m GCC_PMPC -n AMM7

To compile on *archer* using the intel compiler::

  module unload PrgEnv-cray PrgEnv-gnu
  module load PrgEnv-intel
  module load cray-netcdf-hdf5parallel
  export XIOS_HOME #set to the basepath of XIOS-1 (compile before, 2 doesn't work!)
  #typical global:
  ./makenemo -m XC_ARCHER_INTEL_NOSIGNEDZERO -n ORCA2_LIM_FABM

Specific notes on the FABM_-ERSEM_ implementation
=============================================

FABM needs to be compiled separately before the compilation of NEMO with ERSEM enabled.
Usually, the following suffices to achieve this::

   mkdir -p ~/build/nemo && cd ~/build/nemo
   cmake <FABMDIR>/src/ -DFABM_HOST=nemo -DFABM_ERSEM_BASE=<ERSEMDIR> -DFABM_EMBED_VERSION=ON
   make install

In the above, replace `<FABMDIR>` with the directory with the FABM source code, e.g., `~/fabm-git` and `<ERSEMMDIR>` with the directory with the ERSEM_ source code, e.g., `~/ersem-git`.

This will create the library in the standard folder `~/local/fabm/nemo/lib` where NEMO-FABM will look for linking to NEMO.

The FABM coupler for NEMO is added in a sub-folder ``FABM`` in ``NEMOGCM/NEMO/TOP_SRC``.
Changes to existing code in order to accomadate FABM within NEMO are restricted to the ``NEMOGCM/NEMO/TOP_SRC`` and ``NEMOGCM/TOOLS/COMPILE`` folder and shall be marked in the code in the following way:

Additions are encapsulated using the tags::

   ! +++>>> FABM
   ... Some new code...
   ! FABM <<<+++

Removed sections are encapsulated as::

   ! --->>> FABM
   ! ...Some commented code...
   ! FABM <<<---

(In the FCM scripts the ``!`` is replaced by ``#``.)

The initial FABM_ implementation in this repository is carried over from the NEMO-FABM_ repository developed by M. Butenschön and J. Bruggeman, as a single patch commit.

.. _FABM: http://fabm.net
.. _NEMO: http://www.nemo-ocean.eu
.. _ERSEM: http://shelfseasmodelling.org
.. _NEMO-FABM: https://gitlab.ecosystem-modelling.pml.ac.uk/momm/NEMO-FABM

Setting-up XIOS in detached mode using dedicated I/O-servers as mpi tasks
=========================================================================

While it is easiest to run NEMO in attached mode with the I/O operations of XIOS running distributed on the NEMO computation cores, this option is increasingly inefficient with higher degrees of parallelism. In particular, in these cases the maximum number of cores possible is the number of points in j dimension.
It is recomended to set-up the runs in detached mode attributing to XIOS dedicated cores as I/O servers that are not involved in the actual NEMO computations.

To use detached mode the `"using_server"` variable in `iodef.xml` needs to be set to `true`::
  
   <variable id="using_server" type="boolean">true</variable>

(For attached mode the same variable should be set to `false`.)

In addition, it is recommended to set the buffersize for I/O transactions from client to server::

   <variable id="buffer_size" type="integer">10000000</variable>

where the adequate number to use is indicitavely the size of the *local sub-domains* of the NEMO grid times eight (`jpi*jpj*jpk*8`).
Note that the memory buffer required on the server is then given by twice the buffersize times the number of processors the server comunicates with, given by the the number of blocks of the NEMO domain decomposition fitting into a line of XIOS decomposition:
while for NEMO the domain is decomposed in a rectangular i,j structure minimising the surface of the individual blocks and hence the I/O links to neighbouring blocks, for XIOS the domain is decomposed in horizontal lines.

Running in detached mode on archer is achieved by launching::

   aprun -b -n $XIOSCORES -N 1 ./xios_server.exe : -n $NEMOCORES -N 24 ./nemo.exe

where `$XIOSCORES` is the number of I/O-SERVERS and `$NEMOCORES` is the number of compute nodes used for the pure NEMO computations, I/O excluded. `-N` specifies the number of cores used per archer node in the two respective cases.

In addition, the archer architecture consists of nodes with 24 cores on two processors (with 12 nodes each), so if you use more that on server per node it is prudent to specify the distribution on the processors with the `-S` flag giving the number of processes per processor (e.g. running 4 XIOS cores on one node with two servers per processor would require the options `-b -n 4 -N 4 -S 2)`

XIOS-1
==================

Note that NEMO 3.6 stable is incompatible with XIOS-1 versions more recent than September 2015,
due to what is supposed to be a bug-fix, that is incompatible with NEMO 3.6 (XIOS-1 commit of 1st October 2015).

The official NEMO documentation therefore recommends checking out a specific revision (703) of XIOS-1:

http://www.nemo-ocean.eu/Using-NEMO/User-Guides/Basics/XIOS-IO-server-installation-and-use

If you use this official code, you need to add files ``arch/arch-<ARCHITECTURE>.env``, ``arch/arch-<ARCHITECTURE>.fcm``, ``arch/arch-<ARCHITECTURE>.path`` for your computer architecture and OS.
For PML workstations (``<ARCHITECTURE>=GCC_PMPC``), you can base these files on their equivalent for archicture ``GCC_LINUX``;
the only change you need to make is to add ``-DBOOST_DETAIL_NO_CONTAINER_FWD`` to ``BASE_CFLAGS`` in ``arch/arch-GCC_PMPC.fcm``

Note that you can also use the following repository for a customized NEMO 3.6 compatible version of XIOS-1:

https://gitlab.ecosystem-modelling.pml.ac.uk/momm/XIOS1/tree/nemo3.6-fix

This has files for architecture ``GCC_PMPC`` included.

After you obtain the xios code (and optionally, add architecture files), you can compile it on a typical PML workstation with::

   module load mpi #required on fedora
   ./make_xios -arch GCC_PMPC

Troubleshooting
===============

* Missing Perl packages: the fcm compilation system that is used to build xios and nemo depends on several Perl packages including ``URI.pm`` and ``Text/Balanced.pm``. These two packages are not present on all systems. For instance, on the PML Fedora-based workstations they need to be installed through the package manager: ``dnf install install perl-URI``, ``dnf install perl-Text-Balanced``.

* Error building xios: ``.../boost/functional/hash/extensions.hpp:38:33: error: 'template<class T, class A> std::size_t boost::hash_value' conflicts with a previous declaration``. This appears to affect newer versions of GCC. It can be addressed by adding ``-DBOOST_DETAIL_NO_CONTAINER_FWD`` to ``BASE_CFLAGS`` in ``arch/arch-<ARCHITECTURE>.fcm`` (where ``<ARCHITECTURE>`` is the architecture that you provide to ``make_xios`` with ``--arch``.
