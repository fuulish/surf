Build and installation
======================

Build using CMake in a separate directory (out of source build).
That keeps the source tree clean and allows separating multiple builds.
For a basic build, just do:

    $ mkdir build
    $ cd build
    $ cmake .. [options]
    $ make

Use `make VERBOSE=1` to see the actual compilation commands.

Build options
-------------

Option                             | Description
:--------------------------------- | :----------------------------------------
`-DENABLE_OPENMP=ON`               | Enable OpenMP support (on by default)

Installation is currently not supported, just run the binary from the build directory.


Input file
==========

Possible input options:
(have to be given in lower case (I know, I know) then ' = ' and then input value (see also tests directory))

    option          valid keywords
    ------          --------------

    task            surface_distribution - calculates instantaneous surface (string)

    structure       XYZ coordinate input (needed with actual atom names) (string)

    trajectory      either XYZ, or XTC file, if XYZ, each snapshot should be the same number of bytes as every other one

    xdrread         activate XTC reading (integer, 0 - off, everything else - on)

    batch           batch processing mode (three integers: start stop iteration, 0-based)

    pbc             only orthogonal cells, expects three reals (x,y,z), default unit Bohr, if prepended by ANG or NM, angstrom or nano-meters can be given 
    periodic        activate periodicity (integer, 0 - off, everything else - on)

    resolution      expects real, giving grid spacing (only in Bohr!!!)

    refmask         give atoms of one solute of which the distance to the surface will be calculated

    mask            atoms used for surface construction
                        
                        atoms 1 6 8
                        indices 0 1 2 34

                        or similar

    output          string: silent, normal, debug

    surfxyz         turn on output of 2D surface as XYZ file (integer, 0 - off, 1 - two files, surface up and down without atoms, 2 - two files, surface up and down, with atoms)
