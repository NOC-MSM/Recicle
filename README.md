[![DOI](https://zenodo.org/badge/153254595.svg)](https://zenodo.org/badge/latestdoi/153254595)

# RECICLE

## Configuration

The following code was used in the Recicle project:

```
# Checkout configuration directory structure
git clone git@github.com:NOC-MSM/Recicle.git
cd Recicle
```
At this point you can checkout and compile XIOS or use a version you already have. If you're starting from scratch:

```
# Choose an appropriate directory for your XIOS installation
svn co http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0@703
cd xios-1.0
./make_xios --full --prod --arch $ARCH --netcdf_lib netcdf4_par --job 4
```

Where `$ARCH` is the architecture file (see xios documentation for details). Next, compile the NEMO code itself. First we copy the arch files into the appropriate directory.

```
cd $WORK_DIR/NEMO-shelf/NEMOGCM/CONFIG/
```

NB you will either have to write or modify an architecture file in the `./ARCH` directory in order to comile NEMO in your environment.

```
./makenemo -n AMM7_RECICLE -m $ARCH -j 4
```

That should be enough to produce a valid executable. Now to copy the forcing data from JASMIN. 

```
cd AMM7_RECICLE/ENSEMBLE_INPUTS
wget -r -np -nH --cut-dirs=3 -erobots=off --reject="index.html*" http://gws-access.ceda.ac.uk/public/recicle/config/
```

And finally link the XIOS binary to the configuration directory.

```
cd ../ENSEMBLE_CONTROL
rm xios_server.exe
ln -s $XIOS_DIR/xios-1.0/bin/xios_server.exe xios_server.exe
```

Edit and run the ```ensemble.pbs``` script in ```../ENSEMBLE_CONTROL``` accordingly.
