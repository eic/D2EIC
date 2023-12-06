# Work in progress

This is a github repository for the Detector 2 at the EIC (D2EIC).
This conceptual design of the detector is intended for the IP8.

The software infrastructure is based on the existing epic software, which shares a common
ancestry. 

This page is being developed. 
The contact persons now are Bill Li and Kong Tu

Dec 5 2023.


Overview
--------

**Detector geometry viewer:**
- [Viewer only](https://eic.github.io/D2EIC/geoviewer)


Getting Started
---------------

Get a copy of the latest version from this repository:
```bash
git clone https://github.com/eic/D2EIC.git
```

### Compilation

To configure, build, and install the geometry (to the `install` directory), use the following commands:
```bash
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build
cmake --install build
```
To load the geometry, you can use the scripts in the `install` directory:
```bash
source install/setup.sh
```
or
```tcsh
source install/setup.csh
```

### Adding/changing detector geometry

Hint: **Use the CI/CD pipelines**.

To avoid dealing with setting up all the dependencies, we recommend using the continuous integration/continuous deployment (CI/CD) pipelines to make changes and assess their effects. Any feedback to help this process is appreciated.

Here is how to begin:

1. Look at existing detector constructions and reuse if possible. Note that "compact detector descriptions" -> xml files, and "detector construction" -> cpp file.
2. Modify xml file or detector construction.
3. Create a WIP (or draft) merge request or pull request and look at the CI output for debugging. Then go to back to 2 if changes are needed.
4. Remove the WIP/Draft part of the merge request if you would like to see your changes merged into the main.

See:

- [Talk at computing round table](https://indico.jlab.org/event/420/#17-automated-workflow-for-end)

### Compiling (avoid it)

First, see if the use case above is best for you. It most likely is and can save a lot of time for newcomers.
To run the simulation locally, we suggest using the singularity image.
More details can be found at the links below:

- https://dd4hep.web.cern.ch/dd4hep/page/beginners-guide/
- https://eic.phy.anl.gov/tutorials/eic_tutorial/
- https://eicweb.phy.anl.gov/containers/eic_container/


Related useful links
--------------------

- [EIC tutorial](https://eic.phy.anl.gov/tutorials/eic_tutorial)
- [DD4hep repository](https://github.com/AIDAsoft/DD4hep)
- [DD4hep user manual](https://dd4hep.web.cern.ch/dd4hep/usermanuals/DD4hepManual/DD4hepManual.pdf)
- [ACTS DD4hep plugin documentation](https://acts.readthedocs.io/en/latest/plugins/dd4hep.html)
