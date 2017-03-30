# PolyglycerolModeler
A Python tool constructing 3D models of polyglycerol polymers with a user-defined degree of branches ranging from linear and hyperbranched to dendritic. The PDB coordinate file comes with a topology/parameter file according to the AMBER force field and ready for simulations with Gromacs.

Prerequisites
-------------

#### Required libraries

On a Ubuntu 16.04 Linux system, two additional libraries, `python-pygraphviz` and `python-biopython` are required in order to be able to run the Python script. On such a Debian-based system one would execute

`sudo apt install python-pygraphviz python-biopython`

to install both packages and all further dependencies.


#### Modify Gromacs topology files for usage with particular Amber force field

If you want to use the generated PDB coordinates file as input for an AMBER parameterization step using the `gmx pdb2gmx` command, you will need to modify a couple of files of that particular AMBER force field in Gromacs topology directory (usually in `/usr/share/gromacs/top`). You might want to copy the Gromacs topology directory to an own user directory before

`cp -a /usr/share/gromacs/top ~/my-gromacs-top`

and set the corresponding environment variable for Gromacs commands

`export GMXLIB=/home/$USER/my-gromacs-top`

Assuming, we want to use the Amber99sb force field, the following files would be extended by the content of the respective files in the "gmx_topology" directory of PolyglycerolModeler:

`cat ./PolyglycerolModeler/gmx_topology/residuetypes.dat >> ~/my-gromacs-top/residuetypes.dat`<br />
`cat ./PolyglycerolModeler/gmx_topology/specbond.dat >> ~/my-gromacs-top/specbond.dat`<br />
`cat ./PolyglycerolModeler/gmx_topology/aminoacids.hdb >> ~/my-gromacs-top/amber99sb.ff/aminoacids.hdb`<br />
`cat ./PolyglycerolModeler/gmx_topology/aminoacids.rtp >> ~/my-gromacs-top/amber99sb.ff/aminoacids.rtp`<br />

In addition, increase the number in the first line of `~/my-gromacs-top/specbond.dat` by 32 (which is the number of newly added special bonds).




Usage
-----

coming soon

