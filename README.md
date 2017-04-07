# PolymerAssembler

A Python tool constructing 3D models of common polymers yet including various polyglycerol (PG) and soon polyethylen types with a user-defined degree of branches ranging from linear up to hyperbranched and dendritic polymers. The resulting PDB file lacking hydrogen atoms must undergo a parameterization step. For the use of an Amber force field in Gromacs simulations, a set of of Amber topology/parameter files are provided that need to be integrated into the respective force field's topology subdirectory of the Gromacs installation. Afterwards, the final PDB coordinates file constructed by this tool is immediately ready for Amber-parameterization and simulations using Gromacs. Charges were determined using the AM1-BCC method (AmberTool Antechamber) and were only slightly fitted in order to generally yield 0 formal polymer charge. As a consequence, no time-consuming recalculation of high-level partial atomic charges is necessary. The software also generates an image of the polymer graph using the Graphviz library.
 

How PolymerAssembler works
--------------------------

The main settings falling to the user are listed in the YAML configuration file `config_user.yml`. At the top of the file, the user choses the type of polymer. Currently supported types are:

- `branchedPG` (*n=5* units): linear, hyperbranched or dendritic polyglycerol polymers branched to a specified degree
- `linearPG_meth-eth` (*n=6* units): linear (via central glycerol oxygen) PG with methyl or ethyl (in random order) attached to the first oxygen of each monomer. Due to the chiral character of these glycerols central carbon atom, there are two methylized and two ethylized units. 

Depending on the choice of the polymer type, the transition probabilities from one building block to another need to be specified in the same config file. The size of the squared matrix is related to the number *n* of units. Let's consider a typical branched polyglycerol (PG) polymer. For PG, a set of five types of building blocks, {GCR,GCX,GCA,GCB,GCL} or {R,X,A,B,L} has been defined from which a new polymer of size *N* (number of monomers) can be assembled as a directed graph G(V,E) without cycles. This set consists of:

|Unit| function                          | indegree | outdegree |
|----|-----------------------------------|----------|-----------|
|GCR | designated root element of G      | 0        | 3         |
|GCX | branching unit                    | 1        | 2         |
|GCA | linearly extending unit of type A | 1        | 1         |
|GCB | linearly extending unit of type B | 1        | 1         |
|GCL | terminal (leaf) element           | 1        | 0         |

Using that set of five non-physical building blocks derived from glycerol it is possible to build PG polymers with any degree of branching ranging from linear and hyperbranched sequences up to 100% branched dendrimers. Starting with a single root element GCR, the central list L of available binding sites is initialized with three entries (due to indegree 0 and outdegree 3 associated with GCR). With each iteration over that list, another unit is attached to the current one of these sites and the list is updated accordingly. The decision of which unit to choose next is steered through a transition probability matrix P specifying the probability p_ij that block type j is attached to one of the free binding sites of block type i. As consequence, P also affects the extend to which the polymer will be branched. You can produce a linear polymer (apart from the beginning where two of the root element's three sites must be capped by terminal units GCL), various hyperbranched, and fully branched (=dendritic) PG polymers.

As an example, in case of an entirely branched polymer (dendrimer), P might look like this

`[0.0, 1.0, 0.0, 0.0, 0.0]`<br />
`[0.0, 1.0, 0.0, 0.0, 0.0]`<br />
`[0.0, 1.0, 0.0, 0.0, 0.0]`<br />
`[0.0, 1.0, 0.0, 0.0, 0.0]`<br />
`[0.0, 0.0, 0.0, 0.0, 0.0]`<br />

where the row as well as column order corresponds to the order of the units in the tably above. Each field *P_{ij}* sepcifies the probability with which the child building block of column *j* will be attached to the parent unit of row *i*. That is, due to *P_{iX}=1.0* (ones in the second column), any type *i* of the five units is always followed by the second unit type, the branching block X. The first column and last row must always be 0, since the root unit R has no predecessor (first column) and the terminal unit L has no successor (last row).

In case of a somehow hyperbranched PG, one would rather choose values such as
`transmatrix:`<br />
    `[0.00, 0.78, 0.10, 0.10, 0.02]`<br />
    `[0.00, 0.60, 0.15, 0.15, 0.10]`<br />
    `[0.00, 0.70, 0.10, 0.10, 0.10]`<br />
    `[0.00, 0.70, 0.10, 0.10, 0.10]`<br />
    `[0.00, 0.00, 0.00, 0.00, 0.00]`<br />

A linear polymer:
`transmatrix:`<br />
`[0.0, 0.0, 1.0, 0.0, 0.0]`<br />
`[0.0, 0.0, 1.0, 0.0, 0.0]`<br />
`[0.0, 0.0, 1.0, 0.0, 0.0]`<br />
`[0.0, 0.0, 1.0, 0.0, 0.0]`<br />
`[0.0, 0.0, 0.0, 0.0, 0.0]`<br />

The list L of unsatisfied binding sites may be worked off randomly or following the first in-first out principle resulting in highly spheric/symmetric polymers.

Further technical details are available in the following article that should also be used as a reference:

Vedat Durmaz. *Markov model-based polymer assembly from force field-parameterized building blocks*. Journal of Computer-Aided Molecular Design, 29:225-232, 2015.



Prerequisites
-------------

#### Required libraries

On a Ubuntu 16.04 Linux system, three additional libraries, `python-pygraphviz`, `python-biopython`, and `python-yaml` are required in order to be able to run the Python script. On such a Debian-based system one would execute

`sudo apt install python-pygraphviz python-biopython python-yaml`

to install both packages and all further dependencies.


#### Modify Gromacs topology files for usage with particular Amber force field

If you want to use the generated PDB coordinates file as input for an AMBER parameterization step using the `gmx pdb2gmx` command, you will need to modify a couple of files of that particular AMBER force field in Gromacs topology directory (usually in `/usr/share/gromacs/top`). It is then possible to quickly parameterize the polymer on the basis of preparameterized units as in analogy to polypeptides. You might want to copy the Gromacs topology directory to an own user directory before

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
