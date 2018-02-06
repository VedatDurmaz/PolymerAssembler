# PolymerAssembler

PolymerAssembler is a Python-based tool designed for a quick construction of atomistic 3D polymer models useful for molecular mechanics force field simulations. Particularly, in case of branched polymers the degree of branching is influenceable through the user. The resulting PDB file initially lacks hydrogen atoms which can be easily added upon parameterization using the Gromacs command `pdb2gmx` along with its integrated AMBER force field. For this purpose, a set of files with AMBER topology/parameter definitions are provided that need to be added to the (corresponding force field's) topology (sub)directory of your local Gromacs installation. The parameterization procedure based on pre-parameterized building blocks is analogue to that of proteins based on aminoacid blocks. After an energy minimization step, the final PDB file is immediately ready for molecular dynamics simulations using Gromacs. Partial atomic charges of the units were determined using the AM1-BCC method (as implemented in the AmberTool `Antechamber`) and were only negligibly modified in order to generally fit the polymer's overall formal charge to zero. With this strategy, no time-consuming calculation of high-level partial atomic charges is necessary. The program also generates an image of the polymer graph using the Graphviz library. In its current version PolymerAssembler supports the following types of polymers:

- Polyglycerol: ranging from linear to hyperbranched and dendritic polymers
- Polyglycerol: only linear and with methyl and/or ethyl attached to the hydroxy group
- Polyethylene oxide: linear (via central glycerol oxygen) PG with methyl and/or ethyl (in random order) attached to the first oxygen of each monomer



How PolymerAssembler works
--------------------------

The main settings falling to the user are listed in a YAML configuration file called `config_user.yml`. At the top of the file, the user choses the type of polymer by setting the YAML variable `polymer.type` to either

|Type                | #blocks | description |
|:-------------------|:--------|:------------|
|`branchedPG`        |    5    | linear/hyperbranched/dendritic polyglycerol |
|`linearPG-meth-eth` |    6    | methyl/ethyl polyglycerol |
|`linearPEO`         |    3    | linear polyethylene oxide |

Depending on the polymer type choice, the corresponding transition probabilities from one building block to another need to be adjusted in the same config file according to the user's needs. The size of that squared matrix is related to the number *n* of units. Let's consider a typical branched polyglycerol (PG) polymer. For PG, a set of five types of building blocks, {GCR,GCX,GCA,GCB,GCL} or {R,X,A,B,L} has been defined from which a new polymer of size *N* (number of monomers) can be assembled as a directed graph G(V,E) without cycles. This set consists of:

|Unit| function                          | indegree | outdegree |
|----|-----------------------------------|----------|-----------|
|GCR | designated root element of G      | 0        | 3         |
|GCX | branching unit                    | 1        | 2         |
|GCA | linearly extending unit of type A | 1        | 1         |
|GCB | linearly extending unit of type B | 1        | 1         |
|GCL | terminal (leaf) element           | 1        | 0         |

Using that set of five non-physical building blocks derived from glycerol it is possible to build PG polymers with any degree of branching ranging from linear and hyperbranched sequences up to 100% branched dendrimers. Starting with a single root element GCR, the central list L of available binding sites is initialized with three entries (due to indegree 0 and outdegree 3 associated with GCR). With each iteration over that list, another unit is attached to the current one of these sites and the list is updated accordingly. The decision of which unit to choose next is steered through a transition probability matrix P specifying the probability p_ij that block type j is attached to one of the free binding sites of block type i. As consequence, P also affects the extend to which the polymer will be branched. You can produce a linear polymer (apart from the beginning where two of the root element's three sites must be capped by terminal units GCL), various hyperbranched, and fully branched (=dendritic) PG polymers.

As an example, in case of an entirely branched polymer (dendrimer), P might look like this
```
transmatrix:
    - [0.0, 1.0, 0.0, 0.0, 0.0]
    - [0.0, 1.0, 0.0, 0.0, 0.0]
    - [0.0, 1.0, 0.0, 0.0, 0.0]
    - [0.0, 1.0, 0.0, 0.0, 0.0]
    - [0.0, 0.0, 0.0, 0.0, 0.0]
```
where the row as well as column order corresponds to the order of the units in the tably above. Each field *P_{ij}* sepcifies the probability with which the child building block of column *j* will be attached to the parent unit of row *i*. That is, due to *P_{iX}=1.0* (ones in the second column), any type *i* of the five units is always followed by the second unit type, the branching block X. The first column and last row must always be 0, since the root unit R has no predecessor (first column) and the terminal unit L has no successor (last row).

In case of a somehow hyperbranched PG, one would rather choose values such as
```
transmatrix
    - [0.00, 0.78, 0.10, 0.10, 0.02]
    - [0.00, 0.60, 0.15, 0.15, 0.10]
    - [0.00, 0.70, 0.10, 0.10, 0.10]
    - [0.00, 0.70, 0.10, 0.10, 0.10]
    - [0.00, 0.00, 0.00, 0.00, 0.00]
```
A linear polymer:
```
transmatrix:
    - [0.0, 0.0, 1.0, 0.0, 0.0]
    - [0.0, 0.0, 1.0, 0.0, 0.0]
    - [0.0, 0.0, 1.0, 0.0, 0.0]
    - [0.0, 0.0, 1.0, 0.0, 0.0]
    - [0.0, 0.0, 0.0, 0.0, 0.0]
```

The list L of unsatisfied binding sites may be worked off randomly or following the first in-first out principle resulting in highly spheric/symmetric polymers.

Further technical details are available in the following article that should also be used as a reference:

Vedat Durmaz. *Markov model-based polymer assembly from force field-parameterized building blocks*. Journal of Computer-Aided Molecular Design, 29:225-232, 2015.



Prerequisites
-------------

#### Required libraries

On a Ubuntu 16.04 Linux system, three additional libraries, `python-pygraphviz`, `python-biopython`, and `python-yaml` are required in order to be able to run the Python script. On such a Debian-based system one would execute

`sudo apt install python-pygraphviz python-biopython python-yaml`

to install these packages and all of their dependencies.


#### Extend Amber force field in Gromacs topology directory

If you want to use the generated PDB coordinates file as input for an AMBER parameterization step using the `gmx pdb2gmx` command, you will need to modify a couple of files of that particular AMBER force field in Gromacs topology directory (usually in `/usr/share/gromacs/top`). It is then possible to quickly parameterize the polymer on the basis of preparameterized units as in analogy to polypeptides. You might want to copy the Gromacs topology directory to an own user directory before

`cp -a /usr/share/gromacs/top ~/my-gromacs-top`

and set the corresponding environment variable for Gromacs commands

`export GMXLIB=/home/$USER/my-gromacs-top`

Assuming, we want to use the Amber99sb force field, the following files would be extended by the content of the respective files in the "gmx_topology" directory of PolymerAssembler:
```
cd PolymerAssembler
cat ./pdb_units/${poltype}/gmx_topology/residuetypes.dat >> ${GMXLIB}/residuetypes.dat
cat ./pdb_units/${poltype}/gmx_topology/specbond.dat >> ${GMXLIB}/specbond.dat
cat ./pdb_units/${poltype}/gmx_topology/aminoacids.hdb >> ${GMXLIB}/amber99sb.ff/aminoacids.hdb
cat ./pdb_units/${poltype}/gmx_topology/aminoacids.rtp >> ${GMXLIB}/amber99sb.ff/aminoacids.rtp
```
In addition, increase the number in the first line of `${GMXLIB}/specbond.dat` by the number of newly added special bonds (lines), e.g. 32 in case of polyglycerol (polymer type "branchedPG" coming with five building blocks).


Usage
-----

Example command creating a polymer consisting of 50 units 

`../polymer.py 50 polymer.pdb polymer.eps`

where a corresponding pdb file and a graph image in the eps format are generated. 

To be continued ...
