# G-CYPWare-1.0
This repository contains a user friendly web-tool (cypwarecode.html) to be used with G-CYPWare-1.0 program for calculating CYP450 ET parameters
**-- G-CYPWare 1.0 script and HTML Tool --**

Please cite the use of CYPWare 1.0 and the HTML Tools as follows:
CYPWare-1.0 used to calculate the to calculate Marcus ET parameters for Cytochrome P450 BM3 reactions using the AMBER SOFTWARE. Those who are familiar with amber and 
having licence for the amber can use the CYPWare 1.0. here is the liink.
https://github.com/Dixit-s-lab/CYPWare-1.0
Dixit, V. A.; USN Murty, Bajaj, P.; Blumberger, J.; and Sam P. de Visser. Mechanisms of Electron Transfer Rate Modulations in Cytochrome P450 BM3. J. Phys. Chem. B. 2022 126 (47), 9737-9747. DOI: 10.1021/acs.jpcb.2c03967
https://pubs.acs.org/doi/10.1021/acs.jpcb.2c03967

WELCOME to G-CYPWare 1.0 web interface to calculate Marcus ET parameters for Cytochrome P450 BM3 reactions.
This tool is developed by Pulkit Asati (JRF NSM GAP-147) under the supervision of the PI: Dr. Vaibhav A. Dixit, Asst. Prof., Dept. of Med. Chem.,
NIPER Guwahati in the Advanced Center of Computer-Aided Drug Design (A-CADD).
CYPWare software, GUI, websites, tools, layouts, and logos are subject to copyright protection and are the exclusive
property of the Department of Medicinal Chemistry, NIPER Guwahati.
The name and logos of the NIPER Guwahati website must not be associated with publicity or business
promotion without NIPER G's prior written approval.

G-CYPWare 1.0 is available under the creative commons license and is free for academic research groups working
in a degree-granting university/institute.
Any work/report/thesis/research-article/review-article resulting from the use of G-CYPWare 1.0 should properly
cite the software and publication associated with the same.

#-- Dependencies --

G-CYPWare 1.0 makes use of the following opensource tools, thus these need to be installed first from GitHub,
Sourceforge, or as a conda package.
Ensure that dependencies are satisfied before running CYPWare 1.0, else the calculation will not complete as expected.

1) Openbabel 3.1 or higher

Openbabel is available as a conda package and can be installed with one of the following commands.

conda install -c conda-forge openbabel
conda install -c conda-forge/label/cf202003 openbabel

If you don't have conda, then install it from the main website https://www.anaconda.com/
Instructions for conda installation can be found on its website.

2) AmberTools


AmberTools is a freely available software used for the setup and analysis of MD simulations.
It is available from the http://ambermd.org/ website. It can also be installed as a conda package with any one of the following command.
conda install -c conda-forge ambertools
conda install -c conda-forge/label/cf202003 ambertools

3) GROMACS is a molecular dynamics package mainly designed for simulations of proteins, lipids, and nucleic acids.
It was originally developed in the Biophysical Chemistry department of University of Groningen,
and is now maintained by contributors in universities and research centers worldwide.

GROMACS is one of the fastest and most popular software packages available, and can run on central
processing units (CPUs) and graphics processing units (GPUs).[9] It is free, open-source software
released under the GNU General Public License (GPL),[3] and starting with version 4.6, the GNU Lesser
General Public License (LGPL). https://gitlab.com/gromacs/gromacs

For commercial usage of G-CYPWare 1.0, please contact the PI at vaibhavadixit@gmail.com or vaibhav@niperguwahati.in
========================================================================================================================
**Funding Sources**
VAD acknowledges the financial support from the National Supercomputing Mission (NSM), Department of Science and Technology (DST), New Delhi 
