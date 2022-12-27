# TrackFluoNGFvesicles
Scripts in MatLab to track single vesicles containing fluoNGF

Code to execute:
-----------------------------------------------------------------
Fase1L %in the folder with the videos to be analyzed

Fase2newLuin %in the folder with the output files from Fase1L

Fase3Lnew %in the folder with the videos to be analyzed
%load file "Info", e.g. by "load('Info')"

DomenicaFase4

Classifytraj 

------------------------------------------------------------------
The code require the OME Bio-Format package, revision 33bb1150.
Specify in ffpath.m the location of loci_tools.jar.
Part of the code present in the folder is taken or adapted from the Matlab code coming with the OME Bio-Format package.

The folder contains the packet from Raghuveer Parthasarathy, "Rapid, accurate particle tracking by calculation of radial symmetry centers", Nat. Methods 9, 724–6 (2012).

Part of the code in this folder is redistributed under GNU license (gpl-3.0.txt).
You can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

See the single files for more detailed Copyright information.

------------------------------------------------------------------
Part of the code was already used in:
-	“Precursor and mature NGF live tracking: one versus many at a time in the axons”, T. De Nadai, L. Marchetti, C. Di Rienzo, M. Calvello, G. Signore, P. Di Matteo, F. Gobbo, S. Turturro, S. Meucci, A. Viegi, F. Beltram, S. Luin and A. Cattaneo. Scientific Reports 6, 20272 (2016); doi: 10.1038/srep20272.
-	“Graphene promotes axon elongation through local stall of Nerve Growth Factor signaling endosomes”. D. Convertino, F. Fabbri, N. Mishra, M. Mainardi, V. Cappello, G. Testa, S. Capsoni, L. Albertazzi, S. Luin, L. Marchetti, and C. Coletti, Nano Lett. 20(5), 3633–3641 (2020) acs.nanolett.0c00571.
-	“Lysosome Dynamic Properties during Neuronal Stem Cell Differentiation Studied by Spatiotemporal Fluctuation Spectroscopy and Organelle Tracking”. W. Durso, M. Martins, L. Marchetti, F. Cremisi, S. Luin, F. Cardarelli, Int. J. Mol. Sci. 21(9) 3397 (2020). Doi: 10.3390/ijms21093397.

------------------------------------------------------------------
More details are available upon request to Stefano Luin <s.luin@sns.it>
