# MBR_Correlation
Code base for the SWAN data, to generate visibilities from tile data. 

Before we start, let get introduced to the SWAN system, the SWAN setup consists of seven tiles (ason Feb 15 2021), each tile
has a 4X4 bowtie antenna elements (MWA Tiles), the output is coherently added in the Analog beamformer and two Polarization
output (Linear and Vertical) are given as output. This has a collecting area (labda^2*gain of antenna * number of elements) 
16*0.5*1.5^2 = ~15 sq.m area at 150MHz.
