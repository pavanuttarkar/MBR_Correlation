# MBR_Correlation
Code base for the SWAN data, to generate visibilities from tile data. 

Before we start, let get introduced to the SWAN system, the SWAN setup consists of seven tiles (ason Feb 15 2021), each tile
has a 4X4 bowtie antenna elements (MWA Tiles), the output is coherently added in the Analog beamformer and two Polarization
output (Linear and Vertical) are given as output. This has a collecting area (labda^2 x gain of antenna x number of elements) 
16x0.5x1.5^2 = ~15 sq.m area at 150MHz. The two X and Y pols from the tile pass thorough an integrated anmplifier, consisting
of a low pass, high pass and a FM filter. This is sampled with 33 MHz sampling frequency for 16 MHz bandwidth, with 140 MHz
Intermediate Frequency. The rate of data acquisition in the MBR setup is 33MSPS, for individual Tiles, there are voltage sample
that are stored onto the disk.



Software based correlation is made using the SWAN MCU package, there are couple of steps before the correlation takes place,
these steps are af follows,


1. For time stamping, data from two different Tiles, a GPS 1PP signal is embedded in the header of the data (22:26 bytes),
   this is used to synchronize the data set from different tiles. The procedure for this is as follows.
   
   1.a The module \_header_call_to_read_cy, from the source file call_to_read_dev.pyx, consists of function 
