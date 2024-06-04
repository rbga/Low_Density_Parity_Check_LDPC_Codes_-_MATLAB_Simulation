# Low Density Parity Check (LDPC) Codes MATLAB-Simulation

LDPC MATLAB simulation using BPSK + AWGN modulation decoded using Sum Product and Min Sum Algorithm.

## Overview

This repository contains the MATLAB simulation for Low Density Parity Check (LDPC) codes using Binary Phase Shift Keying (BPSK) modulation and Additive White Gaussian Noise (AWGN) channel, decoded using Sum Product Algorithm (SPA) and Min Sum Algorithm (MSA).

## Description

This simulation is part of my Master's Thesis in Electrical and Computer Engineering. The thesis is currently available on [this link](http://hdl.handle.net/10222/81991) from September 2023. Check out the full thesis manuscript.

## Simulation

The MATLAB simulation involves the following steps:

1. Creation of Quasi Cyclic LDPC codes using input parameters `j`, `k`, `m`, `a`, and `b`.
2. Application of MALT form rearrangement function on the generated matrices, which displays graphs of the original and rearranged matrices.
3. SNR values are specified from 0 to 5 dB, and the number of iterations is set as `N2` with the PCM as [N1, N2].
4. For each iteration of SNR and frame, the message bit is encoded, BPSK modulated, passed through AWGN, and then decoded using SPA and MSA decoders in parallel. The Bit Error Rate (BER) for each iteration is accumulated.
5. Once the loop is complete, the BER for SPA and MSA are plotted against the SNR values in dB.

## Performance and Results

The simulation achieves a BER in the range of 10^-6. The software was tested on a Windows 10 Intel Dual Core i5 - 7200 of 2.5GHz, and simulations were run on all available cores simultaneously using Parallel Pool. LDPC codes of lengths in the range of 100s showed faster preprocessing and overall simulation. However, codes of higher lengths took longer processing times. For instance, codes above a length of 1000s required anywhere from 4 to 8 hours to complete. During the thesis period, the simulation was also run on an Intel Xeon E3 with 3.5GHz, resulting in significantly reduced processing time, from 6 hours to approximately 30 minutes in total.

The simulation results are displayed in the repository.

Feel free to explore the MATLAB code and run the simulation using the provided main program for a detailed understanding of LDPC codes and their performance characteristics.

Please check back on the link provided above for the full thesis manuscript, which will contain a comprehensive discussion of the LDPC simulation and its implications.

For any inquiries or feedback, you can contact me at [rbga@dal.ca].


![image](https://user-images.githubusercontent.com/75168756/213565856-6253a0c5-ee36-4f54-ab84-b02fa6483674.png)
