# Low-Density-Parity-Check-Codes---MATLAB-Simulation
LDPC MATLAB simulation using BPSK + AWGN modulation decoded using Sum Product and Min Sum Algorithm.

This is one part of my Masters Thesis in Electrical and Computer Engineering. Thesis on Embargo, Coming soon. Check back on this link for the full thesis manuscript from Sep 2024
http://hdl.handle.net/10222/81991

Open Main Program and Run for easy execution. Read below for subject matter.

For the simulation, a function with input parameters for j, k, m, a and b creates the
specified Quasi Cyclic LDPC code. After which that matrix is run on the MALT
form rearrangement function which uses row column permutations for each of them
iteratively, which will display the graphs of both the original ad rearranged matrices.
The SNR values are specified from 0 to 5 dB and the frame is set for N2 iterations
where the PCM is [N1,N2]. For each iteration of the SNR and Frame the message bit
is encoded, BPSK modulated, AWGN passed and decoded by SPA and MSA decoders
in parallel, as the Bit Error Rate for each iteration is accumulated. Once the loop
is complete, the BER for SPA and MSA are plotted against the SNR values in dB.
This method was able to achieve BER in the range of 10−6.
An entire MATLAB Simulation was created to simulate the transmission of data. The
software was used on a Windows 10 Intel Dual Core i5 - 7200 of 2.5GHz. Simulations
were run on all available cores simultaneously using Parallel Pool. LDPC codes in
the range of 100’s had a faster pre processing run and overall simulation. Codes of
higher length took some processing. For instance, codes above the length of 1000’s
took anywhere from 4 - 8 hours. During the thesis period, the simulation was run on
an Intel Xeon E3 with 3.5GHz, which had a higher speed for the same range, cutting
down the time from 6 hours to mere 30 minutes in total. Results are displayed below.

![image](https://user-images.githubusercontent.com/75168756/213565856-6253a0c5-ee36-4f54-ab84-b02fa6483674.png)
