%Low Density Parity Check Codes (LDPC) MATLAB Simulation
%Binary Phase Shift Keying (BPSK) and Additive White Gaussian Noise Modulation
%(AWGN)
%Encoded by Modified Approximate Lower Triangular (MALT) Pre Processing
%Decoded by Sum Product (SPA) and Minimum Sum (MSA) Algorithms
%- Coded By,
      %Ganeshaanand (Rishi) Balasubramanian
      %MASc. Electrical and Computer Engineering
      %Dalhousie University
      %2018 - 2022
      
%------------------------------%----------------------------------%-----------------------------------%
%------------------------------MAIN PROGRAM-----------------------%

clear; clc;

%Fixes CPM matrix size affecting Parity Matrix Size
%Circular_Permutation_Matrix_Size = 31;

%Create the LDPC Parity Check Matrix
%Ha = Parity_Matrix_Gen(Circular_Permutation_Matrix_Size);

%-------The above are used to create specific design PCMs. For generalised,
%uncomment the below and commment the above.

[matro, Ha] = QCLDPC_PCM_Gen(3, 5, 2, 5, 31);    

%Rearrange PCM to MALT form for Encoding
H = MALT(Ha); 

figure;

%View of Pre Processed Parity matrix
imagesc(H);

%Get its size
[N1, N2] = size(H);

% Setting SNR decibel values for the Sim. You can run it for any values.
EbN0 = 1:1:5;

%Noise Variance Calculation
N0 = (10.^(EbN0/10));

nlen = length(EbN0);

% Size without parity bits per codeword
k = N2-N1;

%Total iterations
iter = 10;

%LDPC Code Rate
rate = k/N2;

%Create containers to store BER values, one for each decoder.
berate1 = zeros(1, length(EbN0));
berate2 = zeros(1, length(EbN0));

%Main Loop Begins, For every value of SNR
for i = 1:length(EbN0)
    
    %BER per SNR and total bits being transmitted
    ber1 = 0; ber2 = 0; nbits=0;
   
    %For each SNR, loop until 500 errors are accumulated. Change as needed.
    while(ber1 < 50)
        
        %Random data insertion
        y = randi([0 1], 1, k);
        
        %Linear Encoding Parity Bits Generation
        e = ldpcencode(H, y);    
        
        %Joining Info and parity to make an encoded code
        msg = [y e];
        
        %Now BPSK
        bpsk = 2 * msg - 1;    
        
        %Now AWGN
        tx = bpsk + sqrt(1 / (2 * rate * N0(i))) * randn(1, length(bpsk));
        
        %Bits are transmitted, accumulate count
        nbits = nbits + k;
        
        %Decode using SUM PRODCT
        vhat1 = SumProduct(tx, H, rate*N0(i), iter);
        ber1 = ber1 + sum(y ~= vhat1(1:k));
        
        %Decode using MIN SUM
        vhat2 = MinimumSum(tx, H, iter);
        ber2 = ber2 + sum(y ~= vhat2(1:k));
        
        %Live BER Gen, use if needed, just gives an indicator that it runs
        fprintf("\n Errors - %d,%d , Bits - %d, BER - %d, %d", ber1,ber2, nbits, ber1/nbits,ber2/nbits);
    end
    
        %Flush and combine values ready to plot
        berate1(i) = ber1/nbits; 
        berate2(i) = ber2/nbits; 
        
end 

%Plot
semilogy(EbN0, berate1, 'o--', EbN0, berate2, 'o-');
grid on;
hold off;

