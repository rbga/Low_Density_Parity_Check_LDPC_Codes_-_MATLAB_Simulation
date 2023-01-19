% Coded by Bagawan S. Nugroho, 2007 
% http://bsnugroho.googlepages.com
%Bagawan Nugroho (2023). LDPC codes BER simulation (https://www.mathworks.com/matlabcentral/fileexchange/14869-ldpc-codes-ber-simulation), MATLAB Central File Exchange.
function vHat = SumProduct(rx, H, N0, iteration)

[M, N] = size(H);

% Prior log-likelihood. 2y/sigma^2 where sigma = sqrt(N0/2)
Lci = (-4*rx/N0);

% Initialization
Lrji = zeros(M, N);
Pibetaij = zeros(M, N);

% Asscociate the L(ci) matrix with non-zero elements of H
Lqij = H.*repmat(Lci, M, 1);
 
% Get non-zero elements
[r, c] = find(H);

% Iteration
for n = 1:iteration
   
  % fprintf('Iteration : %d\n', n);
   
   % Get the sign and magnitude of L(qij)   
   alphaij = sign(Lqij);   
   betaij = abs(Lqij);

   for l = 1:length(r)
      Pibetaij(r(l), c(l)) = log((exp(betaij(r(l), c(l))) + 1)/...
                             (exp(betaij(r(l), c(l))) - 1));
   end
   
   % ----- Horizontal step -----
   for i = 1:M
      
      % Find non-zeros in the column
      c1 = find(H(i, :));
      
      % Get the summation of Pi(betaij))        
      for k = 1:length(c1)

         sumOfPibetaij = 0;
         prodOfalphaij = 1;
         
         % Summation of Pi(betaij)\c1(k)
         sumOfPibetaij = sum(Pibetaij(i, c1)) - Pibetaij(i, c1(k));
         
         % Avoid division by zero/very small number, get Pi(sum(Pi(betaij)))
         if sumOfPibetaij < 1e-20
            sumOfPibetaij = 1e-10;
         end         
         PiSumOfPibetaij = log((exp(sumOfPibetaij) + 1)/(exp(sumOfPibetaij) - 1));
      
         % Multiplication of alphaij\c1(k) (use '*' since alphaij are -1/1s)
         prodOfalphaij = prod(alphaij(i, c1))*alphaij(i, c1(k));
         
         % Update L(rji)
         Lrji(i, c1(k)) = prodOfalphaij*PiSumOfPibetaij;
         
      end % for k
      
   end % for i

   % ------ Vertical step ------
   for j = 1:N

      % Find non-zero in the row
      r1 = find(H(:, j));
      
      for k = 1:length(r1)        
        
         % Update L(qij) by summation of L(rij)\r1(k)
         Lqij(r1(k), j) = Lci(j) + sum(Lrji(r1, j)) - Lrji(r1(k), j);
      
      end % for k
      
      % Get L(Qi)
      LQi = Lci(j) + sum(Lrji(r1, j));
      
      % Decode L(Qi)
      if LQi < 0
         vHat(j) = 1;
      else
         vHat(j) = 0;
      end
                       
   end % for j
   
   %If arrived at proper codeword then end
   cs = mod(vHat*H',2);
    if sum(cs)== 0  
        break;
    end
    
    if vHat==0
        rx(rx>=0)=0;
        rx(rx<0)=1;
        vHat = rx;
    end
    
end % for n
