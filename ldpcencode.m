%- Coded By,
      %Ganeshaanand (Rishi) Balasubramanian
      %MASc. Electrical and Computer Engineering
      %Dalhousie University
      %2018 - 2022
      
%------------------------------%----------------------------------%-----------------------------------%
%------------------------------Linear Encoder for LDPC-----------------------%


function p=ldpcencode(H,msg)
[m, n] = size(H);
% Find the 'gap' length

g = m - find(H(:,end),1, 'first');
% Extracting  the submatrices A, B, C, D, E and T
A = H(1:m-g,1:n-m);
B = H(1:m-g,n-m+1:n-m+g);
T = H(1:m-g,n-m+g+1:end);
C = H(m-g+1:end,1:n-m);
D = H(m-g+1:end,n-m+1:n-m+g);
E = H(m-g+1:end,n-m+g+1:end);

invT = inv(T); 
invD = inv(D);

%From Cu + Dp1 + 0p2 = 0
p1 = mod((-invD * C * msg'), 2);
%p1 = mod((C*msg'), 2);
%From Au + Bp1 + Tp2 = 0
p2 = mod(-invT * (A * (msg') + B * (p1)),2)';

p = [p1' p2];