%- Coded By,
      %Ganeshaanand (Rishi) Balasubramanian
      %MASc. Electrical and Computer Engineering
      %Dalhousie University
      %2018 - 2022
      
%------------------------------%----------------------------------%-----------------------------------%
%------------------------------QC LDPC Specific Type Parity Matrix Generator-----------------------%
%This function generates a particular type of LDPC matrix without short
%cycles.
function Ha = Parity_Matrix_Gen(cpm)

I0 = zeros(cpm, cpm);
I1 = eye(cpm);

 %I0 = circshift(I1, -0);  
 I22 = circshift(I1, -22);  
 I28 = circshift(I1, -28);  
 I2 = circshift(I1, -2);  
 I16 = circshift(I1, -16);  
 I19 = circshift(I1, -19);  
 I30 = circshift(I1, -30);  
 I4 = circshift(I1, -4);  
 %I0 = circshift(I1, -0);  
 I19 = circshift(I1, -19);  
 I11 = circshift(I1, -11);  
 I10 = circshift(I1, -10);  
 I26 = circshift(I1, -26);  
 I6 = circshift(I1, -6);  
 I18 = circshift(I1, -18); 
 
 
 Ha = [I0 I22 I28 I2 I16 
I19 I30 I4 I0 I19 
I11 I10 I26 I6 I18 ];
  


end