%- Coded By,
      %Ganeshaanand (Rishi) Balasubramanian
      %MASc. Electrical and Computer Engineering
      %Dalhousie University
      %2018 - 2022
      
%------------------------------%----------------------------------%-----------------------------------%
%------------------------------Modified Approximate Lower Triangular (MALT)-----------------------%


function ALTF = MALT(H)
    AT = ALT(H); %Makes lower tri
    imagesc(AT);
    ATE = gjetest(AT); %clears bits below lower tri
    imagesc(ATE);
    A = uD(ATE); %Clears zero rows
    imagesc(ATE);
    ALTF = A; %Final
    imagesc(ATE);
    
    function AT = ALT(H)
        [m, n] = size(H);

        count = sum(H);
        count(count==0) = nan;
        [~, mincol] = min(count) ;                              %Find first column with least # of 1s

        temp = H(:,mincol);
        H(:,mincol) = [];
        H = [H(:, 1:end) temp];                                 %Move column to end of matrix

        pero = find(H(:,n)==1);
        tero = H(pero,:);                                       %move rows with 1s in column to end of matrix.
        H(pero,:)=[];
        H = [H(1:end, :); tero];
    
        Tc = find(H(:, n), 1, 'first');                         %current position of diagonal.
        T = Tc-1;                                               %updated diagonal position
        p = n-1;                                                %next column
                                                          %Begin loop
        while T>=1
        
            count = sum(H(1:T, 1:p));
            count(count==0) = nan;
            [~, mincol] = min(count) ;                          %choose next column with least number if 1s 
            temp = H(:, mincol);
            H(:, mincol) = [];
            H = [H(:, 1:p-1) temp H(:, p:end)];                 %move column to index p
     
            [r, ~] = size(find(H(1:T,p)==1));                   %number of ones in column p
            d = H(T,p);                                         %is 1 present in current diagonal
   
            switch true
                case d==1 && r==1                               %if 1 is in current diagonal and number of ones is 1
                    disp("Noice");
            
                case d==1 && r>1                                %if 1 is in current diagonal and number of ones is more than 1
                    col = p;
                    rowsOnly = 1:T-1; 
                    pa = find(H(1:T-1, col)==1);
                    ta = H(pa,:);
                    H(pa,:)=[];
                    H=[H(1:end,:); ta];
            
                case d==0 && r==1                               %if 1 is not in current diagonal and number of ones is 1
                    col = p;
                    rowsOnly = 1:T; 
                    pc = find(H(1:T,col)==1);
                    idxRows = ismember(1:size(H, 1), rowsOnly)'; 
                    idx0 = H(:, col) < 1 & idxRows; 
                    idx1 = H(:, col) > 0 & idxRows;
                    H = [H(idx0, :); H(idx1, :); H(~idxRows, :); ];
            
                case d==0 && r>1                                %if 1 is not in current diagonal and number of ones is more than 1
                    col = p;
                    ro = find(H(:, col), 1, 'first');
                    te = H(ro,:);
                    H(ro,:)=[];
                    H = [H(1:T-1,:); te; H(T:end,:);];
        
                    col = p;
                    rowsOnly = 1:T;
                    pp = find(H(1:T-1, col)==1);
                    tt = H(pp,:);
                    H(pp,:)=[];
                    H=[H(1:end,:); tt];
    
            end
            %imagesc(H)
            T = find(H(:, p), 1, 'first');
            T = T-1;                                            %updated diagonal position
            p = p-1;
        end
    
        AT = H;
    end
    
    function ATE = gjetest(A)
Hi = A;
 
[m, n] = size(Hi);
g = m - (find(Hi(:, n), 1, 'first'));    %Find gap of Matrix

dr = m-g+1;                             %First row of D matrix
dc1 = n-m+1;   
rm = m-g;                               %Row marker
row = m-g;
col = n-m+g+1;

for a = 1:m-g
    
    rta = rm + find(Hi(dr:end, n)==1);
    
    for i = 1:length(rta)
        Hi(rta(i), :) = xor(Hi(rta(i), :), Hi(row, :));     %Add all subsequent rows with current row
    end
    %imagesc(Hi);
    row = row-1;      %next row
    n = n-1;
    
end

ATE = Hi;
    end

    function B = uD(H)
    [m, n] = size(H);
    g = m - (find(H(:, n), 1, 'first'));    %Find gap of PCM
    Mat = H(1:end, 1:n-m+g);   
    
    dr = m-g+1;                             %First row of D matrix
    dc = n-m+1;                            %First col of D matrix
    zeron = 0;
    
    for j = dr : m    
        chk = find(Mat(dr, dc));  
    
        if (chk == 1)
                disp("G");                                      %Subsequent rows only
            
        else
            colo = find(Mat(dr, 1:end), 1, 'first');
            if (length(colo)>0)
                col = Mat(:, colo);                                     %Extract that col
                Mat(:, colo) = Mat(:,dc);                                      %Empty that col in Matrix
                Mat(:,dc)=col;                                             %Move that column to required position
        
            else
                Mat(dr, :) = [];
                zeron = zeron+1;
                m=m-1;
                continue;
            end
        
        end
    
        mark = find(Mat(1:end, dc) == 1);                %Find all rows with 1 at same column
        mark(mark<=dr)=[];
   
        for i = 1:length(mark)
                Mat(mark(i), :) = xor(Mat(mark(i), :), Mat(dr, :));     %Add all subsequent rows with current row
        end
    
        dr=dr+1;
        dc=dc+1;
    
    
    end
    
    rtm = Mat(:, end-zeron+1 : end);  %Identify the redundant columns after Upper Triangular D matrix
    Mat(:, end-zeron+1 : end) = [];   
    Mat = [Mat(:, 1:end-g+zeron) rtm Mat(:, end-g+zeron+1:end)];  %Move Columns
    
    C=Mat;
    [m,n] = size(C);
    
    B = [C H(1:m, n+1:end)]; 

end


end
