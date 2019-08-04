function [ TIMF, R ] = timf( input_image, lmax, c)

%This is the code of FATEMD method
%Cited papers:
%1- Riffi, J., Adnane, M. M., Abbad, A., and Tairi, H. (2014). 3D extension of the fast and adaptive bidimensional empirical mode decomposition. Multidim. Syst. Signal Process. 26, 823–834. doi: 10.1007/s11045-014-0283-6
%2- He, Z., Li, J., Liu, L., and Shen, Y. (2017). Three-dimensional empirical mode decomposition (TEMD): a fast approach motivated by separable filters. Signal Process. 131, 307–319. doi: 10.1016/j.sigpro.2016.08.024
%3- Moradi, N., Dousty, M., Sotero, R. (2019). Spatiotemporal empirical mode decomposition of resting-state fMRI signals: application to global signal regression, Frontiers in Neuroscience, 13, 736. 

%Input arguments:
%  input_image - nifti format image
%  lmax - the maximum number of instrinsic mode functions to find
%  c - neighbourhood size when searching for local minima and maximum is a cube 2*c+1 voxels on a side.


timf_start_time=tic;
R = cast(input_image.img, 'double');
dim = size(R);
   
dim1 = dim(1);
dim2 = dim(2);
dim3 = dim(3);


ONE = ones (dim1 , dim2 , dim3);

%M = round((sqrt((dim1)^2 +(dim2)^2 +(dim3)^2 ))/2)   %margin

%%%
for l = 1: lmax,  %Number of IMF

  l_start_time=tic;
  message=[ 'Starting loop for l = ' num2str(l) ];
  disp(message)
    k = 0;
    g = 0;
    dminq1 = 0; 
    dmaxq1 = 0;       
   
    dminq2 = 0; 
    dmaxq2 = 0; 
    
    distq1 = [];
    distq2 = [];
    MINdistq2 = [1];
    MINdistq1 = [1];
    h2 = zeros ((dim1+2),dim2+2,(dim3+2));  % for max matrix
    h1= h2+999999;  %for min matrix, it should be more than 10^+4 (maximum data)
    
    
%     c=1;
     M = c;
       
    h1((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3)) = R;
    h2((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3)) = R;
    
    % h(x,y,z,1)=original data %R= residue of Timf

    % window size=3 , c=round(c/2)

  message=[ 'Starting local maximum and minimum loops.' ];
  disp(message)

  maxmin_start_time=tic;

    for x = (1+c) : (dim1+c) %90
        for y = (1+c) : (dim2+c)  %108
            for z = (1+c) : (dim3+c) , %90
                if (sum(sum(sum(max(max(max(h2(x-c : x+c , y-c : y+c , z-c: z+c)))))) == h2(x-c : x+c , y-c : y+c , z-c: z+c)))<2,
                     if (max(max(max(h2(x-c : x+c , y-c : y+c , z-c: z+c )))) == h2(x, y, z)), 
                    
                      
                    
                        g = g+1;
                    
                        x2(g) = x;
                        y2(g) = y;
                        z2(g) = z;
                        b=0;
                        dist2q = [];
                        for n = 1:g-1,
                        
                            b = b+1;
                            sqr2 = [x2(g), y2(g) ,z2(g)]-[x2(g-n) , y2(g-n) ,z2(g-n)];
                            dist2q(b) = (sqr2(1))^2 +(sqr2(2))^2 +(sqr2(3))^2;
                       
                        end 
                        if g>1,
                           MINdistq2(g) = min(dist2q); 
                        end
                     end
                end

                if  (sum(sum(sum(min(min(min(h1(x-c : x+c , y-c : y+c , z-c: z+c )))))) == h1(x-c : x+c , y-c : y+c , z-c: z+c)))<2,
                    if min(min(min(h1(x-c : x+c , y-c : y+c , z-c: z+c )))) == h1(x, y, z),
                        
                        
                   
                        k = k+1;
                    
                        x1(k) = x;
                        y1(k) = y;
                        z1(k) = z;
                        v = 0;
                        dist1q = [];
                        for m = 1:k-1,
                        
                            v = v+1;
                            sqr = [x1(k), y1(k) ,z1(k)]-[x1(k-m) , y1(k-m) ,z1(k-m)];
                            dist1q(v) = (sqr(1))^2 +(sqr(2))^2 +(sqr(3))^2 ;
                            
                        end
                        if k>1,
                            MINdistq1(k) = min(dist1q);
                        end
                    end  
                end
                
                
                
            end
        end
    end

    dminq1 = min(MINdistq1(2:end)); 
    dmaxq1 = min(MINdistq2(2:end));       %we can use one of the 4 formulate in the FATEMD paper
    
    d = sqrt(min(dminq1 , dmaxq1));
    
    if isempty(d) ==1,
        break
    else
        C_s = round(d/2)
        W_en = round((2*C_s)+1)
    end
    M = C_s

  toc(maxmin_start_time)

  message=[ 'Starting envelope loops.' ];
  disp(message)
  envelope_start_time1=tic;
    
  
   Zero = zeros((dim1+2*M),(dim2+2*M),(dim3+2*M));
   Zero((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3)) = ONE;

    U_Ej = zeros((dim1+2*M),dim2+2*M,(dim3+2*M));
    L_Ej = zeros ((dim1+2*M),dim2+2*M,(dim3+2*M));
    Min2 = zeros ((dim1+2*M),dim2+2*M,(dim3+2*M));
    Max2 = zeros((dim1+2*M),dim2+2*M,(dim3+2*M));
    %Max = zeros ((dim1+2*M),dim2+2*M,(dim3+2*M));
    
    %Min = zeros ((dim1+2*M),dim2+2*M,(dim3+2*M));
    
    h2 = zeros ((dim1+2*M),(dim2+2*M),(dim3+2*M));  % for max marrix
    h1= h2+999999;  %for min matrix, it should be more than 10^+4 (maximum data)
    
      
    h1((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3)) = R;
    h2((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3)) = R;
   
    
    %% Calculating Min & Max Envelope   
    
    for x = (M+1) : (M+dim1),
        for y = (M+1) : (M+dim2),
            for z= (M+1) : (M+dim3),
                  
                    Max2(x,y,z)=max(max(max(h2(x-C_s : x+C_s , y-C_s : y+C_s , z-C_s: z+C_s)))); % in a cube with size W_en*W_en*W_en
                    Min2(x,y,z)=min(min(min(h1(x-C_s : x+C_s , y-C_s : y+C_s ,z-C_s: z+C_s ))));
                
            end
        end
    end
  toc(envelope_start_time1)
    
  message=[ 'Starting mean envelope loops.' ];
  disp(message)
  envelope_start_time2=tic;
    %% smoothing and calculating the mean matrix(mean envelope) from smoothed MIN and Max envelopes(matrix)
    
    W_sm = W_en;
    Max2(1:(dim1+2*M) , 1:M , 1:(dim3+2*M)) = 0;
    Max2(1:(dim1+2*M) , (M+dim2+1):(dim2+2*M) , 1:(dim3+2*M)) = 0;
    Max2(1:M , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    Max2((M+dim1+1):(dim1+2*M) , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    Max2 ((M+1):(M+dim1), (M+1):(M+dim2), 1:M) =0;
    Max2 ((M+1):(M+dim1), (M+1):(M+dim2), (M+dim3+1):(dim3+2*M)) =0;
    
    Min2(1:(dim1+2*M) , 1:M , 1:(dim3+2*M)) = 0;
    Min2(1:(dim1+2*M) , (M+dim2+1):(dim2+2*M) , 1:(dim3+2*M)) = 0;
    Min2(1:M , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    Min2((M+dim1+1):(dim1+2*M) , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    Min2 ((M+1):(M+dim1), (M+1):(M+dim2), 1:M) =0;
    Min2 ((M+1):(M+dim1), (M+1):(M+dim2), (M+dim3+1):(dim3+2*M)) =0;
    
    
   for x = 1+C_s : (dim1+2*M)-C_s,
        for y = 1+C_s : (dim2+2*M)-C_s,
            for z= 1+C_s : (dim3+2*M)-C_s,
                U_Ej(x, y, z) =( sum(sum(sum((Max2(x-C_s:x+C_s,y-C_s: y+C_s, z-C_s: z+C_s)))))/(sum(sum(sum(Zero(x-C_s:x+C_s,y-C_s: y+C_s,  z-C_s: z+C_s))))));
                L_Ej(x, y, z) = (sum(sum(sum((Min2(x-C_s:x+C_s, y-C_s: y+C_s,  z-C_s: z+C_s)))))/(sum(sum(sum(Zero(x-C_s:x+C_s,y-C_s: y+C_s,  z-C_s: z+C_s))))));
                
            end
        end
   end
     
   
    U_Ej(1:(dim1+2*M) , 1:M , 1:(dim3+2*M)) = 0;
    U_Ej(1:(dim1+2*M) , (M+dim2+1):(dim2+2*M) , 1:(dim3+2*M)) = 0;
    U_Ej(1:M , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    U_Ej((M+dim1+1):(dim1+2*M) , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    U_Ej ((M+1):(M+dim1), (M+1):(M+dim2), 1:M) =0;
    U_Ej ((M+1):(M+dim1), (M+1):(M+dim2), (M+dim3+1):(dim3+2*M)) =0;
    
    
    L_Ej(1:(dim1+2*M) , 1:M , 1:(dim3+2*M)) = 0;
    L_Ej(1:(dim1+2*M) , (M+dim2+1):(dim2+2*M) , 1:(dim3+2*M)) = 0;
    L_Ej(1:M , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    L_Ej((M+dim1+1):(dim1+2*M) , 1:(dim2+2*M) , 1:(dim3+2*M)) = 0;
    L_Ej ((M+1):(M+dim1), (M+1):(M+dim2), 1:M) =0;
    L_Ej ((M+1):(M+dim1), (M+1):(M+dim2), (M+dim3+1):(dim3+2*M)) =0;
    
    
    M_Ej = ( U_Ej+L_Ej)/2;
  message=[ 'End of envelope calculations.' ];
  disp(message)
  toc(envelope_start_time2)
    
    %% making TIMF matrix by subtracting the original data from mean data
    TIMF1 = [];
    TIMF1(:,: ,:)= h2(:,: ,:)-M_Ej(:,: ,:);
    R(:,: ,:) = h2((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3)) - TIMF1((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3));
    TIMF(:,: ,:,l)= h2((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3))-M_Ej((M+1):(M+dim1), (M+1):(M+dim2), (M+1):(M+dim3));
   
  
  
message=[ 'End of loop l = ' num2str(l) ];
disp(message)
toc(l_start_time)
     
end    
   
message=[ 'End of timf function' ];
disp(message)
toc(timf_start_time)
end
