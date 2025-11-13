clear all;
load('Fl_Tib');
data = Fl_Tib; 

%site_out = [3 9 18 30 43 45 49 51 53]; %Sites removed from Muskingum data
%data(:,[site_out]) = [];
[nr,nc]=size(data);

for j= 1:nc
    index = find(data(:,j)>0);
    flows = data(index,j);
    orderedflows = sort(flows,1);
    b0 = mean(orderedflows);
    
    n(j) = size(flows,1);
    aux = 1:n(j);
    aux2 = (aux - 1)/(n(j)*(n(j)-1));
    b1 = sum(aux2'.*orderedflows);
    
    aux3 = max(0,aux - 2);
    b2 = sum(aux2'.*aux3'/(n(j)-2).*orderedflows);
    
    lambda1 = b0;
    lambda2 = 2*b1 - b0;
    lambda3 = 6*b2 - 6*b1 + b0;
    
    tao3i(j) = lambda3/lambda2;
        
    c = (2*b1 - b0)/(3*b2 - b0) - log10(2)/log10(3);
        
    kappa(j) = 7.859*c + 2.9554*c^2;
        
    for jj = 1:nc
        indexequal = find(data(:,j) ~= -999 & data(:,jj) ~= -999); % find years both sites have data
        nij(j,jj) = size(indexequal,1);
    end    
    clear index;
    clear flows;
    clear orderedflows;
    clear aux;
    clear aux2;
    clear aux3;
end
cfij = nij./sqrt(n'*n);
%plot(tao3)
%hold on
plot(kappa)
tao3 = mean(tao3i);
c = 2/(3 + tao3) - log(2)/log(3);
K = 7.859*c + 2.9554*c^2;
Na = mean(n);

% Compute Var(tao3) and Var(kappa). See Chowdhury et al. (1991) and Phill's PhD Thesis (Appendix A).

r = gamma(1 + 2*K)/(gamma(1 + K))^2;
aux = (1 - 2^(-K))^2;
[hyper05] = feval('mhygfx',K,2*K,1+K,-0.5);
[hyper23] = feval('mhygfx',K,2*K,1+K,-2/3);
[hyper13] = feval('mhygfx',K,2*K,1+K,-1/3);

f11 = (r - 1)/aux;
f22 = (r*hyper05 - 1)/(2^(2*K)*aux);
f33 = (r*hyper23 - 1)/(3^(2*K)*aux);
f12 = 0.5*(r - 2^(1 + K) + 2^(2*K))/(2^(2*K)*aux);    
f13 = 0.5*((r - 2*3^K)/(3^(2*K)*aux) - (r*hyper05 - 2^(1 + K))/(2^(2*K)*aux));    
f23 = 0.5*((r*hyper13 - 2*(3/2)^K)/(3^(2*K)*aux) + 1/(2^(2*K)*aux));

var_tao3 = 1/Na*((f11 - 4*f12 + 4*f22)*tao3^2 + 2*(f11 - 8*f12 + 12*f22 + 6*f13 - 12*f23)*tao3 + f11 - 12*f12 + 36*f22 + 12*f13 - 72*f23 + 36*f33);

dK_dtao3 = -(2.872/(tao3 + 3))^2*(2.872/(tao3 + 3) + 1);

var_K = dK_dtao3^2*var_tao3;

var_tao3i = (Na./n)*var_tao3;

result = [tao3i' var_tao3i' var_tao3i'.^0.5];

nvartao3 = Na*var_tao3;

%%Correlation

%%Tibagi
load('TCov0');
for i = 1:length(n)
    for j=1:nc
        if Cov0(i,j) ~= 0
           corr(i,j) = (Cov0(i,j)/cfij(i,j))^(1/2.85);
       end
   end
end   
corr = cfij.*corr.^2.8;

%%Illinois
%load('Cov0Ill');
%for i = 1:length(n)
%    for j=1:nc
%        if Cov0(i,j) ~= 0
%           corr(i,j) = (Cov0(i,j)/cfij(i,j))^(1/3);
%       end
%   end
%end   
%corr = cfij.*corr.^2.8;

%% Muskingum
%corr = cfij.*(0.26)^2.7; %Muskingum

for i = 1:length(n)
   corr(i,i) = 1;
end   
varmatrix = sqrt(var_tao3i')*sqrt(var_tao3i);
varmatrix2 = (var_tao3i').^0.5*(var_tao3i).^0.5;
var12 = varmatrix./varmatrix2;
scovmatrix2 = varmatrix.*corr;
save('TCov0tao3','scovmatrix2','n','tao3i');