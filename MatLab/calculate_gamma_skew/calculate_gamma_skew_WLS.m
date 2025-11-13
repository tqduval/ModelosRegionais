function [checkGamma,gamma,Parameter_WLS,lambda,Ilambda] = calculate_gamma_skew_WLS(gamma,N,scovmatrix,skew,A,p)
diag_gamma = zeros(N);
for i= 1:N
   diag_gamma(i,i) = gamma^2;
end
lambda = diag_gamma + scovmatrix; 
Ilambda = inv(lambda);
Parameter_WLS = inv(A'*Ilambda*A)*A'*inv(lambda)*skew;
checkGamma = (skew-A*Parameter_WLS)'*Ilambda*(skew-A*Parameter_WLS) - (N-p);