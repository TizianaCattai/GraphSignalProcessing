function [struct_TransfParam,struct_J]=Jdivergence_computation(z_0,z_1)

%input:
%z_0: matrix containing data under H0
%z_1: matrix containing data under H1

%output
%struct_TransfParam: structure with all the parameters of the transformation
%struct_J: structure with results associated to Jdiv

%% Transformation
[~,N]=size(z_0); %N: total number of components
eta_0=(mean(z_0,1))'; %mean value under Ho
eta_1=(mean(z_1,1))'; %mean value under H1
K_0=cov(z_0); %covariance matrix under Ho
K_1=cov(z_1); %covariance matrix under H1;
Q_0 = sqrtm(pinv(K_0,1e-5)); %1e-4: tolerance
B=ctranspose(Q_0)*K_1*Q_0;
[V_1,lamda] = eig(B);
T=real(ctranspose(V_1)*ctranspose(Q_0)); %Transformation matrix 
sigma2=lamda;

eta_n=T*(eta_1-eta_0);

%% J divergency computation
sigma=sqrt(sigma2);
sigma_real=real(diag(sigma));
theta=0.11; %threshold
J_div=zeros(N,1);

only_mean=(sigma_real<theta) | (abs(sigma_real-1)<theta);
also_var=not(only_mean);
eta=(real(eta_n));
sigma_real_inverse=(1./sigma_real); 
J_sigma=(sigma_real-sigma_real_inverse).^2+(eta.^2).*(sigma_real_inverse).*(sigma_real+sigma_real_inverse);
J_eta=4*eta.^2;


J_div(also_var)=J_sigma(also_var); % mean and var contributions
J_div(only_mean)=J_eta(only_mean); %only mean contrib

J=sum(J_div); %J divergency

%% output structure preparation

struct_TransfParam.T=T;
struct_TransfParam.eta_n=eta_n';
struct_TransfParam.eta_1=eta_1';
struct_TransfParam.eta_0=eta_0';
struct_TransfParam.sigma_real=sigma_real';
struct_TransfParam.sigma_real_inverse=sigma_real_inverse';
struct_TransfParam.only_mean=only_mean';
struct_TransfParam.also_var=also_var';
struct_TransfParam.theta=theta;

struct_J.J_div=J_div';
struct_J.J_eta=J_eta';
struct_J.J_sigma=J_sigma';
struct_J.J=J;
end