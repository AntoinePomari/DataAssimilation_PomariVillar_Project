%% Checking the convergence rate
clear all
clc

s=wave1d_settings();
L=s.L;
dx=s.dx;
xlocs_waterlevel=[0.0*L,0.25*L,0.5*L,0.75*L,0.99*L];
xlocs_velocity=[0.0*L,0.25*L,0.5*L,0.75*L];
ilocs=[round(xlocs_waterlevel/dx)*2+1,round(xlocs_velocity/dx)*2+2]; %indices of waterlevel locations in x
loc_names={};
names={'Cadzand','Vlissingen','Terneuzen','Hansweert','Bath'};
for i=1:length(xlocs_waterlevel)
    loc_names{i}=sprintf('Waterlevel at x=%f km %s',0.001*xlocs_waterlevel(i),names{i});
end
nn=length(loc_names);
for i=(1+nn):(length(xlocs_velocity)+nn)
    loc_names{i}=sprintf('Velocity at x=%f km %s',0.001*xlocs_velocity(i-nn),names{i-nn});
end
s.xlocs_waterlevel=xlocs_waterlevel;
s.xlocs_velocity=xlocs_velocity;
s.ilocs=ilocs;
s.loc_names=loc_names;
%% EnKF run, FOR the TWIN experiment
load('x_twin2.mat'), load('z_twin3_noisy.mat'), load('v3_noisy.mat'); %load our 'truth' and the perfect observations based on our truth
%rng(1)
results=[];
NN = [50, 100, 200, 500, 1000, 1500, 3000];
for ll=1:length(NN)
N = NN(ll); %size of the ensemble: for N>200 results are decent in fact
[x,t0,s]=wave1d_initialize_enKF(s);

%Initial guess for the ensemble: what to choose? Anything goes basically
ksi = normrnd(0,0.1,N,length(x));
%ksi = zeros(N,length(x));
ksi = ksi';
t=s.t;
times=s.times;
LL=ksi;
H=zeros(5,201);
for ii = 1:length(xlocs_waterlevel)
    H(ii,ilocs(ii)) = 1; 
end
H = sparse(H);
for ii=1:length(t)
    %time step update with AR(1) noise
    for  jj=1:N 
        ksi(:,jj)=wave1d_timestep_enKF(ksi(:,jj),ii,s);
    end
    x = (1/N)*sum(ksi,2);
    series_data(:,ii)=x(ilocs);
    %Getting ready for the measurement update: we create matrix L 
    for kk =1:N
        LL(:,kk) = ksi(:,kk)-x;
        LL(:,kk) = (1/sqrt(N-1))*LL(:,kk);
    end
    LL=sparse(LL);
    %creating PSI for ease of computation
    PSI = H*LL;
    %size(PSI)
    
    %Which matrix R to choose?
      K = (LL*(PSI'))/(PSI*PSI'+0.01*speye(5)); 
%     K = (LL*(PSI'))/(PSI*PSI'); 
    for  jj=1:N
        ksi(:,jj)=ksi(:,jj)+K*(z_twin(1:5,ii)-H*ksi(:,jj)-vector_of_noises(:,ii)); 
    end
end
results((1+(ll-1)*5):5*ll,:)=series_data(1:5,:);
MSE6(ll)=sqrt(sum(sum((results((1+(ll-1)*5):5*ll,:)-z_twin(:,2:end)).^2))/1440);
end
semilogx(NN,MSE6)







