%% Part 2 of Twin experiments: load the truth and observations, run for ensemble size 150
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
load('z_twin3.mat'), %load('v_noisy.mat'); %load our 'truth' and the perfect observations based on our truth

N = 600; %size of the ensemble: for N>230 results are pretty good in fact
[x,t0,s]=wave1d_initialize_enKF(s);
%Initial guess for the ensemble: what to choose? Anything goes basically
%OR NOT?
%N
ksi = normrnd(0,0.2,N,length(x));
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
%     K = (LL*(PSI'))/(PSI*PSI'+0.2*speye(5)); %white-noise type of measurement err, STD err = 0.2, uncorrelated
%   error
    K = (LL*(PSI'))/(PSI*PSI');  %perfect observations
    for  jj=1:N
        ksi(:,jj)=ksi(:,jj)+K*(z_twin3(1:5,ii)-H*ksi(:,jj)); %if we use perfect observations NB why not ii+1 for the observations (like wave1d?): because we take it into account when we generate the truth and our true data
%         ksi(:,jj)=ksi(:,jj)+K*(z_twin(1:5,ii)-H*ksi(:,jj)-vector_of_noises(:,ii)); %if white noise is present
    end
end
%%
% perfect_data=zeros(5,289);
% for ii=1:289
%     %evolving in time, our new truth
%     x_twin=truth_of_model(:,ii);
%     %series_twin(:,ii)=x_twin(ilocs);
% %    v = vector_ %needed to store all the white-noise error, otherwise there is a mismatch during run of the twin model
%     perfect_data(:,ii+1) = H*x_twin; %if we consider perfect obs
% %   z_twin(:,ii+1) = H*x_twin+v; %if we assume white noise
% end


wave1d_plotseries(times,series_data,s,z_twin3)




