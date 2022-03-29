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
%% load observations
% [obs_times,obs_values]=wave1d_read_series('tide_cadzand.txt');
% observed_data=zeros(length(ilocs),length(obs_times));
% observed_data(1,:)=obs_values(:);
[obs_times,obs_values]=wave1d_read_series('waterlevel_vlissingen.txt');
observed_data(2,:)=obs_values(:);
[~,obs_values]=wave1d_read_series('waterlevel_terneuzen.txt');
observed_data(3,:)=obs_values(:);
[~,obs_values]=wave1d_read_series('waterlevel_hansweert.txt');
observed_data(4,:)=obs_values(:);
[~,obs_values]=wave1d_read_series('waterlevel_bath.txt');
observed_data(5,:)=obs_values(:);
%% EnKF run, WITH measurement update
N = 400; %size of the ensemble
[x,t0,s]=wave1d_initialize_enKF(s);
% Initial guess for the ensemble members: NB it's a random guess basically,
% so anything goes
ksi = zeros(N,length(x));
%ksi = normrnd(0,0.2,N,length(x));
% ksi = eye(N,length(x));
ksi = ksi';
t=s.t;
times=s.times;
LL=ksi;
H=zeros(4,201);
H(1,51)=1;
H(2,101) = 1;
H(3,151) = 1;
H(4,199) = 1;
H = sparse(H);
deltat = 0;
for ii=1:165-deltat*6
    %time step update, with AR(1) noise at the left boundary
    for  jj=1:N 
        ksi(:,jj)=wave1d_timestep_enKF(ksi(:,jj),ii,s);
    end
    % forecast data
    x = (1/N)*sum(ksi,2);
    series_data(:,ii)=x(ilocs);
    
    %creating matrix L aka getting ready for the measurement update
    for kk =1:N
        LL(:,kk) = ksi(:,kk)-x;
        LL(:,kk) = (1/sqrt(N-1))*LL(:,kk);
    end
    LL=sparse(LL);
    
    
    %creating PSI for more efficient computations
    PSI = H*LL;    
    
    %Kalman gain: K = (LL*(PSI'))/(PSI*PSI'+R); %Which choice for matrix R?
      K = (LL*(PSI'))/(PSI*PSI'+0.01*speye(4)); %white-noise type of measurement err, STD err = 0.2, uncorrelated
%      K = (LL*(PSI'))/(PSI*PSI'); %considering perfect measurements

%measurement assimilation
    for  jj=1:N 
        ksi(:,jj)=ksi(:,jj)+K*(observed_data(2:5,ii+1)-H*ksi(:,jj)-normrnd(0,0.1,4,1));  %NB why+1? observed_data also contain initial time.
%         ksi(:,jj)=ksi(:,jj)+K*(observed_data(1:5,ii+1)-H*ksi(:,jj)); %line above in case of noisy observations, this line for perfect observations
    end
end

for ii=165-deltat*6+1:length(t)
    for  jj=1:N 
        ksi(:,jj)=wave1d_timestep_enKF(ksi(:,jj),ii,s);
    end
    x = (1/N)*sum(ksi,2);
    series_data(:,ii)=x(ilocs);
end
%% plotting observations vs Kalman predicted w/ measurement update
wave1d_plotseries(times,series_data,s,observed_data(1:5,:))
%save('forecast12.mat','series_data12')