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
[obs_times,obs_values]=wave1d_read_series('tide_cadzand.txt');
observed_data=zeros(length(ilocs),length(obs_times));
observed_data(1,:)=obs_values(:);
[obs_times,obs_values]=wave1d_read_series('tide_vlissingen.txt');
observed_data(2,:)=obs_values(:);
[~,obs_values]=wave1d_read_series('tide_terneuzen.txt');
observed_data(3,:)=obs_values(:);
[~,obs_values]=wave1d_read_series('tide_hansweert.txt');
observed_data(4,:)=obs_values(:);
[~,obs_values]=wave1d_read_series('tide_bath.txt');
observed_data(5,:)=obs_values(:);

%% EnKF run, WITH measurement update included
N = 150; %size of the ensemble
[x,t0,s]=wave1d_initialize_enKF(s);
ksi = zeros(N,length(x));
ksi = ksi';
t=s.t;
times=s.times;
L=ksi;
H=zeros(5,201);
for ii = 1:length(xlocs_waterlevel)
    H(ii,ilocs(ii)) = 1; 
end
H = sparse(H);
for ii=1:length(t)
    %time step update
    for  jj=1:N 
        ksi(:,jj)=wave1d_timestep_enKF(ksi(:,jj),ii,s);
    end
    x = (1/N)*sum(ksi,2);
    series_data(:,ii)=x(ilocs);
    %measurement update: creating matrix L
    for kk =1:N
        L(:,kk) = ksi(:,kk)-x;
        L(:,kk) = (1/sqrt(N-1))*L(:,kk);
    end
    L=sparse(L);
    %creating PSI
    PSI = H*L;
    %size(PSI)
    K = (L*(PSI'))/(PSI*PSI'+speye(5)); 
%     v=zeros(5,1);
%     for dummy =1:5
%         v(dummy) = normrnd(0,0.2);
%     end
    for  jj=1:N 
        ksi(:,jj)=ksi(:,jj)+K*(observed_data(1:5,ii)-H*ksi(:,jj));  %+v we forget about measurement error
    end
end
%%
wave1d_plotseries(times,series_data,s,observed_data)







