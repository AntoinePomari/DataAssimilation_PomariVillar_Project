%% enKF Twin Experiment Routine
%OBJ of this one: run an identical twin experiment, 
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

%% EnKF run for our model
N = 50; %size of the ensemble
[x,t0,s]=wave1d_initialize_enKF(s);
ksi = zeros(N,length(x));
ksi = ksi';
t=s.t;
times=s.times;
for i=1:length(t)
    for  jj=1:N 
        ksi(:,jj)=wave1d_timestep_enKF(ksi(:,jj),i,s);
    end
    x = (1/N)*sum(ksi,2);
    series_data(:,i)=x(ilocs);
end
%% load observations, plot ensemble run (no measurement update) VS observed data
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

wave1d_plotseries(times,series_data,s,observed_data)










