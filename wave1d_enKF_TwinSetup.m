%% enKF Twin Experiment Routine
%OBJ of this one: run an identical twin experiment, and use the stored data
%to compute later in the enKF
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
%% TWIN: 
% run the stochastic model FW to get x "truth", z_obs
%to be used later in the enKF
[x_twin3,t0,s]=wave1d_initialize_enKF(s); %identical twin same stochastic numerical model as for enKF
t=s.t;
%system matrix: will be needed for the twin observations
H=zeros(5,201);
for ii = 1:length(xlocs_waterlevel)
    H(ii,ilocs(ii)) = 1; 
end
H = sparse(H);
z_twin3=zeros(5,length(t));
%z_twin2(:,1)=H*x_twin2;
for ii=1:length(t)
   %evolving in time, our new truth
   x_twin3=wave1d_timestep_enKF(x_twin3,ii,s);
   %series_twin(:,ii)=x_twin(ilocs);
   %basing our new observation on the 'truth', considering PERFECT
   %OBSERVATIONS
   z_twin3(:,ii+1) = H*x_twin3;
end
save('x_twin3.mat','x_twin3'), save('z_twin3.mat','z_twin3');

