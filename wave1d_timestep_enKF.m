function newx=wave1d_timestep_enKF(x,i,settings)
%function newx=timestep(x,i,settings)
%return (h,u, N) one timestep later: NB this is for the Kalman ensemble
%Filter
%take one timestep
    A=settings.A;
    B=settings.B;
    rhs=B*x+[zeros(200,1);normrnd(0,0.2*sqrt(1-settings.alpha^2))];
%     rhs(1)=settings.h_left(i); %left boundary, take at t+dt
%     rhs(end) = 0.2*sqrt(1-settings.alpha^2)*normrnd(0,1); %endpoint: autoregressive update
    rhs(1)=settings.h_left(i)+x(end); %left boundary, take at t+dt
    newx=A\rhs;