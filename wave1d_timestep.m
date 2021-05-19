function newx=wave1d_timestep(x,i,settings)
%function newx=timestep(x,i,settings)
%return (h,u) one timestep later
%take one timestep
    A=settings.A;
    B=settings.B;
    rhs=B*x;
    rhs(1)=settings.h_left(i); %left boundary, take at t+dt
    %rhs(end) = 0.2*sqrt(1-alpha^2)*normrand(0,1)
    newx=A\rhs;
