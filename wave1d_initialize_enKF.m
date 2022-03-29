function [x0,t0,settings_new]=wave1d_initialize_enKF(settings)
% function [x0,t0,settings_new]=wave1d_initialize_enKF(settings); 
%return the extended state at initial time, plus guess for one of the
%ensemble members. To be called for each ensemble member.
    settings_new=settings;
    %compute initial fields and cache some things for speed
    h_0=settings.h_0;
    u_0=settings.u_0;
    n=settings.n;
    x0=zeros(2*n+1,1); %order h[1],u[1],...h[n],u[n], N (our AR incertainty)
    x0(1:2:end-1)=u_0(:); 
    x0(2:2:end-1)=h_0(:);
    x0(end) = 0;
    %time
    t=settings.t;
    reftime=settings.reftime;
    dt=settings.dt;
    times=[];
    second=1./24./60./60.; %time unit is days in matlab
    times=reftime+t*second;
    settings_new.times=times;
    %initialize coefficients
    % create matrices in form A*x_new=B*x+alpha 
    % A and B are tri-diagonal sparse matrices 
    Adata=zeros(3,2*n+1); %order h[1],u[1],...h[n],u[n] N -->  3*2n+1
    Bdata=zeros(3,2*n+1); %3*2n +1
    %left boundary
    Adata(2,1)=1.0;
    %right boundary
    Adata(2,2*n)=1.0;
    % i=2,4,6,... du/dt  + g dh/sx + f u = 0
    % m=1/2,1 1/2, ...
    %  u(n+1,m) + 0.5 g dt/dx ( h(n+1,m+1/2) - h(n+1,m-1/2)) + 0.5 dt f u(n+1,m) 
    %= u(n  ,m) - 0.5 g dt/dx ( h(n  ,m+1/2) - h(n  ,m-1/2)) - 0.5 dt f u(n  ,m)
    g=settings.g;dx=settings.dx;f=settings.f;
    temp1=0.5*g*dt/dx;
    temp2=0.5*f*dt;
    for i=2:2:(2*n-2)
        Adata(1,i-1)= -temp1;
        Adata(2,i  )= 1.0 + temp2;
        Adata(3,i+1)= +temp1;
        Bdata(1,i-1)= +temp1;
        Bdata(2,i  )= 1.0 - temp2;
        Bdata(3,i+1)= -temp1;
    end
    % i=3,5,7,... dh/dt + D du/dx = 0
    % m=1,2,...
    %  h(n+1,m) + 0.5 D dt/dx ( u(n+1,m+1/2) - u(n+1,m-1/2))  
    %= h(n  ,m) - 0.5 D dt/dx ( u(n  ,m+1/2) - u(n  ,m-1/2))
    D=settings.D;
    temp1=0.5*D*dt/dx;
    for i=3:2:(2*n-1)
        Adata(1,i-1)= -temp1;
        Adata(2,i  )= 1.0;
        Adata(3,i+1)= +temp1;
        Bdata(1,i-1)= +temp1;
        Bdata(2,i  )= 1.0;
        Bdata(3,i+1)= -temp1;
    end
    Adata(2,end) = 1;
    Bdata(2,end) = settings.alpha;
    %ALSO: CREATE VECTOR G (zeros everywhere, 0.2*sqrt(1-alpha^2) at the
    %end (do we really need vector G??)
    % build sparse matrix
    A=spdiags(Adata',[-1,0,1],2*n+1,2*n+1);
    B=spdiags(Bdata',[-1,0,1],2*n+1,2*n+1);
    B(:,size(B,1)) = 1;
    size(B)
    size(A)
    settings_new.A=A; %cache for later use
    settings_new.B=B;
    t0=t(1);