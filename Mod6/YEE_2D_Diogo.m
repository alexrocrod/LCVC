% 2D Yee algorithm: TE wave propagating in xy plane
% free space

clear; close all;
hold off;

% constitutive parameters
eps0=8.8541878e-12;
mu0=1.2566371e-6;
c=1/sqrt(eps0*mu0);
eta0=sqrt(mu0/eps0);

% grid size
sizex=101;
sizey=101;

% size plus-minus 1
sizexp1=sizex+1;
sizeyp1=sizey+1;
sizexm1=sizex-1;
sizeym1=sizey-1;

% number of time steps
steps=1000;

% show starts at
% show=1;
% task 2
show=30; 

% space gridding
dx=5.0e-3;
dy=dx;

% time step
% task 4
%m=1.42;
%m=1.41;
% task 5
m=0.9;
% task 8
%m=1.01;
% used time step
dt=m/sqrt((1/dx)^2+(1/dy)^2)/c;

% constants
Cb=dt/dx/eps0;
Db=dt/dx/mu0;

% initialization
Hx=zeros(sizex,sizeyp1);
Hy=zeros(sizexp1,sizey);
Ez=zeros(sizex,sizey);

% task 6
x = floor((1 + sizexp1) / 2);
% task 7
y_min_delta = floor(1 + 35);
y_max_delta = floor(sizey - 35);


% simulation
for t=1:steps
    
   % task 5
   % Eletric field boundary conditions
   %Ez(1, 1:sizey) = eta0 * Hy(1, 1:sizey);
   %Ez(sizex, 1:sizey) = - eta0 * Hy(sizexp1, 1:sizey);
   % Magnetic field boundary conditions
   Hy(1,1:sizey)= (1/eta0) * Ez(1,1:sizey);
   Hy(sizexp1,1:sizey)= -(1/eta0) * Ez(sizex,1:sizey);

%    Hy(x, y_min_delta:y_max_delta) = ones(1, y_max_delta - y_min_delta + 1)*sin(.3*t)/eta0; 
   
   
   % update magnetic field
   Hx(1:sizex,2:sizey)=Hx(1:sizex,2:sizey)+...
      Db*(Ez(1:sizex,1:sizeym1)-Ez(1:sizex,2:sizey));
   Hy(2:sizex,1:sizey)=Hy(2:sizex,1:sizey)+...
      Db*(Ez(2:sizex,1:sizey)-Ez(1:sizexm1,1:sizey));

   % update electric field
   Ez(1:sizex,1:sizey)=Ez(1:sizex,1:sizey)+...
      Cb*(Hy(2:sizexp1,1:sizey)-Hy(1:sizex,1:sizey)+...
      Hx(1:sizex,1:sizey)-Hx(1:sizex,2:sizeyp1));
  
   % excitation default
   %Ez(sizex,1:sizey)=ones(1,sizey)*sin(.3*t);
   % task 3 excitation y component of magnetic field
   %Hy(1, 1:sizey)=ones(1, sizey)*sin(.3*t)/eta0; 
   % task 6 excitation in the middle of the magnetic field domain
   %Hy(x, 1:sizey)=ones(1, sizey)*sin(.3*t)/eta0; 
   % task 7 restricting the excitation
   Hy(x, y_min_delta:y_max_delta) = ones(1, y_max_delta - y_min_delta + 1)*sin(.3*t)/eta0; 
    
   % task 2 start from show time
   if t >= show
       mesh(Ez);
       % task 1 to see the complete domain
       axis([0 sizex 0 sizey -1.5 1.5]);
       xlabel('y'); ylabel('x');
       pause(0.1);
   end
   
end
