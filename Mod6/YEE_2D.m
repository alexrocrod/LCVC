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
show=1;

% space gridding
dx=5.0e-3;
dy=dx;

% time step
m=1.0;
% used time step
dt=m/sqrt((1/dx)^2+(1/dy)^2)/c;

% constants
Cb=dt/dx/eps0;
Db=dt/dx/mu0;

% initialization
Hx=zeros(sizex,sizeyp1);
Hy=zeros(sizexp1,sizey);
Ez=zeros(sizex,sizey);

% simulation
for t=1:steps

   % excitation
   Ez(sizex,1:sizey)=ones(1,sizey)*sin(.3*t);

   % update magnetic field
   Hx(1:sizex,2:sizey)=Hx(1:sizex,2:sizey)+...
      Db*(Ez(1:sizex,1:sizeym1)-Ez(1:sizex,2:sizey));
   Hy(2:sizex,1:sizey)=Hy(2:sizex,1:sizey)+...
      Db*(Ez(2:sizex,1:sizey)-Ez(1:sizexm1,1:sizey));

   % update electric field
   Ez(1:sizex,1:sizey)=Ez(1:sizex,1:sizey)+...
      Cb*(Hy(2:sizexp1,1:sizey)-Hy(1:sizex,1:sizey)+...
      Hx(1:sizex,1:sizey)-Hx(1:sizex,2:sizeyp1));

   mesh(Ez);
   axis([0 50 0 50 -1.5 1.5]);
   xlabel('y'); ylabel('x');

   pause(0.1);

end
