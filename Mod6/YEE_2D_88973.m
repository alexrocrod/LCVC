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
show=30;

% space gridding
dx=5.0e-3;
dy=dx;

% time step
m=1.0;

% question 4
% m = 1.42;
% m = 1.41;

% %question 5
% m = 0.9;
% % question 8
% m = 1.01;

% used time step
dt=m/sqrt((1/dx)^2+(1/dy)^2)/c;

% constants
Cb=dt/dx/eps0;
Db=dt/dx/mu0;

% initialization
Hx=zeros(sizex,sizeyp1);
Hy=zeros(sizexp1,sizey);
Ez=zeros(sizex,sizey);

boundary_condition=1; %bondary_condition=true; 0=false

% simulation
for t=1:steps
    
%    % impedance boundary conditions - question 5
%    Hy(1,1:sizey)=Ez(1,1:sizey)/eta0;
%    Hy(sizexp1,1:sizey)=-Ez(sizex,1:sizey)/eta0;

%    % excitation - question 3
%    Hy(1,1:sizey)=ones(1,sizey)*sin(.3*t)/eta0;
   
%    % question 6
%    Hy(51,1:sizey) = (1/eta0)*ones(1,sizey)*sin(.3*t); 
   
%    % question 7 and 8
%    Hy(sizexp1/2,35:1:sizey-35)=ones(1,length(35:1:sizey-35))*sin(.3*t)/eta0; % ymin+35*dy<y<ymax-35dy
   
   % update magnetic field
   Hx(1:sizex,2:sizey)=Hx(1:sizex,2:sizey)+...
       Db*(Ez(1:sizex,1:sizeym1)-Ez(1:sizex,2:sizey));
   Hy(2:sizex,1:sizey)=Hy(2:sizex,1:sizey)+...
       Db*(Ez(2:sizex,1:sizey)-Ez(1:sizexm1,1:sizey));
   
   % update electric field
   Ez(1:sizex,1:sizey)=Ez(1:sizex,1:sizey)+...
      Cb*(Hy(2:sizexp1,1:sizey)-Hy(1:sizex,1:sizey)+...
      Hx(1:sizex,1:sizey)-Hx(1:sizex,2:sizeyp1));
  
  % question 2
  titulo = sprintf('time = %f s', t);
  if t >= show
   mesh(Ez);
   colorbar;
   % axis([0 s 0 50 -1.5 1.5]);
   axis([0 sizex 0 sizey -1.5 1.5]); % question 1
   xlabel('y'); ylabel('x');
   title(titulo)

   pause(0.1);
  end

end
