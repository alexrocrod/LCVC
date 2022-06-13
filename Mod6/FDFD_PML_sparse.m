% Two-dimensional FDFD: TE wave propagating in the xy plane.
% Stanislav Maslovski, DETI, University of Aveiro, May 2020.

clear; close all;

% constitutive parameters
eps0 = 8.8541878e-12;
mu0 = 1.2566371e-6;
c = 1/sqrt(eps0*mu0);
eta0 = sqrt(mu0/eps0);

% grid size
sizex = 141;
sizey = 241;

% size plus-minus 1
sizexp1 = sizex+1;
sizeyp1 = sizey+1;
sizexm1 = sizex-1;
sizeym1 = sizey-1;

% the frequency
freq = 10e9;
omega = 2*pi*freq;

% the free space wavelength
lambda = c/freq;

% space gridding
dx = lambda/10;
dy = dx;

% grid of E-points
[egrid_y,egrid_x] = meshgrid((0:sizeym1)*dy,(0:sizexm1)*dx);

% grid of Hx points
[hxgrid_y,hxgrid_x] = meshgrid((0:sizey)*dy-0.5*dy,(0:sizexm1)*dx);

% grid of Hy points
[hygrid_y,hygrid_x] = meshgrid((0:sizeym1)*dy,(0:sizex)*dx-0.5*dx);

% material block function
xmin = (sizexm1-40)*dx/4 + 20*dx;
xmax = 3*(sizexm1-40)*dx/4 + 20*dx;
% when using PML in y, decrease the block size so that there is no overlap!
ymin = 0*dy;
ymax = sizeym1*dy;

fun = @(x,y) x > xmin & x < xmax & y > ymin & y < ymax;

% zz component of the electric susceptibility tensor
% (values are relative to eps0)
%chi_ee_zz = zeros(sizex,sizey);
chi_ee_zz = (-2.0-1e-4i)*fun(egrid_x,egrid_y);

% averaging on x cell boundaries
%chi_ee_zz(1:sizexm1,1:sizey) = 0.5*(chi_ee_zz(1:sizexm1,1:sizey)+chi_ee_zz(2:sizex,1:sizey));
% averaging on y cell boundaries
%chi_ee_zz(1:sizex,1:sizeym1) = 0.5*(chi_ee_zz(1:sizex,1:sizeym1)+chi_ee_zz(1:sizex,2:sizey));

% xx and yy components of the magnetic susceptibility tensor
% (values are relative to mu0)
%chi_mm_xx = zeros(sizex,sizeyp1);
%chi_mm_yy = zeros(sizexp1,sizey);
chi_mm_xx = (-2.0-1e-4i)*fun(hxgrid_x,hxgrid_y);
chi_mm_yy = (-2.0-1e-4i)*fun(hygrid_x,hygrid_y);


%%%%%%%%%%%%%%%%%%%%%%% PML %%%%%%%%%%%%%%%%%%%%%%%%%

% how many cells across PML?
ncells = 15;

% PML block function at x = x_min
pml1 = @(x,y) x < ncells*dx;
% PML block function at x = x_max
pml2 = @(x,y) x > sizexm1*dx - ncells*dx;
% PML block function at y = y_min
pml3 = @(x,y) y < ncells*dy;
% PML block function at y = y_max
pml4 = @(x,y) y > sizeym1*dy - ncells*dy;

% PML scale functions
Lpml = (ncells-1)*dx;
smax = 5;
p = 4;
pf = @(x) (1-1j*sin(pi*x/2/Lpml).^2).*(1+smax*(x/Lpml).^p) - 1;
sx = @(x,y) 1 + pml1(x,y).*pf(Lpml-x) + pml2(x,y).*pf(x-sizexm1*dx+Lpml);
%sy = @(x,y) 1 + pml3(x,y).*pf(Lpml-y) + pml4(x,y).*pf(y-sizeym1*dy+Lpml);
sy = @(x,y) 1;

%%%%%%%%%%%%%%% effective material parameters %%%%%%%%%%%%%%%%

epsr_zz = (1+chi_ee_zz).*sx(egrid_x,egrid_y).*sy(egrid_x,egrid_y);
mur_xx = (1+chi_mm_xx).*sy(hxgrid_x,hxgrid_y)./sx(hxgrid_x,hxgrid_y);
mur_yy = (1+chi_mm_yy).*sx(hygrid_x,hygrid_y)./sy(hygrid_x,hygrid_y);

%%%%%%%%%%%%%%%%%%%  Yee scheme factors %%%%%%%%%%%%%%%%%%%%%

dt = 1/(1i*omega);
inv_eps = dt/dx./(eps0*epsr_zz);
inv_mux = dt/dx./(mu0*mur_xx);
inv_muy = dt/dx./(mu0*mur_yy);

% form system matrix for Ez (Hx and Hy will be eliminated)

sys_matr = sparse(sizex*sizey,sizex*sizey);

jj = 1;
for iy=1:sizey
  for ix=1:sizex

   % initial zero fields
   Hx = sparse(sizex,sizeyp1);
   Hy = sparse(sizexp1,sizey);
   Ez = sparse(sizex,sizey);

   % define just one component of Ez
   Ez(ix,iy) = 1;

   % express magnetic fields through Ez
   Hx(1:sizex,2:sizey) =...
      inv_mux(1:sizex,2:sizey).*(Ez(1:sizex,1:sizeym1)-Ez(1:sizex,2:sizey));
   Hy(2:sizex,1:sizey) =...
      inv_muy(2:sizex,1:sizey).*(Ez(2:sizex,1:sizey)-Ez(1:sizexm1,1:sizey));

   % Periodic BCs at y=y_min and y=y_max
   Hx(1:sizex,1) =...
      inv_mux(1:sizex,1).*(Ez(1:sizex,sizey)-Ez(1:sizex,1));
   % Note that inv_mux(1:sizex,1) = inv_mux(1:sizex,sizeyp1) in a truly periodic domain.
   % In the case of discontinuity in mu_x at the periodic boundary, mu_x must be averaged!
   Hx(1:sizex,sizeyp1) =...
      inv_mux(1:sizex,sizeyp1).*(Ez(1:sizex,sizey)-Ez(1:sizex,1));

   % form electric field equation
   Ez(1:sizex,1:sizey) = Ez(1:sizex,1:sizey)-...
      inv_eps(1:sizex,1:sizey).*(Hy(2:sizexp1,1:sizey)-Hy(1:sizex,1:sizey)+...
      Hx(1:sizex,1:sizey)-Hx(1:sizex,2:sizeyp1));

   % store result in a column of the system matrix
   sys_matr(:,jj) = reshape(Ez,sizex*sizey,1);

   % update column index
   jj = jj + 1;

  end
end

% check the structure of the system matrix
figure;
spy(sys_matr);

% set up excitation
Jz = sparse(sizex,sizey);
Jz(20,round(sizeyp1/2)) = -eta0;
%Jz(20,round(sizeyp1/2)+10) = -eta0;
%Jz(20,round(sizeyp1/2)-10) = -eta0;

% define the right-hand side
rhs_vec = reshape(-dx*inv_eps(1:sizex,1:sizey).*Jz,sizex*sizey,1);

% solve the system
sol_vec = sys_matr\rhs_vec;

% get the fields
Ez = reshape(sol_vec,sizex,sizey);
% express magnetic fields through Ez
Hx(1:sizex,2:sizey) =...
   inv_mux(1:sizex,2:sizey).*(Ez(1:sizex,1:sizeym1)-Ez(1:sizex,2:sizey));
Hy(2:sizex,1:sizey) =...
   inv_muy(2:sizex,1:sizey).*(Ez(2:sizex,1:sizey)-Ez(1:sizexm1,1:sizey));
% Periodic BCs at y=y_min and y=y_max
Hx(1:sizex,1) =...
   inv_mux(1:sizex,1).*(Ez(1:sizex,sizey)-Ez(1:sizex,1));
Hx(1:sizex,sizeyp1) =...
   inv_mux(1:sizex,sizeyp1).*(Ez(1:sizex,sizey)-Ez(1:sizex,1));

% plot real part of Ez
figure;
reEz = real(Ez);
surf(egrid_x,egrid_y,reEz,'linestyle','none');
xlabel('x'); ylabel('y');
grid off;
%axis('square');
daspect([1,1,max(max(reEz))/min(sizexm1*dx,sizeym1*dy)]);
view(2);
colorbar;

% make a vector plot of H = (Hx,Hy)
% before plotting, interpolate the magnetic fields at the E-points
figure;
Hx = real(0.5*(Hx(1:sizex,1:sizey)+Hx(1:sizex,2:sizeyp1)));
Hy = real(0.5*(Hy(1:sizex,1:sizey)+Hy(2:sizexp1,1:sizey)));
magH = sqrt(Hx.^2+Hy.^2);
step = 3;
scale = 0.7;
threshold = 0.05*max(max(magH));
quiver(egrid_x(1:step:sizex,1:step:sizey),egrid_y(1:step:sizex,1:step:sizey),...
   Hx(1:step:sizex,1:step:sizey)./(threshold+magH(1:step:sizex,1:step:sizey)),...
   Hy(1:step:sizex,1:step:sizey)./(threshold+magH(1:step:sizex,1:step:sizey)),...
   scale,'LineWidth',1,'MaxHeadSize',0.3);
xlabel('x'); ylabel('y');
%axis('square');
% daspect([sizexm1*dx,sizeym1*dy]);
daspect([sizexm1*dx,sizeym1*dy,0]);
