clear all; clc; clf;close all; % clean the memory, screen, and figure

% Problem definition
f = @(x) -20*x.^3; % f(x) is the source
g = 1.0;           % u    = g  at x = 1
h = 0.0;           % -u,x = h  at x = 0
lnh_table = zeros(1,8);
lnEL2_table = zeros(1,8);
lnEH1_table = zeros(1,8);
for n_el = 2:2:16      % number of elements

% Setup the mesh
pp = 2;                % polynomial degree
n_en = pp + 1;         % number of element or local nodes

n_np = n_el * pp + 1;  % number of nodal points
n_eq = n_np - 1;       % number of equations
n_int = 20;

hh = 1.0 / (n_np - 1); % space between two adjacent nodes
x_coor = 0 : hh : 1;   % nodal coordinates for equally spaced nodes

IEN = zeros(n_el, n_en);

for ee = 1 : n_el
  for aa = 1 : n_en
    IEN(ee, aa) = (ee - 1) * pp + aa;
  end
end

% Setup the ID array for the problem
ID = 1 : n_np;
ID(end) = 0;

% Setup the quadrature rule
[xi, weight] = Gauss(n_int, -1, 1);

% allocate the stiffness matrix
K = spalloc(n_eq, n_eq, (2*pp+1)*n_eq);
F = zeros(n_eq, 1);

% Assembly of the stiffness matrix and load vector
for ee = 1 : n_el
  k_ele = zeros(n_en, n_en); % allocate a zero element stiffness matrix
  f_ele = zeros(n_en, 1);    % allocate a zero element load vector

  x_ele = x_coor(IEN(ee,:)); % x_ele(aa) = x_coor(A) with A = IEN(aa, ee)

  % quadrature loop
  for qua = 1 : n_int    
    dx_dxi = 0.0;
    x_l = 0.0;
    for aa = 1 : n_en
      x_l    = x_l    + x_ele(aa) * PolyShape(pp, aa, xi(qua), 0);
      dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(qua), 1);
    end
    dxi_dx = 1.0 / dx_dxi;

    for aa = 1 : n_en
      f_ele(aa) = f_ele(aa) + weight(qua) * PolyShape(pp, aa, xi(qua), 0) * f(x_l) * dx_dxi;
      for bb = 1 : n_en
        k_ele(aa, bb) = k_ele(aa, bb) + weight(qua) * PolyShape(pp, aa, xi(qua), 1) * PolyShape(pp, bb, xi(qua), 1) * dxi_dx;
      end
    end
  end
 
  % Assembly of the matrix and vector based on the ID or LM data
  for aa = 1 : n_en
    P = ID(IEN(ee,aa));
    if(P > 0)
      F(P) = F(P) + f_ele(aa);
      for bb = 1 : n_en
        Q = ID(IEN(ee,bb));
        if(Q > 0)
          K(P, Q) = K(P, Q) + k_ele(aa, bb);
        else
          F(P) = F(P) - k_ele(aa, bb) * g; % handles the Dirichlet boundary data
        end
      end
    end
  end
end

% ee = 1 F = NA(0)xh
F(ID(IEN(1,1))) = F(ID(IEN(1,1))) + h;

% Solve Kd = F equation
d_temp = K \ F;

disp = [d_temp; g];

% Postprocessing: visualization
%plot(x_coor, disp, '--r','LineWidth',3);

%x_sam = 0 : 0.01 : 1;
%y_sam = x_sam.^5;
%hold on;
%plot(x_sam, y_sam, '-k', 'LineWidth', 3);
%   x_ele = x_coor( IEN(ee, :) );
%   u_ele = disp( IEN(ee, :) );
% 
%   if ee == n_el
%     n_sam_end = n_sam+1;
%   else
%     n_sam_end = n_sam;
%   end
% 
%   for ll = 1 : n_sam_end
%     x_l = 0.0;
%     u_l = 0.0;
%     for aa = 1 : n_en
%       x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
%       u_l = u_l + u_ele(aa) * PolyShape(pp, aa, xi_sam(ll), 0);
%     end
% 
%     x_sam( (ee-1)*n_sam + ll ) = x_l;
%     u_sam( (ee-1)*n_sam + ll ) = u_l;
%     y_sam( (ee-1)*n_sam + ll ) = x_l^5;
%   end
% end
% 
% plot(x_sam, u_sam, '-r','LineWidth',3);
% hold on;
% plot(x_sam, y_sam, '-k','LineWidth',3);

% Calculate error L2 and H1
L2 = 0;
H1 = 0;
nqp = 10;
[xi, weight] = Gauss(nqp,-1,1);
L2_top = 0.0; L2_bot = 0.0; H1_top = 0.0; H1_bot = 0.0;

for ee = 1 : n_el
    x_ele = x_coor( IEN(ee, :) );
    u_ele = disp( IEN(ee, :) );

    for ll = 1 : nqp
        x_l = 0.0;
        uh = 0.0;
        uh_xi = 0.0;
        dx_dxi = 0.0;

        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
            uh = uh + u_ele(aa) * PolyShape(pp, aa, xi(ll), 0);
            dx_dxi = dx_dxi + x_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
            uh_xi = uh_xi + u_ele(aa) * PolyShape(pp, aa, xi(ll), 1);
        end
        dxi_dx = 1/dx_dxi;

        L2_top = L2_top + weight(ll) * (uh - x_l^5)^2 * dx_dxi;
        L2_bot = L2_bot + weight(ll) * (x_l^5)^2 * dx_dxi;

        H1_top = H1_top + weight(ll) * ( uh_xi * dxi_dx - 5 * x_l^4 )^2 * dx_dxi;
        H1_bot = H1_bot + weight(ll) * (5 * x_l^4)^2 * dx_dxi;

    end
end

L2_top = sqrt(L2_top); L2_bot = sqrt(L2_bot);
H1_top = sqrt(H1_top); H1_bot = sqrt(H1_bot);

L2_error = L2_top / L2_bot;
H1_error = H1_top / H1_bot;

lnEL2 = log(L2_error);
lnEH1 = log(H1_error);
lnh = log(hh);
lnEL2_table(n_el/2) = lnEL2;
lnEH1_table(n_el/2) = lnEH1;
lnh_table(n_el/2) = lnh;

end
f1 = figure;
plot(lnh_table,lnEL2_table);
hold on
plot(lnh_table,lnEH1_table);
legend('e_L_2','e_H_1','Location','southeast')
xlabel('ln(h)')
ylabel('ln(e)')
title('log-log plot of Error versus mesh size')
p_L2 = polyfit(lnh_table,lnEL2_table,1);
p_H1 = polyfit(lnh_table,lnEH1_table,1);
slope_eL2 = p_L2(1)
slope_eH1 = p_H1(1)





























% EOF