clear all; clc;

% exact solution
exact = @(x,y) x*(1-x)*y*(1-y);
exact_x = @(x,y) (1-2*x)*y*(1-y);
exact_y = @(x,y) x*(1-x)*(1-2*y);
exact_xx = @(x,y) -2*y*(1-y);
exact_xy = @(x,y) (1-2*x)*(1-2*y);
exact_yy = @(x,y) -2*x*(1-x);

E = 1e9; % Young's Modulus
pr = 0.3; % Poison's ratio
lambda = pr*E/(1+pr)/(1-2*pr);
mu = E/2/(1+pr);
D = [lambda + 2*mu,lambda,0;lambda,lambda + 2*mu,0;0,0,mu];

f_11 = @(x,y) (lambda + 2*mu)*(-2*y*(1-y)) + lambda*(-2*x*(1-x));
f_22 = @(x,y) lambda*(-2*y*(1-y)) + (lambda + 2*mu)*(-2*x*(1-x));

% quadrature rule
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% mesh generation
n_en   = 4;               % number of nodes in an element
n_el_x = 30;               % number of elements in x-dir
n_el_y = 30;               % number of elements in y-dir
n_el   = n_el_x * n_el_y; % total number of elements

n_np_x = n_el_x + 1;      % number of nodal points in x-dir
n_np_y = n_el_y + 1;      % number of nodal points in y-dir
n_np   = n_np_x * n_np_y; % total number of nodal points

x_coor = zeros(n_np, 1);
y_coor = x_coor;

hx = 1.0 / n_el_x;        % mesh size in x-dir
hy = 1.0 / n_el_y;        % mesh size in y-dir

% generate the nodal coordinates
for ny = 1 : n_np_y
  for nx = 1 : n_np_x
    index = (ny-1)*n_np_x + nx; % nodal index
    x_coor(index) = (nx-1) * hx;
    y_coor(index) = (ny-1) * hy;
  end
end

% IEN array
IEN = zeros(n_el, n_en);
for ex = 1 : n_el_x
  for ey = 1 : n_el_y
    ee = (ey-1) * n_el_x + ex; % element index
    IEN(ee, 1) = (ey-1) * n_np_x + ex;
    IEN(ee, 2) = (ey-1) * n_np_x + ex + 1;
    IEN(ee, 3) =  ey    * n_np_x + ex + 1;
    IEN(ee, 4) =  ey    * n_np_x + ex;
  end
end

% ID array
ID = zeros(n_np,2); % 2 dof
counter = 0;
for ny = 2 : n_np_y - 1
  for nx = 2 : n_np_x - 1
    index = (ny-1)*n_np_x + nx;
    counter = counter + 1;
    ID(index,1) = counter;  
    counter = counter + 1;
    ID(index,2) = counter;
  end
end

n_eq = counter;

ID1 = ID(:,1);
ID2 = ID(:,2);
LM1 = ID1(IEN);
LM2 = ID2(IEN);
LM = zeros(n_el,n_en*2);
for aa = 1 : n_en
    LM(:,2*aa-1) = LM1(:,aa);
    LM(:,2*aa) = LM2(:,aa);
end

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
   x_ele = x_coor( IEN(ee, 1:n_en) );
  y_ele = y_coor( IEN(ee, 1:n_en) );

    k_ele = zeros(n_en*2, n_en*2); % element stiffness matrix
    f_ele_11 = zeros(n_en, 1);    % element load vector
    f_ele_22 = zeros(n_en, 1);
    f_ele = zeros(2*n_en,1);

    for ll = 1 : n_int
        x_l = 0.0; y_l = 0.0;
        dx_dxi = 0.0; dx_deta = 0.0;
        dy_dxi = 0.0; dy_deta = 0.0;
        for aa = 1 : n_en
            x_l = x_l + x_ele(aa) * Quad(aa, xi(ll), eta(ll));
            y_l = y_l + y_ele(aa) * Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            dx_dxi  = dx_dxi  + x_ele(aa) * Na_xi;
            dx_deta = dx_deta + x_ele(aa) * Na_eta;
            dy_dxi  = dy_dxi  + y_ele(aa) * Na_xi;
            dy_deta = dy_deta + y_ele(aa) * Na_eta;
        end

        detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

        for aa = 1 : n_en
            Na = Quad(aa, xi(ll), eta(ll));
            [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
            Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
            Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;

            f_ele_11(aa) = f_ele_11(aa) + weight(ll) * detJ * f_11(x_l,y_l) * Na;
            f_ele_22(aa) = f_ele_22(aa) + weight(ll) * detJ * f_22(x_l,y_l) * Na;
            f_ele(2*aa-1) = f_ele_11(aa);
            f_ele(2*aa) = f_ele_22(aa);

            Ba = zeros(3,2);
            Ba(1,2) = Na_x;
            Ba(2,2) = Na_y;
            Ba(3,1) = Na_y;
            Ba(3,2) = Na_x;
            for bb = 1:n_en
                Bb = zeros(3,2);
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                Bb(1,2) = Na_x;
                Bb(2,2) = Na_y;
                Bb(3,1) = Na_y;
                Bb(3,2) = Na_x;
                k_ij = Ba'*D*Bb;
                k_ele(2*aa-1:2*aa,2*bb-1:2*bb) = k_ele(2*aa-1:2*aa,2*bb-1:2*bb)+weight(ll).*detJ.*k_ij;
            end
        end 
    end

    for aa = 1 : n_en
        PP_1 = LM(ee, 2*aa-1);
        PP_2 = LM(ee, 2*aa);
        if PP_1 > 0
            F(PP_1) = F(PP_1) + f_ele(2*aa-1);

            for bb = 1 : n_en
                QQ_1 = LM(ee, 2*bb-1);
                QQ_2 = LM(ee, 2*bb);
                if QQ_1 > 0
                    K(PP_1, QQ_1) = K(PP_1, QQ_1) + k_ele(2*aa-1, 2*bb-1);
                else
                    % modify F with the boundary data 
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
                end
                if QQ_2 > 0
                    K(PP_1, QQ_2) = K(PP_1, QQ_2) + k_ele(2*aa-1, 2*bb);
                else
                    % modify F with the boundary data
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
                end
            end
        end
        if PP_2 > 0
            F(PP_2) = F(PP_2) + f_ele(2*aa);

            for bb = 1 : n_en
                QQ_1 = LM(ee, 2*bb-1);
                QQ_2 = LM(ee, 2*bb);
                if QQ_1 > 0
                    K(PP_2, QQ_1) = K(PP_2, QQ_1) + k_ele(2*aa, 2*bb-1);
                else
                    % modify F with the boundary data 
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
                end
                if QQ_2 > 0
                    K(PP_2, QQ_2) = K(PP_1, QQ_2) + k_ele(2*aa, 2*bb);
                else
                    % modify F with the boundary data
                    % here we do nothing because the boundary data g is zero or
                    % homogeneous
                end
            end
        end
    end

end

d = K\F;
d_x = zeros(length(d)/2,1);
d_y = zeros(length(d)/2,1);
distance = zeros(length(d)/2,1);
for i = 1:length(d)/2
    d_x(i) = d(2*i-1);
    d_y(i) = d(2*i);
    distance(i) = sqrt(d_x(i)^2+d_y(i)^2);
end
disp = zeros(n_np, 2);

for ii = 1 : n_np
  index_x = ID(ii,1);
  index_y = ID(ii,2);
  if index_x > 0 
    disp(ii,1) = d(index_x);
  else
    % modify disp with the g data. Here it does nothing because g is zero
  end
  if index_y > 0
      disp(ii,2) = d(index_y);
  else
      % modify disp with the g data. Here it does nothing because g is zero
  end
end

% plot(x_coor,y_coor,'.b');
% hold on
% plot((disp(:,1)+1).*x_coor,(disp(:,2)+1).*y_coor,'.r')