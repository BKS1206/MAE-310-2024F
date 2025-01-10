clear;clc; close all;


%% Read nodes and IEN for case 1
    % 打开文件
    fid = fopen('quarter-plate-with-hole-quad-Refine3.msh', 'r');

    if fid == -1
        error('无法打开文件');
    end
    
    % 初始化
    i = 1;
    j = 1;
    coor = zeros(1,3);
    ien  = zeros(1,5);
    % 读取文件内容
   while ~feof(fid)
        line = fgetl(fid);
        num  = strsplit(line);
        num = str2double(num(1,:));
        Y = isnan(num);
        num = num (~Y);
        isdata = all(isnan(num));
        if isdata == 0
            if length(num) == 3
                coor(i,:) = num;
                i = i+1;
            end
            if length(num) == 5
                ien(j,:) = num;
                j = j+1;
            end
        end
   end
    node_coor_1 = coor(2:563,:);
    node_coor_1(:,3) = [];
    node_coor_1(:,1) = node_coor_1(:,1) + 1;
    node_coor_1(:,2) = node_coor_1(:,2) + 1;
    IEN_1 = ien(8:519,2:5);
    % 关闭文件
    fclose(fid);

%% Read nodes and IEN for case 2
   %  % 打开文件
   %  fid = fopen('quarter-plate-with-hole-quad-Refine4.msh', 'r');
   % 
   %  if fid == -1
   %      error('无法打开文件');
   %  end
   % 
   %  % 初始化
   %  i = 1;
   %  j = 1;
   %  coor = zeros(1,3);
   %  ien  = zeros(1,5);
   %  % 读取文件内容
   % while ~feof(fid)
   %      line = fgetl(fid);
   %      num  = strsplit(line);
   %      num = str2double(num(1,:));
   %      Y = isnan(num);
   %      num = num (~Y);
   %      isdata = all(isnan(num));
   %      if isdata == 0
   %          if length(num) == 3
   %              coor(i,:) = num;
   %              i = i+1;
   %          end
   %          if length(num) == 5
   %              ien(j,:) = num;
   %              j = j+1;
   %          end
   %      end
   % end
   %  node_coor_2 = coor(2:2147,:);
   %  node_coor_2(:,3) = [];
   %  IEN_2 = ien(8:2055,2:5);
   %  % 关闭文件
   %  fclose(fid);
%% Read nodes and IEN for case 3
   %  % 打开文件
   %  fid = fopen('quarter-plate-with-hole-quad-Refine5.msh', 'r');
   % 
   %  if fid == -1
   %      error('无法打开文件');
   %  end
   % 
   %  % 初始化
   %  i = 1;
   %  j = 1;
   %  coor = zeros(1,3);
   %  ien  = zeros(1,5);
   %  % 读取文件内容
   % while ~feof(fid)
   %      line = fgetl(fid);
   %      num  = strsplit(line);
   %      num = str2double(num(1,:));
   %      Y = isnan(num);
   %      num = num (~Y);
   %      isdata = all(isnan(num));
   %      if isdata == 0
   %          if length(num) == 3
   %              coor(i,:) = num;
   %              i = i+1;
   %          end
   %          if length(num) == 5
   %              ien(j,:) = num;
   %              j = j+1;
   %          end
   %      end
   % end
   %  node_coor_3 = coor(2:8387,:);
   %  node_coor_3(:,3) = [];
   %  IEN_3 = ien(8:8199,2:5);
   %  % 关闭文件
   %  fclose(fid);
%% Manufactured Solution 

sigma_rr = @(Tx,R,r,theta) Tx./2.*(1 - R.^2./r.^2) + Tx./2.*(1 - 4.*R.^2./r.^2 + 3.*R.^4./r.^4).*cos(2.*theta);
sigma_tt = @(Tx,R,r,theta) Tx./2.*(1 + R.^2./r.^2) - Tx./2.*(1 + 3.*R.^4./r.^4).*cos(2.*theta);
sigma_rt = @(Tx,R,r,theta) -Tx./2.*(1 + 2.*R.^2./r.^2 - 3.*R.^4./r.^4).*sin(2.*theta);
% Find the outer boudary points
r_right = find(node_coor_1(:,1) == 1);
r_top   = find(node_coor_1(:,2) == 1);
r_right = r_right(3:length(r_right));
r_top   = r_top(3:length(r_top));
right_coor = node_coor_1(r_right,:);
top_coor   = node_coor_1(r_top,:);
Tx = 1e4; % Boundary condition 
R = 0.3;
r = sqrt((node_coor_1(:,1)-(-1)).^2 + (node_coor_1(:,2)-(-1)).^2);
theta = atan((node_coor_1(:,2)-(-1))./((node_coor_1(:,1)-(-1))));
sigma_rr = sigma_rr(Tx,R,r,theta);
sigma_tt = sigma_tt(Tx,R,r,theta);
sigma_rt = sigma_rt(Tx,R,r,theta);
% f
sigma_11 = (sigma_rr + sigma_tt)./2 + (sigma_rr - sigma_tt)./2.*cos(2.*theta) - sigma_rt.*sin(2.*theta);
sigma_22 = (sigma_rr + sigma_tt)./2 - (sigma_rr - sigma_tt)./2.*cos(2.*theta) + sigma_rt.*sin(2.*theta);
sigma_12 = (sigma_rr - sigma_tt)./2.*sin(2.*theta) + sigma_rt.*cos(2.*theta);




%% Numerical solution

n_np = length(node_coor_1(:,1));
n_en = 4;
n_el = length(IEN_1(:,1));
hh = 2/sqrt(length(IEN_1(:,1))/2);

E = 1e9; % Young's Modulus
pr = 0.3; % Poison's ratio
lambda = pr*E/(1+pr)/(1-2*pr);
mu = E/2/(1+pr);
% D = (E/(1-pr^2)).*[1,pr,0;pr,1,0;0,0,(1-pr)/2];
D = [lambda + 2*mu,lambda,0;lambda,lambda + 2*mu,0;0,0,mu];



% Gauss quadrature
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% Boundary conditions 
g = 0;
f_11 = sigma_11;
f_22 = sigma_22;
f_12 = sigma_12;

% ID array 
ID = zeros(n_np,2); % 2 dof
counter = 0;
for i = 1:n_np
    if i>67 && i<83
        counter = counter + 1;
        ID(i,1) = counter;
        counter = counter + 1;
        ID(i,2) = counter;
    else 
        ID(i,1) = 0;
        ID(i,2) = 0;
    end
    if i>97
        counter = counter + 1;
        ID(i,1) = counter;
        counter = counter + 1;
        ID(i,2) = counter;
    end
end
n_eq = counter;

% LM array
ID1 = ID(:,1);
ID2 = ID(:,2);
LM1 = ID1(IEN_1);
LM2 = ID2(IEN_1);
LM = zeros(n_el,n_en*2);
for aa = 1 : n_en
    LM(:,2*aa-1) = LM1(:,aa);
    LM(:,2*aa) = LM2(:,aa);
end


% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el
    x_ele = node_coor_1(IEN_1(ee, 1:n_en),1);
    y_ele = node_coor_1(IEN_1(ee, 1:n_en),2);

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

            f_ele_11(aa) = f_ele_11(aa) + weight(ll) * detJ * f_11(IEN_1(ee,aa)) * Na;
            f_ele_22(aa) = f_ele_22(aa) + weight(ll) * detJ * f_22(IEN_1(ee,aa)) * Na;
            f_ele(2*aa-1) = f_ele_11(aa);
            f_ele(2*aa) = f_ele_22(aa);

            B = zeros(3,2);
            for bb = 1:n_en
                Na = Quad(aa, xi(ll), eta(ll));
                [Na_xi, Na_eta] = Quad_grad(aa, xi(ll), eta(ll));
                Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                Na_y = (-Na_xi * dx_deta + Na_eta * dx_dxi) / detJ;
                B(1,2) = Na_x;
                B(2,2) = Na_y;
                B(3,1) = Na_y;
                B(3,2) = Na_x;
                k_ij = B'*D*B;
                k_ele(2*aa-1:2*aa,2*bb-1:2*bb) = k_ij;
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




%%

plot(node_coor_1(:,1),node_coor_1(:,2),'o')
hold on
plot(node_coor_1(:,1)+disp(:,1),node_coor_1(:,2)+disp(:,2),'.r')
 
% scatter(node_coor_1(:,1),node_coor_1(:,2),20,disp(:,2),'filled');
% colorbar;
% xlabel("x");
% ylabel("y");
% title("disp");
% colormap jet;
% caxis([0,1e-8]);
% axis equal;
% grid on;

    
 
