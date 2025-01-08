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
%% Manufactured Solution to determine BC on the right surface, top surface fixed u = 0

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
r = sqrt((right_coor(:,1)-(-1)).^2 + (right_coor(:,2)-(-1)).^2);
theta = atan((right_coor(:,2)-(-1))./((right_coor(:,1)-(-1))));
sigma_rr_bc = sigma_rr(Tx,R,r,theta);
sigma_tt_bc = sigma_tt(Tx,R,r,theta);
sigma_rt_bc = sigma_rt(Tx,R,r,theta);
% Boundary condition on the right surface
sigma_11_bc = (sigma_rr_bc + sigma_tt_bc)./2 + (sigma_rr_bc - sigma_tt_bc)./2.*cos(2.*theta) - sigma_rt_bc.*sin(2.*theta);
sigma_22_bc = (sigma_rr_bc + sigma_tt_bc)./2 - (sigma_rr_bc - sigma_tt_bc)./2.*cos(2.*theta) + sigma_rt_bc.*sin(2.*theta);
sigma_12_bc = (sigma_rr_bc - sigma_tt_bc)./2.*sin(2.*theta) + sigma_rt_bc.*cos(2.*theta);




%% Numerical solution

n_np = length(node_coor_1(:,1));
n_en = 4;
n_el = length(IEN_1(:,1));
hh = 2/sqrt(length(IEN_1(:,1))/2);

E = 1e9; % Young's Modulus
mu = 0.3; % Poison's ratio

% Gauss quadrature
n_int_xi  = 3;
n_int_eta = 3;
n_int     = n_int_xi * n_int_eta;
[xi, eta, weight] = Gauss2D(n_int_xi, n_int_eta);

% Boundary conditions 
g = 0;
h_11 = sigma_11_bc;
h_22 = sigma_22_bc;
h_12 = sigma_12_bc;

% ID array 
ID = zeros(n_np,1);
counter = 0;
for i = 1:n_np
    if node_coor_1(i,2) ~= 1
        counter = counter + 1;
        ID(i) = counter;
    else
        ID(i) = 0;
    end
end
n_eq = counter;

% LM array
LM = ID(IEN_1);

% allocate the stiffness matrix and load vector
K = spalloc(n_eq, n_eq, 9 * n_eq);
F = zeros(n_eq, 1);

for ee = 1 : n_el

















end




    
 
