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
    IEN_1 = ien(8:519,2:5);
    % 关闭文件
    fclose(fid);

%% Read nodes and IEN for case 2
    % 打开文件
    fid = fopen('quarter-plate-with-hole-quad-Refine4.msh', 'r');

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
    node_coor_2 = coor(2:2147,:);
    IEN_2 = ien(8:2055,2:5);
    % 关闭文件
    fclose(fid);
%% Read nodes and IEN for case 3
    % 打开文件
    fid = fopen('quarter-plate-with-hole-quad-Refine5.msh', 'r');

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
    node_coor_3 = coor(2:8387,:);
    IEN_3 = ien(8:8199,2:5);
    % 关闭文件
    fclose(fid);





