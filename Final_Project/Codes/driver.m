clear;clc; close all;


    % 打开文件
    fid = fopen('quarter-plate-with-hole-quad-Refine3.msh', 'r');

    if fid == -1
        error('无法打开文件');
    end
    
    % 初始化
    nodes = [];
    elements = [];
    level = 0;
    i = 1;
    coor = zeros(562,3);
    % 读取文件内容
   while ~feof(fid)
        line = fgetl(fid);
        num  = strsplit(line);
        num = str2double(num(1,:));
        isdata = all(isnan(num));
        if isdata == 0
            if length(num) == 3
                coor(i,:) = num;
                i = i+1;
            end
        end
    end
    node_coor = coor(2:length(coor)-7,:);
    % 关闭文件
    fclose(fid);





