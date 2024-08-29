% 读取BMP图像  
img = imread('atlas-horizons.bmp'); 

% 读取OBJ文件  
fid= fopen('SAM2.obj','r');  
uvs = [];
vertices = [];  
facets = [];

while ~feof(fid)
    str = fgetl(fid);
    
    s = regexp(str,' ','split');
    if s{1} == 'vt'    
        uvss1 = str2double(s{2});
        uvss2 = str2double(s{3});
        uvs = [uvs;uvss1 uvss2];
    elseif s{1} == 'v'
        verticess1 = str2double(s{2});
        verticess2 = str2double(s{3});
        verticess3 = str2double(s{4});
        vertices = [vertices;verticess1 verticess2 verticess3];
    elseif s{1} == 'f'
        facet = zeros(1,6);
        for i=2:1:4
            facets1=regexp(s{i},'\/','split');
            facet(i-1) = str2double(facets1{1});
            facet(i+2) = str2double(facets1{2});
        end
        facets =[facets;facet];
    end
end
% 提取顶点坐标和纹理坐�?  

fclose(fid);

facets_size = size(facets,1);
max_fl = 0;
for i= 1:1:facets_size
    max_fl_new = max([norm(uvs(facets(i,4),:)-uvs(facets(i,5),:)) norm(uvs(facets(i,4),:)-uvs(facets(i,6),:)) norm(uvs(facets(i,5),:)-uvs(facets(i,6),:))]);
    if (max_fl_new>max_fl)
        max_fl = max_fl_new;
    end
end
  
vsize = size(uvs,1);
dis_uv = zeros(vsize,1);

save('uv2xyz_info.mat','vertices','uvs','dis_uv','vsize','facets','img','max_fl');

% % 随机生成�?�?0-1之间的二维UV坐标  
% random_uv = rand(1, 2);  
% 
% v3d = uv2xyz(random_uv,vertices,uvs,dis_uv,vsize,facets,img,max_fl);
% 
% %遍历uv图生成带颜色的点�?
% bmp_size1 = size(img,1);
% bmp_size2 = size(img,2);
% xyz_rgb = [];
% tag = 0;
% pro = 0;
% t1 = tic();
% for i = 1:2:bmp_size1
%     for j = 1:2:bmp_size2
%         v3d = uv2xyz([(j-1)/(bmp_size2-1),1-(i-1)/(bmp_size1-1)],vertices,uvs,dis_uv,vsize,facets,img,max_fl);
%         if iscell(v3d)
%             xyz_rgb = [xyz_rgb;v3d{1} cast(v3d{2},'double')./255];
%         end
%     end
%     tag = tag + 1;
%     if tag >= ceil(bmp_size1/200)
%         pro = pro + 1;
%         disp(num2str(pro));
%         tag = 0;
% %         break;
%     end
% end
% t2 = toc(t1);
% save('SAM2_xyz_rgb_3.mat','xyz_rgb');
