M = 4; %天线数量
y_depth = 0.5e-3; %7倍趋肤深度
xigema_28GHz_SAM_liquid = 27.5704;%电导率
A_area = input('积分面积(mm2)=');%积分面积400mm2
dx = 0.5e-3;
dy = 0.5e-3;
dz = 0.5e-3;
% distance = 5e-3;
lx = 50e-3;
ly = 50e-3;
repeats = 1;
global PD_infos;
PD_infos = [];


load('E1_mask_resample_05.mat');
x_num = size(Axis0,2)-1;
y_num = size(Axis1,2)-1;
z_num = size(Axis2,2)-1;

x_min = 1;
x_max = size(Axis0,2)-1;
y_min = 1;
y_max = size(Axis1,2)-1;
z_min = 1;
z_max = size(Axis2,2)-1;

x_axis =zeros(x_num,1);y_axis = zeros(y_num,1);z_axis = zeros(z_num,1);
for n=1:1:x_num
    x_axis(n) = (Axis0(n)+Axis0(n+1))*0.5;
end
for n=1:1:y_num
    y_axis(n) = (Axis1(n)+Axis1(n+1))*0.5;
end
for n=1:1:z_num
    z_axis(n) = (Axis2(n)+Axis2(n+1))*0.5;
end

E1 = Snapshot0;

load('E2_mask_resample_05.mat');
E2 = Snapshot0;

load('E3_mask_resample_05.mat');
E3 = Snapshot0;

load('E4_mask_resample_05.mat');
E4 = Snapshot0;

e = [E1 E2 E3 E4];
e = reshape(e,[],3,4);
e_3d = reshape(e,x_num,y_num,z_num,3,4);

load('H1_mask_resample_05.mat');
H1 = Snapshot0;

load('H2_mask_resample_05.mat');
H2 = Snapshot0;

load('H3_mask_resample_05.mat');
H3 = Snapshot0;

load('H4_mask_resample_05.mat');
H4 = Snapshot0;

h = [H1 H2 H3 H4];
h = reshape(h,[],3,4);
h_3d = reshape(h,x_num,y_num,z_num,3,4);

t1 = tic;
% surf_xyz = [];
surf_xyz = zeros(x_num*y_num*z_num,3);
surf_ijk = zeros(x_num*y_num*z_num,3);
ijk_surfn = zeros(x_num,y_num,z_num);
n = 0;
for i = 2:1:x_num-1
    for j= 1:1:y_num-1
        for k = 1:1:z_num
            if ~isnan(e_3d(i,j,k,1,1))
               if isnan(e_3d(i-1,j,k,1,1)) || isnan(e_3d(i+1,j,k,1,1)) || isnan(e_3d(i,j-1,k,1,1)) || isnan(e_3d(i,j+1,k,1,1)) || isnan(e_3d(i,j,k-1,1,1)) || isnan(e_3d(i,j,k+1,1,1))
%                     surf_xyz = [surf_xyz; i j k];
                    n = n + 1;
                    surf_xyz(n,:) = [x_axis(i) y_axis(j) z_axis(k)];
                    surf_ijk(n,:) = [i j k];
                    ijk_surfn(i,j,k) = n;
                    continue;
                end
            end
        end
    end
end
surf_xyz = surf_xyz(1:n,:);
surf_ijk = surf_ijk(1:n,:);
t2 = toc(t1);
% load('surf_shell.mat');
surf_num = size(surf_xyz,1);

load('surf_fai.mat');
%计算A_im
t1 = tic;
y_depth_num = ceil(y_depth/dy); %先考虑沿y轴向内积分
% A_const = 0.5*xigema_28GHz_SAM_liquid*dy; %电导率数值待确认
A_const = 0.5; 
A_im = zeros(surf_num,M,M);
for n = 1:1:surf_num
    i = surf_ijk(n,1);
    j = surf_ijk(n,2);
    k = surf_ijk(n,3);
%     if n==75393
%         disp('81960');
%     end
    while isnan(e_3d(i,j,k,1,1))
        j=j+1;
        if j > y_max
            break;
        end
    end
    
    for j_plus = 1:y_depth_num
        if j+y_depth_num-1 > y_max
            break;
        end
        e_peaks = squeeze(e_3d(i,j+j_plus-1,k,:,:));
        h_peaks = squeeze(h_3d(i,j+j_plus-1,k,:,:));
        h_peaks = h_peaks';
        for i = 1:M
            for j = 1:M
                poin_peaks(i,j,:) = cross(e_peaks(:,i),h_peaks(j,:));
            end
        end
%         if any(any(any(isnan(poin_peaks))))
%             disp('nan');
%         end
%         poin_peaks = real(poin_peaks);
        A_im_zslice = dot(poin_peaks,repmat(reshape(surf_fai(n,:)',1,1,3),M),3);
%         if any(any(isnan(A_im_zslice)))
%             disp('nan');
%         end
        A_im_zslice = (A_im_zslice + A_im_zslice') * 0.5;
%         if j_plus==1 && A_im(i,k,1,1)~=0
%             disp(num2str([i j k]));
%         end

        A_im(n,:,:) = A_im(n,:,:) + reshape(A_im_zslice,[1,M,M]);
    end
    A_im(n,:,:) = A_im(n,:,:) .* A_const;
end
t3 = toc(t1);
save('PD_info.mat','A_im','M');

L_k = sqrt(A_area)*1e-3;
R_k = sqrt(A_area/pi)*1e-3;
M_k = fix(L_k / dx);
N_k = fix(L_k / dy);

save('infos','surf_xyz','surf_ijk','ijk_surfn','x_axis','y_axis','z_axis','L_k','R_k','M_k','N_k');

% load('PD_info.mat')
num_x = size(A_im,1);
num_z = size(A_im,2);

%PSO 优化
func = @get_APD;
nvars = 2;
lb = [0 0];
ub = [1 1];

options = optimoptions('particleswarm','DisplayInterval',1,'Display','iter','MaxStallIterations',10);

times=zeros(repeats,1);
APDs=zeros(repeats,1);
APD_infos=cell(repeats,1);
for i=1:1:repeats
    disp(['repeat:',num2str(i)]);
    t1 = tic;
    [APD_loc,APD]= particleswarm(func,nvars,lb,ub,options);
    t2 = toc(t1);
    times(i)=t2;
    APD_infos{i}=SDR_w(APD_loc(:));
    APDs(i) = 1 / APD;
end
save(['PSO_SDR_22_patch_',num2str(A_area*1e-2),'cm2_ra_EH.mat'],'APDs','APD_infos','times');
save('PD_infos_EH.mat','PD_infos');
% APD = get_APD(64,100,0);

% disp(APD);


function APD = get_APD(p_Positon) 
    PD_info =  SDR_w(p_Positon);
    APD = 1 / PD_info{end};
end

function PD_info = SDR_w(p_Positon)
    load('PD_info.mat','A_im','M');
    load('infos','surf_xyz','ijk_surfn','x_axis','y_axis','z_axis','R_k','M_k','N_k');
    global PD_infos;
    
%     surf_num = size(surf_xyz,1);
    x_min = min(surf_xyz(:,1));
    x_max = max(surf_xyz(:,1));
    y_min = min(surf_xyz(:,2));
    y_max = max(surf_xyz(:,2));
    z_min = min(surf_xyz(:,3));
    z_max = max(surf_xyz(:,3));

    u = p_Positon(1);v = p_Positon(2);
    ind = 0;
    A_im_n = zeros(M_k*N_k,M,M);
    
    pre_v = uv2xyz([u,v]).*1e-3;
    if(~any(pre_v)||pre_v(1)<=x_min||pre_v(1)>=x_max||pre_v(2)<=y_min||pre_v(2)>=y_max||pre_v(3)<=z_min||pre_v(3)>=z_max)
        PD_info = {0,pre_v,[u,v],0};
%         PD_infos = [PD_infos;PD_info];
    else
%         for nn = 1:1:surf_num
%             ab_v = surf_xyz(nn,:);
%             if norm(ab_v - pre_v) <= R_k
%                 ind = ind + 1;
%                 A_im_n(ind,:,:) = A_im(nn,:,:);
%             end
%         end
        st_x = find(x_axis>=(pre_v(1)-R_k),1);
        en_x = find(x_axis>=(pre_v(1)+R_k),1)-1;
        st_y = find(y_axis>=(pre_v(2)-R_k),1);
        en_y = find(y_axis>=(pre_v(2)+R_k),1)-1;
        st_z = find(z_axis>=(pre_v(3)-R_k),1);
        en_z = find(z_axis>=(pre_v(3)+R_k),1)-1;
        
        if isempty(en_x)
            en_x = size(x_axis,1);
        end
        if isempty(en_y)
            en_y = size(y_axis,1);
        end
        if isempty(en_z)
            en_z = size(z_axis,1);
        end

        

        for ii = st_x:en_x
            for jj = st_y:en_y
                for kk = st_z:en_z
                    if ijk_surfn(ii,jj,kk)==0
                        continue;
                    else
                        nn = ijk_surfn(ii,jj,kk);
                        ab_v = surf_xyz(nn,:);
                        if norm(ab_v - pre_v) <= R_k
                            ind = ind + 1;
                            A_im_n(ind,:,:) = A_im(nn,:,:);
                        end
                    end
                end
            end
        end
        A_im_n = A_im_n(1:ind,:,:);
        A_im_sub_sum = squeeze(sum(A_im_n,1));
        A_im_sub_sum = A_im_sub_sum./size(A_im_n,1);
%         A_im_sub_sum = complex(A_im_sub_sum);


        %半定规划SDR
        cvx_begin quiet
%         cvx_solver sdpt3
            variable W(M,M) hermitian semidefinite
            maximize (trace(A_im_sub_sum*W))
            subject to
                trace(W)==1;
                W(1,1) <= 0.5;
                W(2,2) <= 0.5;
                W(3,3) <= 0.5;
                W(4,4) <= 0.5;
        cvx_end

        W_best = W;
        [V,D]=eigs(W_best);
        max_la_num = find(max(D)==max(max(D)),1);
        max_la = D(max_la_num,max_la_num);
        q = V(:,max_la_num);
        w = sqrt(max_la)*q;%天线复权重结果

        %求mpsPD结果
        APD = w'*A_im_sub_sum*w;

        PD_info = {w,pre_v,[u,v],APD};
        PD_infos = [PD_infos;PD_info];
    end
end