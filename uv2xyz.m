function v3d = uv2xyz(uv)
    load('uv2xyz_info_li.mat','vertices','uvs','dis_uv','vsize','facets','img','max_fl');

    v3d = [nan,nan,nan];
    v_rgb = is_in_uv(uv,img);
    if (~v_rgb == 0)
        for i = 1:1:vsize
            dis_uv(i) = norm(uvs(i,:)-uv);
        end
        dis_uv = [dis_uv(dis_uv<max_fl),find(dis_uv<max_fl)];
        dis_uv = sortrows(dis_uv,1);
        
        indd = 1;
        facets_ind =zeros(size(facets,1),1);
        for i = 1:1:size(dis_uv,1)
            v_uv_ind = dis_uv(i,2);
            facet_ind = find(facets(:,4)==v_uv_ind|facets(:,5)==v_uv_ind|facets(:,6)==v_uv_ind);
            for n = 1:1:size(facet_ind,1) 
                if(~ismember(facet_ind(n),facets_ind(1:indd)))
                    v_tri = uvs(facets(facet_ind(n),4:6),:);
                    in_tri = is_in_tri(uv,v_tri);
                    if (~in_tri == 0)
                        break;
                    end
                    facets_ind(indd) = facet_ind(n);
                    indd = indd + 1;
                end
            end
             if (~in_tri == 0)
                break;
             end
        end
        if (~in_tri == 0)
            x_tri = vertices(facets(facet_ind(n),1:3),:);
            A = x_tri(1,:);
            B = x_tri(2,:);
            C = x_tri(3,:);
            AC = C-A;
            AB = B-A;
            v3d = AC.*in_tri(1) + AB.*in_tri(2) + A;
        end
    end
end