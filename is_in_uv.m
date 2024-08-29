function in_uv = is_in_uv(uv,img)
    %此次默认uv结果原点为左下角
    j = round(uv(1)*2047) + 1;
    i = round((1-uv(2))*2047) + 1;
    if(all(squeeze(img(i,j,:))==[50 50 50]'))
        in_uv = 0;
    else
        in_uv = squeeze(img(i,j,:));
    end
end