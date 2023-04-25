function dat = f_Clean_Mathots_way_v2(data, rej_criteria)
    
    data_dnoise = movmean(data,11,2); %smoothing
    dy_dt = zeros(size(data_dnoise));
    data_all = zeros(size(data_dnoise));
    datai = [];  %zeros(size(data_dnoise));
    recon_first = zeros(size(data_dnoise));
    removed = cell(1,size(data_dnoise,1));
    percent_rem = zeros(1,size(data_dnoise,1));
    selected_trial = [];

    %bounds
    ub = 3;
    lb = -3;

    %cutoff
    continuity = 20;
        
    %main manipulation
    for windw = 1:size(data_dnoise,1)

        %slopes
        dy_dt(windw,:) = gradient(data_dnoise(windw,:))...
                        ./gradient(1:size(data_dnoise,2));

        temp = data_dnoise(windw,:);
        temp(find(dy_dt(windw,:) < lb | dy_dt(windw,:) > ub)) = nan;
        temp(find(data_dnoise(windw,:) == 0)) = nan;

        zs1 = find(isnan(temp));

        if isempty(zs1) == 1
            temp2 = temp;
            data_all(windw,:) = temp2;
            datai(end+1,:) = temp2;
            selected_trial(end+1,:) = windw;
            continue;
        else
            %Ckeck gaps
            time_rng = f_discontinuties(zs1,continuity);
        end 
        
        %time_rng_new = zeros(size(time_rng));
        
        % +- 20 curoff
        for j= 1:size(time_rng,1)
            r1 = time_rng(j,1) - continuity;
            r2 = time_rng(j,2) + continuity;

            if r1 < 1
               r1 = 1;
            end

            if r2 > length(temp)
               r2 = length(temp);
            end
            %time_rng_new(j,1:2) = [r1,r2];

            temp(r1:r2) = nan;
        end


         if isnan(temp(1))
           %temp(1) = nanmean(temp);
           var = temp(~isnan(temp));

           if isempty(var) == 1 %the signal is a strightline (maybe a faulty signal)
              data_all(windw,:) = temp;
              continue;
           end

           temp(1) = var(1); %replace NaN with first legitimate val
         end


        if isnan(temp(end))
           %temp(end) = nanmean(temp);
           var = temp(~isnan(temp));
           temp(end) = var(end); %replace NaN with last legitimate val
        end


        xxi = 1:size(data,2);
        xx = 1:size(data,2);
        zs = find(isnan(temp));
        temp1 = temp;
        temp1(zs)=[];
        xx(zs)=[];
        temp2 = interp1(xx, temp1, xxi); %recon


        data_all(windw,:) = temp2;
        
        removed{windw} = find(data_dnoise(windw,:) ~= temp2);
        percent_rem(windw) = (size(removed{windw},2)/size(data_dnoise(windw,:),2))*100;

        if percent_rem(windw) < rej_criteria
           datai(end+1,:) = temp2;
           selected_trial(end+1,:) = windw;
        end

    end

trl_raus = ones(size(data_dnoise,1),1);
trl_raus(selected_trial) = 0; % 1=rejected, 0=selected

dat.data_all = data_all;
dat.datai = datai;
dat.trl_raus = trl_raus;
dat.auto_select = selected_trial;
dat.removed = removed;
dat.percent_rem = percent_rem;

end %function