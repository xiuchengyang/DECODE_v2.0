function [rec_cg,clrx,clry] = DECODEPixel(sdate, line_t, qa,pos, T_cg,Tmax_cg,conse,num_c,varargin)
%
% Inputs:
% stk_n='stack'; stack image name
% ncols = 8021; % number of pixels processed per line
% nrows=1; % the nrowsth lines
% for example    1 2 3 4 5
%                6 7 8 9 10
%
% Outputs:
%
% rec_cg RECord information about all curves between ChanGes
% rec_cg(i).t_start record the start of the ith curve fitting (julian_date)
% rec_cg(i).t_end record the end of the ith curve fitting (julian_date)
% rec_cg(i).t_break record the first observed break time (julian_date)
% rec_cg(i).coefs record the coefficients of the ith curve
% rec_cg(i).pos record the position of the ith pixel (pixel id)
% rec_cg(i).magnitude record the change vector of all spectral bands
% rec_cg(i).category record what fitting procudure and model is used
% cateogry category 5x: persistent snow    4x: Fmask fails
% cateogry category 3x: modified fit       2x: end fit
% category category 1x: start fit           x: normal procedure
% cateogry category x1: mean value         x4: simple model
% category category x6: advanced model     x8: full model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  defining variables
    warning('off', 'all'); % no warnings for linear fit
    %% Constants
    p = inputParser;
    addParameter(p,'recordSR', 1); % number of consec
    addParameter(p,'tide',[]); % Tide influence or not;
    addParameter(p,'numBandGreen',2);
    addParameter(p,'numBandSWIR1',9);
    addParameter(p,'bandNames',[]);
    addParameter(p,'B_detect',1:size(line_t,2)); % Use all bands
    parse(p,varargin{:});
    recordSR = p.Results.recordSR; % Whether to record the observation values
    tide = p.Results.tide;
    numBandGreen = p.Results.numBandGreen;
    numBandSWIR1 = p.Results.numBandSWIR1;
    bandNames = p.Results.bandNames;
    B_detect = p.Results.B_detect;
    nbands = size(line_t,2);

    if isempty(tide)
        isTide = false;  % If empty, then dont consider the tide fluctuation
    else
        isTide = true;
    end

    % maximum number of coefficient required
    % 2 for tri-modal; 2 for bi-modal; 2 for seasonality; 2 for linear; 1 for
    % water level
    if isTide
        min_num_c = 5;
        mid_num_c = 7;
        max_num_c = 9;
    else
        min_num_c = 4;
        mid_num_c = 6;
        max_num_c = 8;
    end
    % number of clear observation / number of coefficients
    n_times = 3;
    % update frequency
    prct_update_model = 0.03; % 3% of the previous model fit
    % initialize NUM of Functional Curves
    num_fc = 0;
    % number of days per year
    num_yrs = 365.25;
    % Threshold for cloud, shadow, and snow detection.
    T_const = 4.42;
    % minimum year for model intialization
    mini_yrs = 1;
    % no change detection for permanent snow pixels
    t_sn = 0.75;
    % threshold (degree) of mean included angle
    nsign = 45;
    
    % Tmasking of noise
    Tmax_cg = chi2inv(Tmax_cg,length(B_detect));
    % adjust threshold based on chi-squared distribution;
    T_cg = chi2inv(T_cg,length(B_detect));
    
    % initialize the struct data of RECording of ChanGe (rec_cg)
    rec_cg = struct('t_start',[],'t_end',[],'t_break',[],'coefs',[],'rmse',[],...
        'pos',[],'change_prob',[],'num_obs',[],'category',[],'magnitude',[],'durchange',[],'data',[]);

    
    % mask data
    idrange = qa<254; % Exclude the saturated pixels
    
    % # of clear observatons
    idclr = qa<=1; % Clear land and water


    % snow pixels
    idsn = qa == 3;

    % percent of snow observations
    sn_pct = sum(idsn)/(sum(idclr)+sum(idsn)+eps);
    
    % not enough clear observations for change detection
    if sum(idclr) < n_times*max_num_c
        % permanent snow pixels
        if sn_pct > t_sn
            % snow observations are "good" now
            idgood = idsn|idclr;
            % number of snow pixel within range
            n_sn = sum(idgood);
            
            if n_sn < n_times*min_num_c % not enough snow pixels
                % Xs & Ys for computation
                clrx = sdate(idgood);
                % bands 1-5,7,6
                clry = line_t(idgood,1:end);
                clry = double(clry);
                
                return
            else
                % Xs & Ys for computation
                clrx = sdate(idgood);
                % bands 1-5,7,6
                clry = line_t(idgood,1:end);
                clry = double(clry);
                            
                % find repeated ids
                [clrx,uniq_id,~] = unique(clrx);
                % mean of repeated values
                tmp_y = zeros(length(clrx),nbands);
                
                % get the mean values
                for i = 1:nbands
                    tmp_y(:,i) = clry(uniq_id,i);
                end
                clry = tmp_y;
                
                if isTide
                    % tide info 
                    tide = double(tide(idgood));
                    tmp_wl = tide(uniq_id);
                    tide = tmp_wl;
                end

                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands);
                % rmse for each band
                rmse = zeros(nbands,1);
                % snow qa = 50
                qa = 50;
                
                for i_B=1:nbands
                    i_span = length(clrx);
                    if i_span < min_num_c*n_times % fill value for frequently saturated snow pixels
                        fit_cft(1,i_B) = 10000; % give constant value
                    else % fit for enough snow pixels
                        if isTide
                            [fit_cft(:,i_B),rmse(i_B)]=autoTSFitWL(clrx,clry(:,i_B),min_num_c,tide);
                        else
                            [fit_cft(:,i_B),rmse(i_B)]=autoTSFit(clrx,clry(:,i_B),min_num_c);
                        end
                    end
                end
                
                % updating information at each iteration
                % record time of curve start
                rec_cg(num_fc).t_start=clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0; % no break at the moment
                % record postion of the pixel
                rec_cg(num_fc).pos = pos;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = n_sn;
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands);
                % record durchange
                rec_cg(num_fc).durchange = zeros(2,nbands);

                if recordSR
                    if isTide
                        data = [clrx,clry,tide];
                        if ~isempty(bandNames)
                            rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                        else
                            rec_cg(num_fc).data = data;
                        end
                    else
                        data = [clrx,clry];
                        if ~isempty(bandNames)
                            rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                        else
                            rec_cg(num_fc).data = data;
                        end
                    end
                end
            end
        else % no change detection for clear observations
            % within physical range pixels
            idgood =  idrange;
            
            % Xs & Ys for computation
            clrx = sdate(idgood);
            % bands 1-5,7,6
            clry = line_t(idgood,:);
            clry = double(clry);
             
            % find repeated ids
            [clrx,uniq_id,~] = unique(clrx);
            % mean of repeated values
            tmp_y = zeros(length(clrx),nbands);
            
            % get the mean values
            for i = 1:nbands
                tmp_y(:,i) = clry(uniq_id,i);
            end
            clry = tmp_y;

            if isTide
                tide = double(tide(idgood));
                tmp_wl = tide(uniq_id);
                tide = tmp_wl;
            end

            idclr = clry(:,numBandGreen) < median(clry(:,numBandGreen)) + 400;
            n_clr = sum(idclr);
            
            if n_clr < n_times*min_num_c % not enough clear pixels
                return
            else
                % Xs & Ys for computation
                clrx = clrx(idclr);
                clry = clry(idclr,:);

                if isTide
                    tide = tide(idclr);
                end
                % the first observation for TSFit
                i_start = 1;
                % identified and move on for the next curve
                num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
                
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands);
                % rmse for each band
                rmse = zeros(nbands,1);
                % Fmask fail qa = 40
                qa = 40;
                
                for i_B = 1:nbands
                    % fit basic model for all within range snow pixels
                    if isTide
                        [fit_cft(:,i_B),rmse(i_B)] = autoTSFitWL(clrx,clry(:,i_B),min_num_c,tide);
                    else
                        [fit_cft(:,i_B),rmse(i_B)] = autoTSFit(clrx,clry(:,i_B),min_num_c);
                    end
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end = clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = pos;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx);
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands);
                % record during change
                rec_cg(num_fc).durchange = zeros(2,nbands);

                if recordSR
                    if isTide
                        data = [clrx,clry,tide];
                        if ~isempty(bandNames)
                            rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                        else
                            rec_cg(num_fc).data = data;
                        end
                    else
                        data = [clrx,clry];
                        if ~isempty(bandNames)
                            rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                        else
                            rec_cg(num_fc).data = data;
                        end
                    end
                end
            end
        end
    else % normal CCDC procedure
        % clear and within physical range pixels
        idgood = idclr & idrange;

        % Xs & Ys for computation
        clrx = sdate(idgood);
        % bands 1-5,7,6
        clry = line_t(idgood,1:end);
        clry = double(clry);  
        
        % find repeated ids
        [clrx,uniq_id,~] = unique(clrx);
        
        % continue if not enough clear pixels
        if length(clrx) < n_times*min_num_c+conse
            return
        end
        
        % mean of repeated values
        tmp_y = zeros(length(clrx),nbands);
        % tmp_wl = zeros(length(clrx),1);
        
        % get the mean values
        for i = 1:nbands
            tmp_y(:,i) = clry(uniq_id,i);
        end
        clry = tmp_y;
        clear tmp_y;
        if isTide
            tide = double(tide(idgood));
            tmp_wl = tide(uniq_id);
            tide = tmp_wl;
            clear tmp_wl;
        end
        % caculate median variogram
        var_clry = clry(2:end,:)-clry(1:end-1,:);
        adj_rmse = median(abs(var_clry),1);
        
        % start with the miminum requirement of clear obs
        i = n_times*min_num_c;
        
        % initializing variables
        % the first observation for TSFit
        i_start = 1;
        % record the start of the model initialization (0=>initial;1=>done)
        BL_train = 0;
        % identified and move on for the next curve
        num_fc = num_fc + 1; % NUM of Fitted Curves (num_fc)
        % record the num_fc at the beginning of each pixel
        rec_fc = num_fc;
        % initialize i_dense
        i_dense = 1;
        
        if recordSR
            % Initilization of for the clear observation recording
            i_start_segment = 1;
        end

        % while loop - process till the last clear observation - conse
        while i<= length(clrx)-conse
            % span of "i"
            i_span = i-i_start+1;
            % span of time (num of years)
            time_span = (clrx(i)-clrx(i_start))/num_yrs;
            % max time difference
            time_dif = max(clrx(i_start+1:i) - clrx(i_start:i-1));
            
            % basic requrirements: 1) enough observations; 2) enough time
            if i_span >= n_times*min_num_c && time_span >= mini_yrs
                % initializing model
                if BL_train == 0
                    % check max time difference
                    if time_dif > num_yrs
                        i = i+1;
                        i_start = i_start+1;
                        % i that is dense enough
                        i_dense = i_start;
                        continue
                    end
                    % Tmask: noise removal (good => 0 & noise => 1)
                    blIDs = autoTmask(clrx(i_start:i+conse),clry(i_start:i+conse,[numBandGreen,numBandSWIR1]),...
                        (clrx(i+conse)-clrx(i_start))/num_yrs,adj_rmse(numBandGreen),adj_rmse(numBandSWIR1),T_const);
                    
                    % IDs to be removed
                    IDs = i_start:i+conse;
                    rmIDs = IDs(blIDs(1:end-conse) == 1);
                    
                    % update i_span after noise removal
                    i_span = sum(~blIDs(1:end-conse));
                    
                    % check if there is enough observation
                    if i_span < n_times*min_num_c
                        % move forward to the i+1th clear observation
                        i = i+1;
                        % not enough clear observations
                        continue;
                        % check if there is enough time
                    else
                        % copy x & y
                        cpx = clrx;
                        cpy = clry;
                        
                        % remove noise pixels between i_start & i
                        cpx(rmIDs) = [];
                        cpy(rmIDs,:) = [];
                        if isTide
                            cpwl = tide;
                            cpwl(rmIDs) = [];
                        end
                        % record i before noise removal
                        % This is very important as if model is not initialized
                        % the multitemporal masking shall be done again instead
                        % of removing outliers in every masking
                        i_rec = i;
                        
                        % update i afer noise removal (i_start stays the same)
                        i = i_start+i_span-1;
                        % update span of time (num of years)
                        time_span = (cpx(i)-cpx(i_start))/num_yrs;
                        
                        % check if there is enough time
                        if time_span < mini_yrs
                            % keep the original i
                            i = i_rec;
                            % move forward to the i+1th clear observation
                            i = i+1;
                            % not enough time
                            continue;
                            % Step 2: model fitting
                        else
                            % remove noise
                            clrx = cpx;
                            clry = cpy;
                            if isTide
                                tide = cpwl;
                            end
                            % Step 2: model fitting
                            % initialize model testing variables
                            % defining computed variables
                            fit_cft = zeros(max_num_c,nbands);
                            % rmse for each band
                            rmse = zeros(nbands,1);
                            % value of differnce
                            v_dif = zeros(nbands,1);
                            % record the diference in all bands
                            rec_v_dif = zeros(i-i_start+1,nbands);
                            
                            % update number of coefficients
                            update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                            
                            for i_B = 1:nbands
                                % initial model fit
                                if isTide
                                    [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                        autoTSFitWL(clrx(i_start:i),clry(i_start:i,i_B),update_num_c,tide(i_start:i));
                                else
                                    [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                        autoTSFit(clrx(i_start:i),clry(i_start:i,i_B),update_num_c);
                               
                                end
                            end
                            
                            
                            % normalized to z-score
                            for i_B = B_detect
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                
                                % compare the first clear obs
                                v_start = rec_v_dif(1,i_B)/mini_rmse;
                                % compare the last clear observation
                                v_end = rec_v_dif(end,i_B)/mini_rmse;
                                % anormalized slope values
                                v_slope = fit_cft(2,i_B)*(clrx(i)-clrx(i_start))/mini_rmse;
                                
                                % differece in model intialization
                                v_dif(i_B) = abs(v_slope) + max(abs(v_start),abs(v_end));
                            end
                            v_dif = norm(v_dif(B_detect))^2;
                            
                            % find stable start for each curve
                            if v_dif > T_cg
                                % start from next clear obs
                                i_start = i_start + 1;
                                % move forward to the i+1th clear observation
                                i = i + 1;
                                % keep all data and move to the next obs
                                continue
                            else
                                % model ready!
                                BL_train = 1;
                                % count difference of i for each iteration
                                i_count = 0;
                                % label 0 for computing the update frequency of model fit
                                num_skip = 0;
                                
                                % find the previous break point
                                if num_fc == rec_fc
                                    % first curve
                                    i_break = 1;
                                else
                                    % after the first curve
                                    i_break = find(clrx >= rec_cg(num_fc-1).t_break);
                                    i_break = i_break(1);
                                end
                                
                                if i_start > i_break
                                    % model fit at the beginning of the time series
                                    for i_ini = i_start-1:-1:i_break
                                        if i_start - i_break < conse
                                            ini_conse = i_start - i_break;
                                        else
                                            ini_conse = conse;
                                        end
                                        % value of difference for conse obs
                                        v_dif = zeros(ini_conse,nbands);
                                        % record the magnitude of change
                                        v_dif_mag = v_dif;
                                        % record the date
                                        v_dif_clrx = zeros(ini_conse,1);
                                        % chagne vector magnitude
                                        vec_mag = zeros(ini_conse,1);
                                        
                                        for i_conse = 1:ini_conse
                                            for i_B = 1:nbands
                                                % absolute difference
                                                if isTide
                                                    v_dif_mag(i_conse,i_B) = clry(i_ini-i_conse+1,i_B)-autoTSPredWL(clrx(i_ini-i_conse+1),fit_cft(:,i_B),tide(i_ini-i_conse+1));
                                                else
                                                    v_dif_mag(i_conse,i_B) = clry(i_ini-i_conse+1,i_B)-autoTSPred(clrx(i_ini-i_conse+1),fit_cft(:,i_B));
                                            
                                                end
                                                % normalized to z-scores
                                                if sum(i_B == B_detect)
                                                    % minimum rmse
                                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                                    
                                                    % z-scores
                                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/mini_rmse;
                                                end
                                            end
                                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                                            v_dif_clrx(i_conse) = clrx(i_ini-i_conse+1);
                                        end
                                        
                                        % get the vec sign
                                        % vec_sign = abs(sum(sign(v_dif(:,B_detect)),1));
                                        max_angl = mean(angl(v_dif(:,B_detect)));
                                        
                                        % change detection
                                        if min(vec_mag) > T_cg && max_angl < nsign% change detected
                                            break
                                        elseif vec_mag(1) > Tmax_cg % false change
                                            % remove noise
                                            clrx(i_ini) = [];
                                            clry(i_ini,:) = [];
                                            if isTide
                                                tide(i_ini) = [];
                                            end
                                            i=i-1; % stay & check again after noise removal
                                        end
                                        
                                        % update new_i_start if i_ini is not a confirmed break
                                        i_start = i_ini;
                                    end
                                end
                                % only fit first curve if 1) have more than
                                % conse obs 2) previous obs is less than a year
                                if num_fc == rec_fc && i_start - i_dense >= conse
                                    % defining computed variables
                                    fit_cft = zeros(max_num_c,nbands);
                                    % rmse for each band
                                    rmse = zeros(nbands,1);
                                    % start fit qa = 10
                                    qa = 10;
                                    
                                    update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                                    
                                    for i_B=1:nbands
                                        if isTide
                                            [fit_cft(:,i_B),rmse(i_B)] = ...
                                                autoTSFitWL(clrx(i_dense:i_start-1),clry(i_dense:i_start-1,i_B),update_num_c,tide(i_dense:i_start-1));
                                        else
                                            [fit_cft(:,i_B),rmse(i_B)] = ...
                                                autoTSFit(clrx(i_dense:i_start-1),clry(i_dense:i_start-1,i_B),update_num_c);
                                       
                                        end
                                    end
                                    % record time of curve end
                                    rec_cg(num_fc).t_end = clrx(i_start-1);
                                    % record postion of the pixel
                                    rec_cg(num_fc).pos = pos;
                                    % record fitted coefficients
                                    rec_cg(num_fc).coefs = fit_cft;
                                    % record rmse of the pixel
                                    rec_cg(num_fc).rmse = rmse;
                                    % record break time
                                    rec_cg(num_fc).t_break = clrx(i_start);
                                    % record change probability
                                    rec_cg(num_fc).change_prob = 1;
                                    % record time of curve start
                                    rec_cg(num_fc).t_start = clrx(1);
                                    % record fit category
                                    rec_cg(num_fc).category = qa + min_num_c;
                                    % record number of observations
%                                     rec_cg(num_fc).num_obs = i_start - i_dense;
                                    rec_cg(num_fc).num_obs = i_start - 1;
                                    % record change magnitude
                                    rec_cg(num_fc).magnitude = - median(v_dif_mag,1);
                                    % record during-change
                                    rec_cg(num_fc).durchange = durchange(nbands, v_dif_clrx, - v_dif_mag);

                                    if recordSR
                                        if isTide
                                            data = [clrx,clry,tide];
                                            data = data(1:i_start-1,:);
                                            if ~isempty(bandNames)
                                                rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                                            else
                                                rec_cg(num_fc).data = data;
                                            end
                                        else
                                            data = [clrx,clry];
                                            data = data(1:i_start-1,:);
                                            if ~isempty(bandNames)
                                                rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                                            else
                                                rec_cg(num_fc).data = data;
                                            end
                                        end
                                    end
                                    % identified and move on for the next functional curve
                                    num_fc = num_fc + 1;


                                    if recordSR
                                        % Initilization of for the clear observation recording
                                        i_start_segment = i_start;
                                    end
                                end
                            end
                        end
                    end
                end
                
                % continuous monitoring started!!!
                if BL_train == 1
                    % all IDs
                    IDs = i_start:i;
                    % span of "i"
                    i_span = i-i_start+1;
                    % # of observations skipping to fit model
                    if num_skip == 0 % if is 0, not ready to skip for this observation, and thus to compute when we can start to skip
                        num_skip = fix(i_span*prct_update_model); % close to lower value
                    end
                    
                    % determine the time series model
                    update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                    
                    % initial model fit when there are not many obs
                    if  i_count == 0 || i_span <= max_num_c*n_times
                        % update i_count at each interation
                        i_count = clrx(i)-clrx(i_start);
                        
                        % record previous end i
                        pre_end = i;
                        % label as 0, indicating for being not ready to skip next observations
                        num_skip = 0; 
                        
                        % defining computed variables
                        fit_cft = zeros(max_num_c,nbands);
                        % rmse for each band
                        rmse = zeros(nbands,1);
                        % record the diference in all bands
                        rec_v_dif = zeros(length(IDs),nbands);
                        % normal fit qa = 0
                        qa = 0;
                        
                        for i_B=1:nbands
                            if isTide
                                [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                    autoTSFitWL(clrx(IDs),clry(IDs,i_B),update_num_c,tide(IDs));
                            else
                                 [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                    autoTSFit(clrx(IDs),clry(IDs,i_B),update_num_c);
                            end
                        end
                        
                        % updating information for the first iteration
                        % record time of curve start
                        rec_cg(num_fc).t_start = clrx(i_start);
                        % record time of curve end
                        rec_cg(num_fc).t_end = clrx(i);
                        % record break time
                        rec_cg(num_fc).t_break = 0; % no break at the moment
                        % record postion of the pixel
                        rec_cg(num_fc).pos = pos;
                        % record fitted coefficients
                        rec_cg(num_fc).coefs = fit_cft;
                        % record rmse of the pixel
                        rec_cg(num_fc).rmse = rmse;
                        % record change probability
                        rec_cg(num_fc).change_prob = 0;
                        % record number of observations
                        rec_cg(num_fc).num_obs = i-i_start+1;
                        % record fit category
                        rec_cg(num_fc).category = qa + update_num_c;
                        % record change magnitude
                        rec_cg(num_fc).magnitude = zeros(1,nbands);
		                % record durchange
		                rec_cg(num_fc).durchange = zeros(2, nbands);

                        if recordSR
                            if isTide
                                data = [clrx,clry,tide];
                                data = data(i_start_segment:i,:);
                                if ~isempty(bandNames)
                                    rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                                else
                                    rec_cg(num_fc).data = data;
                                end
                            else
                                data = [clrx,clry];
                                data = data(i_start_segment:i,:);
                                if ~isempty(bandNames)
                                    rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                                else
                                    rec_cg(num_fc).data = data;
                                end
                            end
                        end

                        % detect change
                        % value of difference for conse obs
                        v_dif = zeros(conse,nbands);
                        % record the magnitude of change
                        v_dif_mag = v_dif;
                        v_dif_clrx = zeros(conse,1);
                        
                        vec_mag = zeros(conse,1);
                        
                        for i_conse = 1:conse
                            for i_B = 1:nbands
                                % absolute difference
                                if isTide
                                    v_dif_mag(i_conse,i_B) = clry(i+i_conse,i_B)-autoTSPredWL(clrx(i+i_conse),fit_cft(:,i_B),tide(i+i_conse));
                                else
                                    v_dif_mag(i_conse,i_B) = clry(i+i_conse,i_B)-autoTSPred(clrx(i+i_conse),fit_cft(:,i_B));
                               
                                end
                                % normalized to z-scores
                                if sum(i_B == B_detect)
                                    % minimum rmse
                                    mini_rmse = max(adj_rmse(i_B),rmse(i_B));
                                    
                                    % z-scores
                                    v_dif(i_conse,i_B) = v_dif_mag(i_conse,i_B)/(mini_rmse+eps);
                                end
                            end
                            v_dif_clrx(i_conse) = clrx(i+i_conse);
                            vec_mag(i_conse) = norm(v_dif(i_conse,B_detect))^2;
                        end
                        % IDs that haven't updated
                        IDsOld = IDs;
                    else
                        if i - pre_end > num_skip % time to update the time series model
%                         if i - pre_end >= 3% % every 3 observations, update a model, version 13.04 recalled
                            % update i_count at each interation
                            i_count = clrx(i)-clrx(i_start);
                            
                            % record previous end i
                            pre_end = i;
                            % label as 0, indicating for being not ready to skip next observations
                            num_skip = 0; 
                        
                           % defining computed variables
                            fit_cft = zeros(max_num_c,nbands);
                            % rmse for each band
                            rmse = zeros(nbands,1);
                            % record the diference in all bands
                            rec_v_dif = zeros(length(IDs),nbands);
                            % normal fit qa = 0
                            qa = 0;
                                                       
%                             % moving window with nyr = 5 start
%                             nyr = 5;
%                             id_nyr = find(clrx(IDs(end))-clrx(IDs) < nyr*num_yrs);
%                             if ~isempty(id_nyr) && length(id_nyr) >= max_num_c*n_times
%                                 IDs = IDs(id_nyr(1):end);
%                             end
%                             % end of moving window
                                                       
                            for i_B = 1:nbands
                                if isTide
                                    [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                        autoTSFitWL(clrx(IDs),clry(IDs,i_B),update_num_c,tide(IDs));
                                else
                                    [fit_cft(:,i_B),rmse(i_B),rec_v_dif(:,i_B)] = ...
                                        autoTSFit(clrx(IDs),clry(IDs,i_B),update_num_c);

                                end
                            end
                            
                            % record fitted coefficients
                            rec_cg(num_fc).coefs = fit_cft;
                            % record rmse of the pixel
                            rec_cg(num_fc).rmse = rmse;
                            % record number of observations
                            rec_cg(num_fc).num_obs = i-i_start+1;
                            % record fit category
                            rec_cg(num_fc).category = qa + update_num_c;

                            % IDs that haven't updated
                            IDsOld = IDs;
                        end
                        
                        % record time of curve end
                        rec_cg(num_fc).t_end = clrx(i);
                        
                        if recordSR
                            if isTide
                                data = [clrx,clry,tide];
                                data = data(i_start_segment:i,:);
                                if ~isempty(bandNames)
                                    rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                                else
                                    rec_cg(num_fc).data = data;
                                end
                            else
                                data = [clrx,clry];
                                data = data(i_start_segment:i,:);
                                if ~isempty(bandNames)
                                    rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                                else
                                    rec_cg(num_fc).data = data;
                                end
                            end

                        end

                        % use fixed number for RMSE computing
                        n_rmse = n_times*rec_cg(num_fc).category;
                        tmpcg_rmse = zeros(nbands,1);
                        % better days counting for RMSE calculating
                        % relative days distance
                        d_rt = clrx(IDsOld) - clrx(i+conse);
                        d_yr = abs(round(d_rt/num_yrs)*num_yrs-d_rt);
                        
                        [~,sorted_indx] = sort(d_yr);
                        sorted_indx = sorted_indx(1:n_rmse);
                        
                        for i_B = B_detect
                            % temporally changing RMSE
                            tmpcg_rmse(i_B) = norm(rec_v_dif(IDsOld(sorted_indx)-IDsOld(1)+1,i_B))/...
                                sqrt(n_rmse-rec_cg(num_fc).category);
                        end
                        
                        % move the ith col to i-1th col
                        v_dif(1:conse-1,:) = v_dif(2:conse,:);
                        % only compute the difference of last consecutive obs
                        v_dif(conse,:) = 0;
                        % move the ith col to i-1th col
                        v_dif_mag(1:conse-1,:) = v_dif_mag(2:conse,:);
                        % record the magnitude of change of the last conse obs
                        v_dif_mag(conse,:) = 0;
                        % record the during change
                        v_dif_clrx(1:conse-1) = v_dif_clrx(2:conse);
                        % move the ith col to i-1th col
                        vec_mag(1:conse-1) = vec_mag(2:conse);
                        % change vector magnitude
                        vec_mag(conse) = 0;
                        
                        for i_B = 1:nbands
                            % absolute difference
                            if isTide
                                v_dif_mag(conse,i_B) = clry(i+conse,i_B)-autoTSPredWL(clrx(i+conse),fit_cft(:,i_B),tide(i+conse));
                            else
                                v_dif_mag(conse,i_B) = clry(i+conse,i_B)-autoTSPred(clrx(i+conse),fit_cft(:,i_B));
                           
                            end
                                % normalized to z-scores
                            if sum(i_B == B_detect)
                                % minimum rmse
                                mini_rmse = max(adj_rmse(i_B),tmpcg_rmse(i_B));
                                
                                % z-scores
                                v_dif(conse,i_B) = v_dif_mag(conse,i_B)/(mini_rmse+eps);
                            end
                        end
                        v_dif_clrx(conse) = clrx(i+conse);
                        vec_mag(conse) = norm(v_dif(end,B_detect))^2;
                    end
                    
                    % sign of change vector
                    % vec_sign = abs(sum(sign(v_dif(:,B_detect)),1));
                    max_angl = mean(angl(v_dif(:,B_detect)));
                    
                    % change detection
                    if min(vec_mag) > T_cg && max_angl < nsign% change detected
                        % record break time
                        rec_cg(num_fc).t_break = clrx(i+1);
                        % record change probability
                        rec_cg(num_fc).change_prob = 1;
                        % record change magnitude
                        rec_cg(num_fc).magnitude = median(v_dif_mag,1);
                         
                        % update during change
                        rec_cg(num_fc).durchange = durchange(nbands, v_dif_clrx, v_dif_mag);
                        

                        % identified and move on for the next functional curve
                        num_fc = num_fc + 1;
                        % start from i+1 for the next functional curve
                        i_start = i + 1;
                        % start training again
                        BL_train = 0;
                        
                        if recordSR
                            % Initilization of for the clear observation recording
                            i_start_segment = i_start;
                        end

                    elseif vec_mag(1) > Tmax_cg % false change
                        % remove noise
                        clrx(i+1) = [];
                        clry(i+1,:) = [];
                        if isTide
                            tide(i+1) = [];
                        end
                        i=i-1; % stay & check again after noise removal
                    end
                end % end of continuous monitoring
            end % end of checking basic requrirements
            
            % move forward to the i+1th clear observation
            i=i+1;
        end % end of while iterative
        
        % Two ways for processing the end of the time series
        if BL_train == 1
            % 1) if no break find at the end of the time series
            % define probability of change based on conse
            for i_conse = conse:-1:1
                % sign of change vector
                % vec_sign = abs(sum(sign(v_dif(i_conse:conse,B_detect)),1));
                max_angl = mean(angl(v_dif(i_conse:conse,B_detect)));
                
                if vec_mag(i_conse) <= T_cg || max_angl >= nsign
                    % the last stable id
                    id_last = i_conse;
                    break;
                end
            end
            
            % update change probability
            rec_cg(num_fc).change_prob = (conse-id_last)/conse;
            % update end time of the curve
            rec_cg(num_fc).t_end=clrx(end-conse+id_last);
            
            if recordSR   
                if isTide
                    data = [clrx,clry,tide];
                    data = data(i_start_segment:end,:);
                    if ~isempty(bandNames)
                        rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                    else
                        rec_cg(num_fc).data = data;
                    end
                else
                    data = [clrx,clry];
                    data = data(i_start_segment:end,:);
                    if ~isempty(bandNames)
                        rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                    else
                        rec_cg(num_fc).data = data;
                    end
                end
            end
            

            if conse > id_last % > 1
                % update time of the probable change
                rec_cg(num_fc).t_break = clrx(end-conse+id_last+1);
                % update magnitude of change
                rec_cg(num_fc).magnitude = median(v_dif_mag(id_last+1:conse,:),1);
                % update during change
                rec_cg(num_fc).durchange = durchange(nbands, v_dif_clrx, v_dif_mag);
            end
            

        elseif BL_train == 0
            % 2) if break find close to the end of the time series
            % Use [conse,min_num_c*n_times+conse) to fit curve
            
            if num_fc == rec_fc
                % first curve
                i_start = 1;
            else
                i_start = find(clrx >= rec_cg(num_fc-1).t_break);
                i_start = i_start(1);
            end
            
            % Tmask
            if length(clrx(i_start:end)) > conse
                blIDs = autoTmask(clrx(i_start:end),clry(i_start:end,[numBandGreen,numBandSWIR1]),...
                    (clrx(end)-clrx(i_start))/num_yrs,adj_rmse(numBandGreen),adj_rmse(numBandSWIR1),T_const);
                
                % update i_span after noise removal
                i_span = sum(~blIDs); %#ok<NASGU>
                
                IDs = i_start:length(clrx); % all IDs
                rmIDs = IDs(blIDs(1:end-conse) == 1); % IDs to be removed
                
                % remove noise pixels between i_start & i
                clrx(rmIDs) = [];
                clry(rmIDs,:) = [];
                if isTide
                    tide(rmIDs) = [];
                end
            end
            
            % enough data
            if length(clrx(i_start:end)) >= conse
                % defining computed variables
                fit_cft = zeros(max_num_c,nbands);
                % rmse for each band
                rmse = zeros(nbands,1);
                % end of fit qa = 20
                qa = 20;
                
                update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c);
                for i_B = 1:nbands
                    if isTide
                        [fit_cft(:,i_B),rmse(i_B)] = ...
                            autoTSFitWL(clrx(i_start:end),clry(i_start:end,i_B),update_num_c,tide(i_start:end));
                    else
                        [fit_cft(:,i_B),rmse(i_B)] = ...
                            autoTSFit(clrx(i_start:end),clry(i_start:end,i_B),update_num_c);
                    end
                end
                
                % record time of curve start
                rec_cg(num_fc).t_start = clrx(i_start);
                % record time of curve end
                rec_cg(num_fc).t_end=clrx(end);
                % record break time
                rec_cg(num_fc).t_break = 0;
                % record postion of the pixel
                rec_cg(num_fc).pos = pos;
                % record fitted coefficients
                rec_cg(num_fc).coefs = fit_cft;
                % record rmse of the pixel
                rec_cg(num_fc).rmse = rmse;
                % record change probability
                rec_cg(num_fc).change_prob = 0;
                % record number of observations
                rec_cg(num_fc).num_obs = length(clrx(i_start:end));
                % record fit category
                rec_cg(num_fc).category = qa + min_num_c;
                % record change magnitude
                rec_cg(num_fc).magnitude = zeros(1,nbands);
                % record durchange
                rec_cg(num_fc).durchange = zeros(2, nbands);

                if recordSR
                    if isTide
                        data = [clrx,clry,tide];
                        data = data(i_start:end,:);
                        if ~isempty(bandNames)
                            rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames,'tide']);
                        else
                            rec_cg(num_fc).data = data;
                        end
                    else
                        data = [clrx,clry];
                        data = data(i_start:end,:);
                        if ~isempty(bandNames)
                            rec_cg(num_fc).data = array2table(data,'VariableNames',['sdate',bandNames]);
                        else
                            rec_cg(num_fc).data = data;
                        end
                    end
                end
            end
        end
    end % end of if sum(idgood) statement    
   
end % end of function

% FUNCTION of updating the number of coeffs of the time series model
function update_num_c = update_cft(i_span,n_times,min_num_c,mid_num_c,max_num_c,num_c)

% determine the time series model
if i_span < mid_num_c*n_times
    % start with 5 coefficients model
    update_num_c = min(min_num_c,num_c);
elseif i_span < max_num_c*n_times
    % start with 7 coefficients model
    update_num_c = min(mid_num_c,num_c);
else
    % start with 9 coefficients model
    update_num_c =  min(max_num_c,num_c);
end

end

% FUNCTION of caculating included angle between ajacent pair of change vectors
function y = angl(v_dif)

[row,~] = size(v_dif);
y = zeros(row-1,1);

if row > 1
    for i = 1:row-1
        a = v_dif(i,:);
        b = v_dif(i+1,:);
        % y measures the opposite of cos(angle)
        y(i) = acos(a*b'/(norm(a)*norm(b)));
    end
else
    y = 0;
end

% convert angle from radiance to degree
y = y*180/pi;

end

% FUNCION of calculting dur-change trend with OLS fit
function v_dif_linear = durchange(nbands, clrx, v_dif_mag)
    v_dif_linear = zeros(2,nbands);
    for i_B = 1:nbands
        p = polyfit(clrx, v_dif_mag(:,i_B),1);
        v_dif_linear(1, i_B) = p(1); % yfit =  p(1) * x + p(2)
        v_dif_pred = p(1) .* clrx + p(2);
        v_dif_linear(2, i_B) = sqrt(mean((v_dif_mag(:,i_B) - v_dif_pred).^2));  % Root Mean Squared Error
    end
end
