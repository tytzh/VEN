%% Step2: Get prediction
% This prediction step is used to speed up computation in Step 3. 
% It can also be integrated into Step 3 for joint prediction and decision-making.
clc; clear; close all;
load('LSTM_model_64.mat');

%% 2. Test
[full_pred1_supply,full_true1_supply] = get_full_predictions(list_supply,idx_test,net_supply,supply_avg_norm,test_supply_norm,windows,supply_avg,scaler_supply);
[full_pred1_flow,full_true1_flow] = get_full_predictions(list_flow,idx_test,net_flow,flow_avg_norm,test_flow_norm,windows,flow_avg,scaler_flow);
save('prediction_64.mat');

%% ==== Function ====
function [full_pred1_supply,full_true1_supply] = get_full_predictions(list_supply,idx_test,net_supply,supply_avg_norm,test_supply_norm,windows,supply_avg,scaler_supply)
    [~,pred_list]  = Testing_LSTM2(net_supply,supply_avg_norm,test_supply_norm,windows,scaler_supply);
    pred1 = pred_list(:,2); 
    full_pred1_supply = stitch_windows(pred1, supply_avg, windows);
    full_true1_supply =  list_supply(idx_test);
end

function full_data_list = stitch_windows(delta_list, supply_avg, windows)
    num_windows = size(windows, 1);
    num_samples = numel(delta_list{1});
    num_nodes = size(delta_list{1}{1}, 2);
    T_max = windows(end, 2);  
    full_data_list = cell(num_samples, 1);
    w0 = windows(1, 1); w1 = windows(1, 2);
    for s = 1:num_samples
        full_data_list{s} = nan(T_max, num_nodes);
        full_data_list{s}(w0:w1,:) = supply_avg(w0:w1,:);
    end
    for w = 2:num_windows
        t1 = windows(w, 1);
        t2 = windows(w, 2);
        len = t2 - t1 + 1;
        for s = 1:num_samples
            window_data = delta_list{w-1}{s};
            full_data_list{s}(t1:t2, :) = window_data; 
        end
    end
    for i = 1:num_samples
        temp = full_data_list{i};
        temp = process_windows(temp,windows);
         full_data_list{i} = temp;
    end
end

function [delta_list,pred_list] =  Testing_LSTM2(net_supply,supply_avg_norm,test_supply_norm,windows,scaler)
    num_samples = length(test_supply_norm);
    num_windows = size(windows, 1);
    delta_list = cell(num_windows-1,1); 
    pred_list = cell(num_windows-1,2); 
    for k = 1:num_windows-1
        temp = {}; 
        temp1 = {}; temp2 = {};
        t1 = windows(k,1);   t2 = windows(k,2);    
        t3 = windows(k+1,1); t4 = windows(k+1,2);  
        for i = 1:num_samples
            supply = test_supply_norm{i};
            input_avg  = supply_avg_norm(t1:t2,  :)'; 
            input_real = supply(t1:t2, :)'; 
            input_avg_next = supply_avg_norm(t3:t4,  :)'; 
            XTest = [input_avg; input_real; input_avg_next];
            YTrue = supply(t3:t4, :)';
            YPred = predict(net_supply, {XTest});
            YPred = YPred{1};
            YPred_stored = YPred' .* (scaler.max - scaler.min + 1e-8) + scaler.min;
            YTrue_stored = YTrue' .* (scaler.max - scaler.min + 1e-8) + scaler.min;
            delta_i = YTrue_stored - YPred_stored ; 
            temp{i} = delta_i;
            temp1{i} = YTrue_stored;
            temp2{i} = YPred_stored;
        end
        delta_list{k,1} = temp;
        pred_list{k,1} = temp1;
        pred_list{k,2} = temp2;
    end
end

function supply2 = process_windows(supply2,windows)
    ll = round((windows(1,2)-windows(2,1))*0.0); 
    idx_know = [];
    idx_missing = [];
    for n = 1:size(windows,1)-2
       w0 = windows(n+1,1);
       w1 = windows(n,2);
       idx_missing_n = [w0:w1];
       idx_missing = [idx_missing,idx_missing_n];
    end
    idx_know = setdiff([1:size(supply2,1)],idx_missing);
    for nsu = 1:size(supply2,2)
        temp = supply2(:,nsu);
        temp (idx_missing) = interp1(idx_know , temp(idx_know), idx_missing, 'pchip'); 
        supply2(:,nsu) = temp;
    end
end
