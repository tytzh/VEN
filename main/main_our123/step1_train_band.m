%% Step1: Training model and obtain banding
clc; clear; close all;

%% 1. Load data
load('input_our3.mat');

%% 2. Models
% 1. Dataset Construction
[XTrainSupply,YTrainSupply] = get_dataset_combined2(train_supply_norm, supply_avg_norm, windows);
[XTrainFlow,YTrainFlow] = get_dataset_combined2(train_flow_norm,flow_avg_norm, windows);

% 2. Network
numHiddenUnits = 64;
net_supply = Training_LSTM(XTrainSupply,YTrainSupply,numHiddenUnits);
net_flow = Training_LSTM(XTrainFlow,YTrainFlow,numHiddenUnits);

%% 3. Banding
[supply_upper,supply_lower] = get_robust_band2(list_supply,idx_train,net_supply,supply_avg_norm,train_supply_norm,windows,supply_avg,scaler_supply);
[flow_upper,flow_lower] = get_robust_band2(list_flow,idx_train,net_flow,flow_avg_norm,train_flow_norm,windows,flow_avg,scaler_flow);

%% Save training outcomes：
save('LSTM_model_64.mat');

%% ==== Function ====
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

function [supply_upper,supply_lower] = get_robust_band2(list_supply,idx_test,net_supply,supply_avg_norm,test_supply_norm,windows,supply_avg,scaler_supply)
    [~,pred_list]  = Testing_LSTM2(net_supply,supply_avg_norm,test_supply_norm,windows,scaler_supply);
    pred1 = pred_list(:,2);
    full_pred1_supply = stitch_windows(pred1, supply_avg, windows);
    full_true1_supply = list_supply(idx_test);
    [supply_upper,supply_lower] = get_upper_lower2(full_pred1_supply,full_true1_supply,windows);
end

function [error_upper,error_lower] = get_upper_lower2(full_pred1_supply,full_true1_supply,windows)
    n = length(full_pred1_supply);
    list_error = cell(n,1);
    for i = 1:n
        errori = full_true1_supply{i} - full_pred1_supply{i} ;
        errori(find(abs(errori)<1e-4)) = 0;
        list_error{i} =  errori;
    end
    error_upper = -Inf(size(errori));
    error_lower = +Inf(size(errori));
    for k = 1:n
        current_error = list_error{k};
        error_upper = max(error_upper, current_error);
        error_lower = min(error_lower, current_error);
    end
    error_upper(windows(1,1):windows(1,2),:) = 0;
    error_lower(windows(1,1):windows(1,2),:) = 0;

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

function net = Training_LSTM(XTrain,YTrain,numHiddenUnits)
    inputSize = size(XTrain{1},1);
    outputSize = size(YTrain{1},1);
    layers = [ ...
        sequenceInputLayer(inputSize)
        lstmLayer(numHiddenUnits)
        fullyConnectedLayer(outputSize)
        regressionLayer];
    options = trainingOptions('adam', ...
        'MaxEpochs', 50, ...
        'MiniBatchSize', 64, ...
        'Shuffle','every-epoch', ...
        'Verbose', false, ...
        'Plots','training-progress', ...
        'ExecutionEnvironment','auto');
    net = trainNetwork(XTrain, YTrain, layers, options);
end

function [XTrain, YTrain] = get_dataset_combined2(list_supply, supply_avg, windows)   
    XTrain = {};
    YTrain = {};
    for u = 1:length(list_supply)
        supply_real = list_supply{u};
        for w = 1:(size(windows,1) - 1)
            t1 = windows(w,1);   t2 = windows(w,2);
            t3 = windows(w+1,1); t4 = windows(w+1,2);
            if t2 <= size(supply_real,1) && t4 <= size(supply_real,1)
                input_avg  = supply_avg(t1:t2,  :)';  
                input_real = supply_real(t1:t2, :)'; 
                input_avg_next = supply_avg(t3:t4,  :)'; 
                XTrain{end+1} = [input_avg; input_real; input_avg_next]; 
                YTrain{end+1} = supply_real(t3:t4, :)';
            end
        end
    end
end
