%% The proposed method (Our12): Route-guided time-expanded graph construction
clc;clear;close all;

%% Load data and basic setting
load('input_our2.mat'); 
adj = double(adj);
dis = double(dis);
zc = .95;
zd = .95;
yc=0.97;
yd = 0.97;
W = 1;  % kWh
T = 300; % Unit：s
num_period = 20; % simulation time
inoutcap = 1000; % same as demand
storage_cap = 1000*T*2; %store 2 periods worth of demand

%% Calculation
% 1. base
tic
[supplyG1,forward_star_rep1, look_up1] = specify_G_6(supply2, veh_routes2, zc, zd, W, T, yc, yd, inoutcap, storage_cap);
[optimal_flow1, cost1] = find_opt_flow_via_simplex_alg(supplyG1, forward_star_rep1);
t1 = toc;
n_arcs1 = size(forward_star_rep1,1);
n_nodes1 = size(look_up1,1);

% 2. our2
tic
[supplyG2,forward_star_rep2, look_up2] = specify_G_7(supply2, veh_routes2, zc, zd, W, T, yc, yd, inoutcap, storage_cap);
[optimal_flow2, cost2] = find_opt_flow_via_simplex_alg(supplyG2, forward_star_rep2);
t2 = toc;
n_arcs2 = size(forward_star_rep2,1);
n_nodes2 = size(look_up2,1);

% 3. our 12
n_trans = 1;
p_trans = 0.2;
tic
[veh_routes3,supply3] = get_routes_extend(supply2,n_trans,veh_routes2,p_trans,num_period);
[supplyG3,forward_star_rep3, look_up3] = specify_G_7(supply3, veh_routes3, zc, zd, W, T, yc, yd, inoutcap, storage_cap);
[optimal_flow3, cost3] = find_opt_flow_via_simplex_alg(supplyG3, forward_star_rep3);
t3 = toc;
n_arcs3 = size(forward_star_rep3,1);
n_nodes3 = size(look_up3,1);

fprintf('============ Base ============\n');
fprintf('cost: %.2f\n', cost1);
fprintf('time: %.5f\n', t1);
fprintf('number of arcs: %.2f\n', n_arcs1);
fprintf('number of nodes: %.5f\n', n_nodes1);

fprintf('============ Our2 ============\n');
fprintf('cost: %.2f\n', cost2);
fprintf('time: %.5f\n', t2);
fprintf('number of arcs: %.2f\n', n_arcs2);
fprintf('number of nodes: %.5f\n', n_nodes2);
error_absolute = abs(cost1-cost2)/cost1;
fprintf('Error between base and proposed (our2): %.2f\n', error_absolute);

fprintf('============ Our12 ============\n');
fprintf('cost: %.2f\n', cost3);
fprintf('time: %.5f\n', t3);
fprintf('number of arcs: %.2f\n', n_arcs3);
fprintf('number of nodes: %.5f\n', n_nodes3);
error_absolute = abs(cost1-cost3)/cost1;
fprintf('Error between base and proposed (our12): %.2f\n', error_absolute);


%% ======== Function ======== 
function [veh_routes,supply] = get_routes_extend(supply,n_extend,veh_routes,p_exchange,num_period)
    nodes_demand = find(supply(end,:)<0);
    nodes_supply = find(supply(end,:)>0);
    n_routes = size(veh_routes,2);  
    idx_routes = [];  
    idx_nodes = [nodes_demand,nodes_supply]; 
    for t = 1:num_period
        veh_routes_t = cell(2,n_routes);
        veh_routes_t(1,:) = veh_routes(1,:);
        veh_routes_t(2,:) = cellfun(@(x) x(t), veh_routes(2, :), 'UniformOutput', false);
        for e = 1:n_extend
            selected_routes_idx = select_routes(veh_routes_t,nodes_supply,nodes_demand);
            selected_nodes_idx = select_nodes(selected_routes_idx,veh_routes_t,p_exchange);   
            idx_routes = [idx_routes,selected_routes_idx];
            idx_nodes = [idx_nodes,selected_nodes_idx];
            nodes_supply = selected_nodes_idx;
        end
    end
    idx_routes = unique(idx_routes);
    idx_nodes = unique(idx_nodes);
    supply = supply(:,idx_nodes);
    veh_routes = veh_routes(:, idx_routes);
    for i = 1:size(veh_routes, 2)
        original_route = veh_routes{1, i};  
        original_T = veh_routes{3,i};
        idx_i = ismember(original_route, idx_nodes);
        new_route = original_route(idx_i);
        new_T = original_T(idx_i);
        veh_routes{1, i} = new_route; 
        veh_routes{3, i} = new_T; 
    end
    veh_routes(1, :) = cellfun(@(route) arrayfun(@(x) find(idx_nodes == x, 1), route), veh_routes(1, :), 'UniformOutput', false); 
end

function selected_nodes_idx = select_nodes(selected_routes_idx,veh_routes,p_exchange)
    selected_routes = veh_routes(1,selected_routes_idx);  
    selected_flows = cell2mat(veh_routes(2,selected_routes_idx));
    all_nodes = cell2mat(selected_routes(:)'); 
    all_flows = repelem(selected_flows, cellfun(@numel, selected_routes));  
    [unique_nodes, ~, idx] = unique(all_nodes);  
    total_flows = accumarray(idx, all_flows); 
    [~, sort_idx] = sort(total_flows, 'descend');
    important_nodes_by_flow = unique_nodes(sort_idx);
    num_top_nodes = ceil(p_exchange * length(important_nodes_by_flow));
    selected_nodes_idx = important_nodes_by_flow(1:num_top_nodes);
end

function selected_routes_idx = select_routes(veh_routes,nodes_supply,nodes_demand)
    routes = veh_routes(1, :);
    selected_routes_idx = find(cellfun(@(route) check_order(route, nodes_supply, nodes_demand), routes));   
end

function is_valid = check_order(route, nodes_supply, nodes_demand)
    supply_idx = find(ismember(route, nodes_supply), 1, 'first');
    demand_idx = find(ismember(route, nodes_demand), 1, 'first');
    is_valid = ~isempty(supply_idx) && ~isempty(demand_idx) && (supply_idx < demand_idx);
end

