%% The proposed method (Our1): Flow-Guided Graph Reduction
clc
clear
close all
rng('default')
rng_num = 0.8;  


%% Load data and basic setting
load('input_our1.mat');
zc = .95;  % efficiency of charging
zd = .95; 
W = 1; 
buffer_ratio = 0.8;
n_trans = 1;
p_trans = 0.2;
adj = double(adj);
dis = double(dis);
 

%% Calculation
tic
extend_veh_routes = get_routes_extend(supply,n_trans,veh_routes,p_trans);
[supplyG, forward_star_rep, ~, ~] = specify_G_1(supply, extend_veh_routes, zc, zd, W);
t2_graph = toc;
tic
[x_G_2, cost_2, Aeq2 ,beq2 , ~, ~, ~] = find_opt_flow_via_simplex_alg(supplyG, forward_star_rep);
t2_lp = toc;
n_nodes2 = length(supplyG);
n_arcs2 = size(forward_star_rep,1);


%% 3. Show Result：  
fprintf('============ Proposed(Our1) ============\n');
fprintf('cost: %.2f\n', cost_2);
fprintf('graph time: %.5f\n', t2_graph);
fprintf('LP time: %.5f\n', t2_lp);
fprintf('Total time: %.5f\n', t2_lp+t2_graph );
fprintf('number of arcs: %.2f\n', n_arcs2);
fprintf('number of nodes: %.5f\n', n_nodes2);




%% ======== Function ======== 
function veh_routes = get_routes_extend(supply,n_extend,veh_routes,p_exchange)

    nodes_demand = find(supply<0);
    nodes_supply = find(supply>0); 
    idx_routes = [];  
    idx_nodes = [nodes_demand,nodes_supply];
   
    for e = 1:n_extend
    
        selected_routes_idx = select_routes(veh_routes,nodes_supply,nodes_demand);
        selected_nodes_idx = select_nodes(selected_routes_idx,veh_routes,p_exchange);
    
        idx_routes = [idx_routes,selected_routes_idx];
        idx_nodes = [idx_nodes,selected_nodes_idx];
    
        nodes_supply = selected_nodes_idx;
    end
    idx_routes = unique(idx_routes);
    idx_nodes = unique(idx_nodes);
    
    veh_routes = veh_routes(:, idx_routes);
    veh_routes(1, :) = cellfun(@(route) route(ismember(route, idx_nodes)), veh_routes(1, :), 'UniformOutput', false);

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
    demand_idx = find(ismember(route, nodes_demand), 1, 'last');
    is_valid = ~isempty(supply_idx) && ~isempty(demand_idx) && (supply_idx < demand_idx);
end


