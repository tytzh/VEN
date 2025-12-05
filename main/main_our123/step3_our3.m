%% Proposed method (our123)
clc;clear;close all;

%% 1. Load and setting
load('prediction_64.mat');load('input_our2.mat');
ns = 200;
T = 100; num_period = 800;
n_0 = 30;
zc = .95; zd = .95; yc=0.97; yd = 0.97;
inoutcap = 1000; storage_cap = 600000; W = 1;
gamma = 1; n_trans = 3; p_trans = 0.5;

%% 2. Processing
[supply2,veh_routes2] = get_vehroute_supply(windows,veh_routes2,gamma,ns,flow_lower,supply_lower,n_0,full_pred1_flow, full_pred1_supply);
[idx_routes,idx_nodes] = get_extend_idx(supply2,n_trans,veh_routes2,p_trans,num_period);
[supply2,veh_routes2] = update_route_supply(idx_routes,idx_nodes,supply2, veh_routes2);

%% 3. Planning
nw = size(windows,1);
result_mpc_robust = cell(nw,5); 
w1 = windows(1,1); w2 = windows(1,2);  
w3 = windows(2,1);
n_extra = 0;  p_extra_demand = 0; 
[supplyG_t,forward_star_rep_t, look_up_t,visited_R_T] = specify_G_9(n_extra,p_extra_demand,w1, w2,supply2, veh_routes2, zc, zd, W, T, yc, yd, inoutcap, storage_cap);
[optimal_flow_t, cost_t, Aeq_t, beq_t, lb_t, ub_t, f_t, ~] = find_opt_flow_via_simplex_alg2(supplyG_t, forward_star_rep_t,[]);
for i = 2:nw
    w0 = windows(i-1,1); w1 = windows(i-1,2); w2 = windows(i,1);  w3 = windows(i,2); 
    [edge_overlap_last,flow_overlap_last] = get_edge_flow_last_window(w2,look_up_t,forward_star_rep_t,optimal_flow_t);
    [supplyG_t,forward_star_rep_t, look_up_t,junctions_nodes_times] = specify_G_9(n_extra,p_extra_demand,w2, w3,supply2, veh_routes2, zc, zd, W, T, yc, yd, inoutcap, storage_cap);
    [supplyG_t_update,idx_matched] = update_supplyG_t(w1,supplyG_t,edge_overlap_last,flow_overlap_last,look_up_t,forward_star_rep_t,junctions_nodes_times);
    [optimal_flow_t, cost_t,~,~,~,~,f_t,~] = find_opt_flow_via_simplex_alg2(supplyG_t_update, forward_star_rep_t,idx_matched);
end


%% ======== Function ======== 
function [supply,veh_routes] = update_route_supply(idx_routes,idx_nodes,supply, veh_routes)
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

function [idx_routes,idx_nodes] = get_extend_idx(supply,n_extend,veh_routes,p_exchange,num_period)
    nodes_demand = find(supply(end,:)<-1e-3);
    nodes_supply = find(supply(end,:)>1e-3);
    n_routes = size(veh_routes,2);
    n_nodes = size(supply,2);
    idx_routes = [];
    idx_nodes = [nodes_demand,nodes_supply];
    for t = 1:num_period
        veh_routes_t = cell(2,n_routes);
        veh_routes_t(1,:) = veh_routes(1,:);
        veh_routes_t(2,:) = cellfun(@(x) x(t), veh_routes(2, :), 'UniformOutput', false);   
        for e = 1:n_extend
            selected_routes_idx = selct_routes(veh_routes_t,nodes_supply,nodes_demand);
            idx_routes = [idx_routes,selected_routes_idx];
            if length(idx_nodes)~=n_nodes
                selected_nodes_idx = select_nodes(selected_routes_idx,veh_routes_t,p_exchange,idx_nodes);
                idx_nodes = [idx_nodes,selected_nodes_idx];
                nodes_supply = selected_nodes_idx;
            end
        end
    end
    idx_routes = unique(idx_routes);
    idx_nodes = unique(idx_nodes);  
end

function [supply2,veh_routes2] = get_vehroute_supply(windows,veh_routes0,gamma,ns,flow_lower,supply_lower,n_0,full_pred1_flow, full_pred1_supply)    
    flow2 = full_pred1_flow{ns};
    supply2 = full_pred1_supply{ns};
    supply2 = supply2 + gamma*supply_lower;
    flow2 = flow2 + gamma*flow_lower;
    idx_demand = find(supply2(1,:)<0); 
    supply2(1:n_0,idx_demand) = 0; 
    flow2(flow2<0) = 0;
    ll = round((windows(1,2)-windows(2,1))*0.0);  
    idx_know = [];
    idx_missing = [];
    for n = 1:size(windows,1)-2
       w0 = windows(n+1,1);  
       w1 = windows(n,2);
       idx_missing_n = [w0:w1];
       idx_missing = [idx_missing,idx_missing_n];
    end
    idx_know = setdiff([1:size(flow2,1)],idx_missing);  
    for nsu = 1:size(supply2,2)
        temp = supply2(:,nsu);
        temp (idx_missing) = interp1(idx_know , temp(idx_know), idx_missing, 'pchip');
        temp(abs(temp)<1e-3)=0;
        supply2(:,nsu) = temp;
    end
    for nfl = 1:size(flow2,2)
        temp = flow2(:,nfl);
        temp (idx_missing) = interp1(idx_know , temp(idx_know), idx_missing, 'linear');
        flow2(:,nsu) = temp;   
    end
    n_route = size(veh_routes0,2);
    veh_routes2 = veh_routes0;
    for n = 1:n_route
        veh_routes2{2,n} = flow2(:,n);
    end
end

function  [supplyG_t_update,idx_matched] = update_supplyG_t(w1,supplyG_t,edge_overlap_last,flow_overlap_last,look_up_t,forward_star_rep_t,junctions_nodes_times)

    node1 = edge_overlap_last(:,1:4);
    node2 = edge_overlap_last(:,5:8);

    [~, idx_node1] = ismember(node1, look_up_t(:,1:4), 'rows'); 
    [~, idx_node2] = ismember(node2, look_up_t(:,1:4), 'rows');
    idx_node_st = [idx_node1,idx_node2];
    [is_match, ~] = ismember(idx_node_st, forward_star_rep_t(:,1:2), 'rows');
    idx_matched = find(is_match);
    
    is_target1 = node1(:,3) == 0 & node1(:,2) == 0;  
    is_target2 = node2(:,3) == 0 & node2(:,2) == 0;
    
    outflow = [node1(is_target1,1), node1(is_target1,4), -ones(sum(is_target1),1).*flow_overlap_last(is_target1)];
    inflow  = [node2(is_target2,1), node2(is_target2,4), ones(sum(is_target2),1).*flow_overlap_last(is_target2)];
    result_flow = [outflow; inflow];  
    [is_match, idx_in_target] = ismember(result_flow(:,1:2), junctions_nodes_times, 'rows');
    supply_accumulated = accumarray(idx_in_target(is_match), result_flow(is_match,3), ...
                              [size(junctions_nodes_times,1), 1], @sum, 0);
    supply_accumulated(end+1 : length(supplyG_t)) = 0;
    supplyG_t_update = supplyG_t + supply_accumulated';

    idx_le_w1 = find(look_up_t(:,4) <= w1);
    supplyG_t_update(idx_le_w1(supplyG_t_update(idx_le_w1) < 0)) = 0;

end

function [edge_overlap,flow_overlap]=get_edge_flow_last_window(w2,look_up_t,forward_star_rep_t,optimal_flow_t)
    idx_node_overlap = find( ...
        look_up_t(:,4) >= w2 & ...
        look_up_t(:,2) == 0 & ...
        look_up_t(:,3) == 0);
    idx_edge_overlap = find( ...
        ismember(forward_star_rep_t(:,1), idx_node_overlap) | ...
        ismember(forward_star_rep_t(:,2), idx_node_overlap));
    nodes_edges_overlap = forward_star_rep_t(idx_edge_overlap,1:2); 
    edge_overlap = [look_up_t(nodes_edges_overlap(:,1),1:4), look_up_t(nodes_edges_overlap(:,2),1:4)];  % 具体的edge里node的数据
    flow_overlap = optimal_flow_t(idx_edge_overlap,3);
end

function selected_nodes_idx = select_nodes(selected_routes_idx,veh_routes,p_exchange,idx_nodes)    
    selected_routes = veh_routes(1,selected_routes_idx);  
    selected_flows = cell2mat(veh_routes(2,selected_routes_idx));
    all_nodes = cell2mat(selected_routes(:)'); 
    all_flows = repelem(selected_flows, cellfun(@numel, selected_routes));  
    [unique_nodes, ~, idx] = unique(all_nodes);  
    total_flows = accumarray(idx, all_flows);
    [~, sort_idx] = sort(total_flows, 'descend');
    important_nodes_by_flow = unique_nodes(sort_idx);
    mask = ~ismember(important_nodes_by_flow, idx_nodes);
    filtered_nodes = important_nodes_by_flow(mask);
    num_top_nodes = min(length(filtered_nodes),ceil(p_exchange * length(important_nodes_by_flow)));
    selected_nodes_idx = filtered_nodes(1:num_top_nodes);
end

function selected_routes_idx = selct_routes(veh_routes,nodes_supply,nodes_demand)
    routes = veh_routes(1, :);
    selected_routes_idx = find(cellfun(@(route) check_order(route, nodes_supply, nodes_demand), routes)); 
end

function is_valid = check_order(route, nodes_supply, nodes_demand)
    supply_idx = find(ismember(route, nodes_supply), 1, 'first');
    demand_idx = find(ismember(route, nodes_demand), 1, 'first');
    is_valid = ~isempty(supply_idx) && ~isempty(demand_idx) && (supply_idx < demand_idx);
end

