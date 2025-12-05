function [supplyG, forward_star_rep, look_up, A] = specify_G_1(supply, veh_routes, zc, zd, W)
% this function specify_G_1 creates a graph G with supply and
% forward_star_rep that specifies the supplies of each node in graph G
% and arcs in graph G respectively

% Note that the self-loop arcs for nodes with positive supplies are created
% in find_opt_flow_via_simplex_alg, not here

% It takes as input 
% 1)the supply of road junction in the road network
% 2) veh_routes
    % veh_routes is a cell structure that contains the vehicular routes and
    % flows on each route. It is a 2 by |R| structure where |R| is the number
    % of vehicular routes. Each route is represented as a cell in the top row
    % of veh_routes, with the cell containing a vector of the list of nodes
    % that represents the path of the route. The flow on that route is recorded
    % as a scalar or vector (for time-varying VEN) in the cell in the same column as the route in veh_routes
% 3,4) zc, zd, which are the charging and discharging effciencies respectively
% 5)W = amount of energy each EV can carry



R = size(veh_routes,2); % R = number of routes
n_D = size (supply,2); %number of road junctions

%% look_up
% look_up is a matrix that stores the numerical node assignment of the
% non-junction nodes that are created for G
% Each non-junction node is a row in look_up.

% Each non-junction node has nomaclature J_m^(n).n.m where
% n = the n^th route
% J_m^n = the junction node that the non-junction node is connected
% to
% m = the number of node along the route (which is a path) n that
% J_m^n is

% look_up has 4 columns: 1=J_m^n, 2=n, 3=m, 4
% number assigned to the non-junction node

% the original junction nodes are also included in look up, their 4
% columns:
% 1. junction number 2=0, 3=0, 4=number assigned = junction number


look_up = zeros(n_D,4);
look_up(:,1) = 1:n_D;
look_up(:,4) = 1:n_D;

% calculate the totoal number of non-junction nodes 
sum = 0;
for n = 1:R %for the n^th route
    route = cell2mat(veh_routes(1,n));% the route considered in this iteration
    r_len = size(route,2);
    sum = sum+r_len;  
end
%pre-allocate entries in look_up
look_up = [look_up; zeros(sum,4)];
n_nodes = size(look_up,1);  
A =  zeros(n_nodes,n_nodes);

num_assigned = n_D; 
% at any time, num_assigned stores the number of nodes that have been assigned
% since node assignment starts from 1, the next assignment of node is always
% num_assigned + 1 
% now there are n_D nodes with number assigned (from 1 to num_assigned), which 
% are the road junction nodes
% so non-junction node starts with n_D+1

%% forward_star_rep
% forward_star_rep is the matrix that stores all the arcs in G
% five columns are tail, head, cost, cap, gain
forward_star_rep = [];
B = inf; %capacity in many arcs

for i = 1:R
    %disp(['-----------considering ' num2str(i) 'th route--------']);
    route = cell2mat(veh_routes(1,i));% the route considered in this iteration
    r_len = size(route,2); %length of the route, = no. of nodes visited
    flow_i = cell2mat(veh_routes(2,i));%flow on route i
    
    %create two nodes corr to start and end of route  
    supply = [supply 0 0]; 
    look_up(num_assigned+1,:) = [route(1,1) i 1  num_assigned+1];
    look_up(num_assigned+2,:) = [route(1,end) i r_len  num_assigned+2];
    %store num_assigned2last_discharging_node to use at the end of for loop
    num_assigned2last_discharging_node = num_assigned+2;
    
    A(route(1,1),num_assigned+1) = flow_i;
    A(num_assigned+2,route(1,end)) = flow_i;

    % let   num_assigned2charging_node_of_prev_it    denote
    % the node number assigned to teh charging node of the previous
    % iteration, to be used in for loop below
    num_assigned2nonjunction_node_of_prev_it = num_assigned+1;
    
    %add directed arc to forward_star_rep
    forward_star_rep = [forward_star_rep; ...
        route(1,1) num_assigned+1 1-zc B zc;...
        num_assigned+2 route(1,end) 1-zd B zd];
    
    % adjust num_assigned since we have added two nodes
    num_assigned = num_assigned +2; 
    
    for j = 2:r_len-1
        %disp(['---' num2str(j) 'th node along route----']);
        junc_node = route(j); %the intermediate node of route (junction node)
        %create 1 nodes corr to each intermediate node in route
        supply = [supply 0];
        look_up(num_assigned+1,:) = [junc_node i j num_assigned+1];
        
        %create two arcs
        forward_star_rep = [forward_star_rep; ...
            num_assigned+1 junc_node 1-zd B zd;...
            junc_node num_assigned+1 1-zc B zc];
        %create an arc to connect current node to node
        %from previous iteration
        forward_star_rep = [forward_star_rep; ...
            num_assigned2nonjunction_node_of_prev_it    num_assigned+1   0   W*flow_i   1];
        
        A(num_assigned+1,junc_node) = flow_i;  
        A(junc_node,num_assigned+1) = flow_i;
        A(num_assigned2nonjunction_node_of_prev_it,num_assigned+1) = flow_i; 


        %adjust num_assigned_to_last_charging_node for use in next iteration
        num_assigned2nonjunction_node_of_prev_it = num_assigned+1;
        
        % adjust num_assigned for use in next iteration since we have added two nodes
        num_assigned = num_assigned +1; 
          
    end
   
    %create one arc to link discharging node of the last node in route to 
    % the charging node of the second to last node in route
    
    forward_star_rep = [forward_star_rep; ...
        num_assigned2nonjunction_node_of_prev_it    num_assigned2last_discharging_node   0   W*flow_i   1];
    A(num_assigned2nonjunction_node_of_prev_it,num_assigned2last_discharging_node ) = flow_i;
    
end

%% supply 
supplyG = supply;

end

