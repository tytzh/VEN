
function [supplyG ,forward_star_rep, look_up,junctions_nodes_times] = specify_G_9(num_extra,p_extra_demand,t_start,t_end,supply, veh_routes, zc, zd, W, T, yc, yd, inoutcap, storage_cap)
R = size(veh_routes,2); % R = number of routes
n_D = size (supply,2); %number of road junctions
num_period = t_end - t_start + 1;
windowt = [t_start:t_end];
visited_R_T = cell(R ,num_period); 
visited_C = cell(R ,num_period); 

%% 
t_diff = t_start-1;
 num_period_all = min(size(supply,1),t_end+num_extra) ;
for t_start = windowt
    for r = 1:R       
        nodes = veh_routes{1, r};    
        rel_times = veh_routes{3, r}; 
        abs_times = t_start + rel_times - 1; 
        flows = veh_routes{2, r}; 

        abs_times_clean = abs_times(abs_times <= num_period_all);
        nodes_clean = nodes (abs_times <= num_period_all);

        capi = arrayfun(@(j) W*T*sum(flows(abs_times_clean(j):abs_times_clean(j+1))), 1:length(abs_times_clean)-1);

        temp =[nodes_clean ;abs_times_clean ];
        visited_R_T{r,t_start - t_diff} =  temp;
        visited_C{r,t_start - t_diff} = capi;
    end
end

%% 1.
look_up = zeros(0,5);
forward_star_rep = zeros(0,5);
l_up_used = 0;
for r = 1:R
    for t = 1:num_period
        visit_data = visited_R_T{r, t};

        if isempty(visit_data)
            continue;
        end

        nodes = visit_data(1, :);
        times = visit_data(2, :);

        for i = 1:length(nodes)
            entry = [nodes(i), 0, 0, times(i)];

            if ~ismember(entry, look_up(:,1:4), 'rows')
                look_up(end+1, :) = [entry, l_up_used+1];
                l_up_used = l_up_used + 1;
            end
        end
    end
end
n_junctions = l_up_used;
f_star_used = 0;  


%% 2.
B = inf;
for r = 1:R
    for t = 1:num_period
        R_Ti = visited_R_T{r,t};
        routei = R_Ti(1,:);
        Ti = R_Ti(2,:);
        Ci = visited_C{r,t};

        for i = 1:length(routei)-1
            i1 = i;
            i2 = i+1;
            nodes = routei(i1);
            nodet = routei(i2);
            times = Ti(i1);
            timet = Ti(i2);
           
            if i == 1
                look_up(l_up_used+1,:)= [nodes,r,i1,times,l_up_used+1];
                look_up(l_up_used+2,:)= [nodet,r,i2,timet,l_up_used+2];
                l_up_used = l_up_used  + 2;
            else
                look_up(l_up_used+1,:) = [nodet,r,i2,timet,l_up_used+1];
                l_up_used = l_up_used  + 1;
            end

            idx_nodes_xuni = look_up(all(look_up(:,1:4)==[nodes,r,i1,times],2),5);
            idx_nodet_xuni = look_up(all(look_up(:,1:4)==[nodet,r,i2,timet],2),5);           
            idx_nodes_zhenshi = look_up(all(look_up(:,1:4)==[nodes,0,0,times],2),5);
            idx_nodet_zhenshi = look_up(all(look_up(:,1:4)==[nodet,0,0,timet],2),5);
            
            % 1) nodes
            if i1 == 1
                forward_star_rep(f_star_used+1,:) = [idx_nodes_zhenshi,idx_nodes_xuni,1-zc,B,zc];
                f_star_used  = f_star_used+1;
            else
                forward_star_rep(f_star_used+1,:) =  [idx_nodes_zhenshi,idx_nodes_xuni,1-zc,B,zc];
                forward_star_rep(f_star_used+2,:) =  [idx_nodes_xuni,idx_nodes_zhenshi,1-zc,B,zc];

                f_star_used  = f_star_used+2;
            end

            % 2) nodet
            if i2 == length(routei)
                forward_star_rep(f_star_used+1,:) =  [idx_nodet_xuni,idx_nodet_zhenshi,1-zc,B,zc];
                f_star_used  = f_star_used+1;
            else
                forward_star_rep(f_star_used+1,:) =  [idx_nodet_zhenshi,idx_nodet_xuni,1-zc,B,zc];
                forward_star_rep(f_star_used+2,:) = [idx_nodet_xuni,idx_nodet_zhenshi,1-zc,B,zc];
                f_star_used  = f_star_used+2;
            end

            % 3) nodes --> nodet
            capi = Ci(i);
            forward_star_rep(f_star_used+1,:) = [idx_nodes_xuni,idx_nodet_xuni,0,capi,1];
            f_star_used  = f_star_used+1;

        end
    end
end


%% for arcs for storage
for t = windowt
for i = 1:n_D

    %create storage node
    look_up(l_up_used+1,:) = [i, -1, -1, t,l_up_used+1] ; 

    node_at_t = look_up(all(look_up(:,1:4)==[i,0,0,t],2),5);  % the number assinged to junction node i at t
    if isempty(node_at_t)
        continue;
    end

    if t==windowt(1)
        forward_star_rep(f_star_used+1,:) = [node_at_t l_up_used+1 1-yc inoutcap*T yc];
        
    elseif t == windowt(end)
        forward_star_rep(f_star_used+1,:) = [l_up_used+1 node_at_t 1-yd inoutcap*T yd];
        forward_star_rep(f_star_used+2,:) = [num_of_storage_at_prev_iteration l_up_used+1 0 storage_cap 1];
        forward_star_rep(f_star_used+3,:) = [l_up_used+1 l_up_used+1 0 inf 0.5]; %self-loop arc to allow energy to be left in storage
    
        f_star_used = f_star_used+3;
    else
         
        forward_star_rep(f_star_used+1,:) = [node_at_t l_up_used+1 1-yc inoutcap*T yc];
        forward_star_rep(f_star_used+2,:) = [l_up_used+1 node_at_t 1-yd inoutcap*T yd];
        forward_star_rep(f_star_used+3,:) = [num_of_storage_at_prev_iteration l_up_used+1 0 storage_cap 1];
        
        f_star_used = f_star_used+3;
    end
    num_of_storage_at_prev_iteration = l_up_used+1;
    l_up_used = l_up_used+1; % we have assigned one node in this iteration
    
end
end

%% 3. supplyG
supplyG = zeros(1,size(look_up,1));
times = look_up(1:n_junctions, 4);  
nodes = look_up(1:n_junctions, 1);
junctions_nodes_times = [nodes ,times];
supplyG(1:n_junctions) = T * supply(sub2ind(size(supply), times, nodes))';
mask_time = look_up(1:n_junctions, 4) > t_end;
mask_negative = supplyG(1:n_junctions) < 0;
supplyG(mask_time'&mask_negative) = supplyG(mask_time'&mask_negative).*p_extra_demand;  


end


