% This function normalizes the projection map of probability flow
% according to Mitra et al (2020).
%
% The input pf_proj_map should be a double, not a mtrix. 

% Qunjun Liang 2022/01/15 

function [pf_proj_map_out] = probability_flow_norm(pf_proj_map)

   % initial the output
   pf_proj_map_out = pf_proj_map;

   % find positive values
   pos_value_ind = find(pf_proj_map>0);
   pos_value = pf_proj_map(pos_value_ind);
   pos_value_norm = pos_value./nansum(pos_value);
   
   % find negative values
   neg_value_ind = find(pf_proj_map<0);
   neg_value = pf_proj_map(neg_value_ind);
   neg_value_norm = - neg_value./nansum(neg_value);
	
   % re-value the projection map
   pf_proj_map_out(pos_value_ind) = pos_value_norm;
   pf_proj_map_out(neg_value_ind) = neg_value_norm;
   
   
end


