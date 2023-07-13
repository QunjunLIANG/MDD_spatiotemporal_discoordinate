% This function calculates the probability flow (PF) base on the TD and 
% Peak correlation according to Mitra et al (2020).
%
% The function uses the TD and peak correlation estimated generated
% in time delay computation by Ryan's code (https://github.com/ryraut/lag-code). 

% Qunjun Liang 2022/01/15 

function [probflow] = probability_flow_esti(time_delay,peak_corr)

	s = size(time_delay);
   
    % linearize
	time_delay = reshape(time_delay,[s(1)*s(2),1])';
    peak_corr = reshape(peak_corr,[s(1)*s(2),1])';
    
    % initial the output 
    probflow = zeros(1,s(1)*s(2));
    
    % identify the correlation under threshold and fitlered out missing
    % value
    thresh_ind = find(abs(peak_corr)> 0.1);
    nan_ind = find(~isnan(time_delay));
    use_ind = intersect(thresh_ind, nan_ind);
    
    % normalize the filtered data
    time_delay(use_ind) = zscore(time_delay(use_ind));
    peak_corr(use_ind) = zscore(peak_corr(use_ind));

    % find the positive tau value
    tau_pos_ind = find(time_delay >= 0);
    pro_flow_pos_ind = intersect(use_ind, tau_pos_ind);
    probflow(pro_flow_pos_ind) = sqrt(time_delay(pro_flow_pos_ind).^2 + peak_corr(pro_flow_pos_ind).^2);
    
    % find indication of the negative tau value
    tau_neg_ind = find(time_delay < 0);
    pro_flow_neg_ind = intersect(use_ind, tau_neg_ind);
    probflow(pro_flow_neg_ind) = -sqrt(time_delay(pro_flow_neg_ind).^2 + peak_corr(pro_flow_neg_ind).^2);
    
    % export the results
	probflow = reshape(probflow,[s(1) s(2)]);
    probflow = probflow - diag(diag(probflow)); % it should be zero to the probability flow matrix
    
	
end


