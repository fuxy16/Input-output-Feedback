function c = cost(P_K,sigma_0)
% calculate cost with model
c=trace(P_K*sigma_0);
end

