function [best_delay best_paths min_sse] = optimal_delay_mediation_search(x, y, m, min_delay, max_delay)
    if(min_delay >= max_delay)
        error('Min must be less than max delay');
    end
    if(~isscalar(min_delay) || ~isscalar(max_delay))
        error('Delay values must be scalars');
    end
    if(mod(min_delay, 1) ~= 0 || mod(max_delay, 1) ~= 0)
        error('Delay values must be integers');
    end
    
    delay_range = min_delay:max_delay;
    min_sse = Inf;
    intcpt = ones(size(x));
    
    for i=1:length(delay_range)
        for j=1:length(delay_range)
            [total_sse paths] = mediation_shift_sse([delay_range(i) delay_range(j)], x, y, m, intcpt);
            if(total_sse < min_sse)
                min_sse = total_sse;
                best_paths = paths;
                best_delay = [delay_range(i) delay_range(j)];
            end
        end
    end
end