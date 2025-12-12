function dist = calculate_accel_distance(accel)
%CALCULATE_ACCEL_DISTANCE Calculate horizontal distance covered during acceleration
%
% Level acceleration - all velocity contributes to horizontal distance

if isfield(accel, 'x') && ~isempty(accel.x)
    dist = accel.x(end);
elseif isfield(accel, 't') && isfield(accel, 'V')
    t = accel.t;
    V = accel.V;
    n = length(t);
    dist = 0;

    for i = 2:n
        dt = t(i) - t(i-1);
        V_avg = (V(i) + V(i-1)) / 2;
        dist = dist + V_avg * dt;
    end
else
    % Estimate: average of initial and final velocity
    V_avg = (accel.V(1) + accel.V(end)) / 2;
    dist = V_avg * accel.total_time;
end

end
