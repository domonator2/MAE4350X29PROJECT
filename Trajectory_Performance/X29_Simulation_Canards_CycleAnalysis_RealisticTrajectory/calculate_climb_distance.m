function dist = calculate_climb_distance(climb)
%CALCULATE_CLIMB_DISTANCE Calculate horizontal distance covered during climb
%
% Uses climb profile velocity and time to integrate horizontal distance

if isfield(climb, 'x') && ~isempty(climb.x)
    % If distance is directly available
    dist = climb.x(end);
elseif isfield(climb, 't') && isfield(climb, 'V') && isfield(climb, 'h')
    % Calculate from trajectory
    t = climb.t;
    V = climb.V;
    h = climb.h;

    n = length(t);
    dist = 0;

    for i = 2:n
        dt = t(i) - t(i-1);
        dh = h(i) - h(i-1);
        V_avg = (V(i) + V(i-1)) / 2;

        % Horizontal distance from climb angle
        if V_avg > 0 && dt > 0
            sin_gamma = dh / (V_avg * dt);
            sin_gamma = max(-1, min(1, sin_gamma));
            cos_gamma = sqrt(1 - sin_gamma^2);
            dx = V_avg * cos_gamma * dt;
            dist = dist + dx;
        end
    end
else
    % Estimate based on climb profile
    % Average climb: 10 deg angle, M~0.6, 172s
    V_avg = 250;  % m/s average
    t_climb = climb.total_time;
    gamma_avg = 10;  % degrees
    dist = V_avg * cosd(gamma_avg) * t_climb;
end

end
