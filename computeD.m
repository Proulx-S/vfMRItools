function diam = computeD(area)
    % Diameter proxy from an area proxy, assuming a circular cross-section:
    %   D = 2*sqrt(A/pi)
    % Elementwise; diam has the same size as area.
    diam = 2*sqrt(area/pi);
end
