function draw_annulus_ring(vals, colors, r_in, r_out)
% draw_annulus_ring  Draw a two-slice annular (donut) ring as patches.
%   vals   : [v1, v2] nonnegative values (WUI, Wildland)
%   colors : 2x3 RGB rows in [0,1]
%   r_in   : inner radius
%   r_out  : outer radius

    vals = double(vals);
    vals(~isfinite(vals)) = 0;
    if all(vals <= 0)
        return
    end

    fracs = vals / sum(vals);
    angs  = 2*pi * [0, cumsum(fracs)];   % cumulative angles
    theta0 = pi/2;                       % start at 12 o'clock (upwards)

    for k = 1:numel(fracs)
        a0 = theta0 - angs(k);
        a1 = theta0 - angs(k+1);
        th = linspace(a0, a1, 80);

        % outer arc
        xo = r_out * cos(th);  yo = r_out * sin(th);
        % inner arc (reverse to close the patch)
        xi = r_in  * cos(fliplr(th)); yi = r_in  * sin(fliplr(th));

        patch([xo xi], [yo yi], colors(k,:), ...
              'EdgeColor','w','LineWidth',1);
    end
end
