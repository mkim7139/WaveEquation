clear all;
close all;
clf;
clc;

% calculate weights for initial conditions
syms x;
weights=[];
for i=[1: 10]
    weights = [weights, int(0.2*x*sin(i*x), x, 0, pi/2)];
end

% generate data
x_domain = [0:0.1:2*pi];

% create movie struct
frames = 100;
movieVector(frames) = struct('cdata',[],'colormap',[]);

% draw figure and capture video
for t=[1:frames]
    clf
    f = HeatSoln(x_domain, t/8, weights, 0.3);

    % control axis size
    hold on;
    xlim([0 2*pi])
    ylim([-0.7 0.7])
    plot(x_domain, f, color='k');

    % save frame to video
    movieVector(t) = getframe;
end

% write video
myWrighter = VideoWriter('guitar_string', 'MPEG-4');
myWrighter.FrameRate = 30;
open(myWrighter);
writeVideo(myWrighter, movieVector);
close(myWrighter);

% fourier functions
function [f_0, f] = SuperSinSum(x_domain, weights)
    f_0 = zeros(1, length(x_domain));
    for i=[1:length(weights)]
        f_0 = f_0 + weights(i) * sin(i*x_domain);
    end
end

function [f] = HeatSoln(x_domain, t, weights, damp_coeff)
    f = zeros(1, length(x_domain));
    for i=[1:length(weights)]
        f = f + weights(i) * exp(-damp_coeff*t) * cos(i*t) * sin(i*x_domain);
    end
end