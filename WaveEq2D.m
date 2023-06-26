clear all;
close all;
clc;
clf;


% generate data
x_domain = [0:0.03:1];
y_domain = [0:0.03:1];
t_domain = [0:0.1:10];
[X,Y] = meshgrid(x_domain, y_domain);
Z = HeatSoln(X, Y, t_domain);

% create video struct
frames = length(t_domain);
movieVector(frames) = struct('cdata',[],'colormap',[]);

% draw figure and save frames
fig = figure('Name', 'myfig');
for t=[1:frames]
    clf
    hold on;
%     fig.Position = [600 600 600 600];

    xlim([0 1])
    ylim([0 1])
    zlim([-0.2 0.2])
    view([30, 35])
    surf(X,Y,Z{t})
    movieVector(t) = getframe;
end

% write video
myWrighter = VideoWriter('drum', 'MPEG-4');
myWrighter.FrameRate = 10;
open(myWrighter);
writeVideo(myWrighter, movieVector);
close(myWrighter);


% fourier function
function [Z] = HeatSoln(X, Y, t_domain)
    % initialize data
    spacial_domain = size(X);
    time_domain = length(t_domain);
    Z = cell(1, time_domain);
%     Z = zeros(spacial_domain(1), spacial_domain(2));

    % weights and initial conditions
    syms x;
    syms y;
    m=5;
    n=5;

    % calculate weights for initial conditions
    % generate surface
    for t=[1:time_domain]
        % initialize frame before generating surface
        Z{t} = zeros(spacial_domain(1), spacial_domain(2));

        for i=[1:m]
            for j=[1:n]

                % position initial condition
                p_weight = (4/1) * int(int((0.1*sin(3*pi*x) + 0.3*sin(pi*x))*0.5*y * sin(i*pi*x)*sin(j*pi*y), y, 0, 1), x, 0, 1);

                % generate surface for frame t
                Z{t} = Z{t} + exp(-0.2*t_domain(t)) * double(p_weight*sin(i*pi*X).*sin(j*pi*Y) * cos(pi*sqrt((i/2)^2 + (j/3)^2)*t_domain(t)));
            end
        end
    end
end

