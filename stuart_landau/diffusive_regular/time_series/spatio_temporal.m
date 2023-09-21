% Sample data
clear all;
x=load('L_omega=2.0_ts.out');
x1=x(:,2);
x2=x(:,3);
x3=x(:,4);
% Combine the three time series into a single 3D array
all_time_series = [x1,x2,x3];

% Create the spatiotemporal plot using imagesc
imagesc(x(:,1), 3:-1:1, all_time_series');

% Set axis labels and title
%xlabel('Time');
%ylabel('Time Series');
%title('Spatiotemporal Plot of Three Time Series');

cmap=colormap((inferno));  
colorbar;
yticks(1:3); % Set ticks at positions 1, 2, and 3
yticklabels({'3', '2', '1'}); % Set custom tick labels