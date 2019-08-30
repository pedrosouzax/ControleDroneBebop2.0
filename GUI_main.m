clear all
close all
clc


path(path,'.\Interface')
path(path,'.\trajectory')

%% Graphical User Interface
gui = guidata(GUI_interface);

%% GUI settings
% Buttons settings

set(gui.Plotbutton,'Callback','plot_gui_graphics');
set(gui.start_button,'Callback','Main');


% Axes settings
set(gcf,'CurrentAxes',gui.graph_xyz)
plot3(0,0,0)
dist = 3;
set(gui.graph_xyz,'Xlim',[-dist dist],'Ylim',[-dist dist],'Zlim',[0 dist])
hold on;
grid;
[quad1_obj] = load_object_f; % Object of the quadrotor in the Axes 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gcf,'CurrentAxes',gui.graph_xy)
dist = 3;
set(gui.graph_xy,'Xlim',[-dist dist],'Ylim',[-dist dist],'Zlim',[0 dist])
hold on;
grid;
[quad2_obj] = load_object_f; % Object of the quadrotor in the Axes 2



