%%
%
clc
close all
clearvars
%
%% To load the data set of eckart compressors
% digitalised from (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
%
load("eckardt_impeller_O.mat")
%
data_O = data_impeller_O;
%
load("eckardt_impeller_A.mat")
%
data_A = data_impeller_A;
%
load("eckardt_impeller_B.mat")
%
data_B = data_impeller_B;
%
%% To set the operating conditions
%
p_01 = 1.01 * 1e5;
T_01 = 288.15;
%
fluid = 'air';
%
cmp = centrifugal_compressor();
cmp = cmp.set_inlet_conditions(p_01,T_01,fluid);
%
%% Impeller O
%
cmp_O = cmp.set_geometry_eckardt_impeller_O(); %It sets the Eckardt geometry for impeller O
%
for i = 1 : height(data_O)
    %
    cmp_O = cmp_O.set_operating_conditions(data_O.mdot(i),data_O.rpm(i));
    %
    cmp_O = cmp_O.simulation_stage();
    %
    pr_tt_model_O(i) = cmp_O.PR_tt_stage;
    eta_is_tt_model_O(i) = cmp_O.eta_is_tt_stage;
    %
end
%
rmse_pr_O = rmse(pr_tt_model_O',data_O.pr);
rmse_eta_O = rmse(eta_is_tt_model_O',data_O.eta_is_tt);
%
%% Impeller A
%
cmp_A = cmp.set_geometry_eckardt_impeller_A(); %It sets the Eckardt geometry for impeller A
%
for i = 1 : height(data_A)
    %
    cmp_A = cmp_A.set_operating_conditions(data_A.mdot(i),data_A.rpm(i));
    %
    cmp_A = cmp_A.simulation_stage();
    %
    pr_tt_model_A(i) = cmp_A.PR_tt_stage;
    eta_is_tt_model_A(i) = cmp_A.eta_is_tt_stage;
    %
end
%
rmse_pr_A = rmse(pr_tt_model_A',data_A.pr);
rmse_eta_A = rmse(eta_is_tt_model_A',data_A.eta_is_tt);
%
%% Impeller B
%
cmp_B = cmp.set_geometry_eckardt_impeller_B(); %It sets the Eckardt geometry for impeller B
%
for i = 1 : height(data_B)
    %
    cmp_B = cmp_B.set_operating_conditions(data_B.mdot(i),data_B.rpm(i));
    %
    cmp_B = cmp_B.simulation_stage();
    %
    pr_tt_model_B(i) = cmp_B.PR_tt_stage;
    eta_is_tt_model_B(i) = cmp_B.eta_is_tt_stage;
    %
end
%
rmse_pr_B = rmse(pr_tt_model_B',data_B.pr);
rmse_eta_B = rmse(eta_is_tt_model_B',data_B.eta_is_tt);
%
%% Figures
% It requires "brewermap" package for the plot colors (It is just a
% graphical thing)
% download it from: https://it.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
%
colors = brewermap(4,'set2');
red = brewermap(9,'reds');
red = red(7,:);
blue = brewermap(9,'blues');
blue = blue(7,:);
green = brewermap(9,'greens');
green = green(7,:);
%
figure
t1 = tiledlayout(2,4,"TileSpacing","compact","Padding","compact");
%
% Impeller O
%
nexttile(t1,1)
%
for i = 1 : 4
    %
    idx2plot = 1 + (10 * (i-1)):10 * (1 + i-1);
    %
    plot(data_O.mdot(idx2plot),data_O.pr(idx2plot),...
        "Marker","o","Color",colors(i,:),"LineWidth",2,"LineStyle","none")
    %
    hold on
    %
    plot(data_O.mdot(idx2plot),pr_tt_model_O(idx2plot),...
    "Color",colors(i,:),"LineWidth",2,"LineStyle","-")
    %
end
%
grid on
xlim([2 8])
xticks([2 4 6 8])
ylim([1 2.8])
yticks([1 1.5 2 2.5])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$\dot{m}\;(kg/s)$","Interpreter","latex")
ylabel("$PR_{tt}\;(-)$","Interpreter","latex")
text(0.05,0.95,"(a)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
%
nexttile(t1,5)
%
for i = 1 : 4
    %
    idx2plot = 1 + (10 * (i-1)):10 * (1 + i-1);
    %
    plot(data_O.mdot(idx2plot),data_O.eta_is_tt(idx2plot),...
        "Marker","o","Color",colors(i,:),"LineWidth",2,"LineStyle","none")
    %
    hold on
    %
    plot(data_O.mdot(idx2plot),eta_is_tt_model_O(idx2plot),...
    "Color",colors(i,:),"LineWidth",2,"LineStyle","-")
    %
end
%
grid on
ylim([0.6 0.9])
yticks([0.5 0.6 0.7 0.8 0.9])
xlim([2 8])
xticks([2 4 6 8])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$\dot{m}\;(kg/s)$","Interpreter","latex")
ylabel("$\eta_{tt}\;(-)$","Interpreter","latex")
text(0.05,0.1,"(e)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
%
dummy1 = plot(nan,nan, ...
    "Marker","o","Color",'k',"LineWidth",2,"LineStyle","none");
dummy2 = plot(nan,nan, ...
    "Color",'k',"LineWidth",2,"LineStyle","-");
dummy3 = plot(nan,nan, ...
    "Color",colors(1,:),"LineWidth",10,"LineStyle","-");
dummy4 = plot(nan,nan, ...
    "Color",colors(2,:),"LineWidth",10,"LineStyle","-");
dummy5 = plot(nan,nan, ...
    "Color",colors(3,:),"LineWidth",10,"LineStyle","-");
dummy6 = plot(nan,nan, ...
    "Color",colors(4,:),"LineWidth",10,"LineStyle","-");
legend([dummy1,dummy2,dummy3,dummy4,dummy5,dummy6], ...
    ["Data","Model","10 kRPM","12 kRPM","14 kRPM","16 kRPM"], ...
    "FontSize",20,'Location','southeast','Interpreter','latex')
%
% Impeller A
%
nexttile(t1,2)
%
for i = 1 : 4
    %
    idx2plot = 1 + (5 * (i-1)):5 * (1 + i-1);
    %
    plot(data_A.mdot(idx2plot),data_A.pr(idx2plot),...
        "Marker","o","Color",colors(i,:),"LineWidth",2,"LineStyle","none")
    %
    hold on
    %
    plot(data_A.mdot(idx2plot),pr_tt_model_A(idx2plot),...
    "Color",colors(i,:),"LineWidth",2,"LineStyle","-")
    %
end
%
grid on
xlim([2 8])
xticks([2 4 6 8])
ylim([1 2.8])
yticks([1 1.5 2 2.5])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$\dot{m}\;(kg/s)$","Interpreter","latex")
ylabel("$PR_{tt}\;(-)$","Interpreter","latex")
text(0.05,0.95,"(b)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
%
nexttile(t1,6)
%
for i = 1 : 4
    %
    idx2plot = 1 + (5 * (i-1)):5 * (1 + i-1);
    %
    plot(data_A.mdot(idx2plot),data_A.eta_is_tt(idx2plot),...
        "Marker","o","Color",colors(i,:),"LineWidth",2,"LineStyle","none")
    %
    hold on
    %
    plot(data_A.mdot(idx2plot),eta_is_tt_model_A(idx2plot),...
    "Color",colors(i,:),"LineWidth",2,"LineStyle","-")
    %
end
%
grid on
ylim([0.6 0.9])
yticks([0.5 0.6 0.7 0.8 0.9])
xticks([2 4 6 8])
xlim([2 8])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$\dot{m}\;(kg/s)$","Interpreter","latex")
ylabel("$\eta_{tt}\;(-)$","Interpreter","latex")
% title("(f)")
text(0.05,0.1,"(f)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
%
% Impeller B
%
nexttile(t1,3)
%
for i = 1 : 4
    %
    idx2plot = 1 + (5 * (i-1)):5 * (1 + i-1);
    %
    plot(data_B.mdot(idx2plot),data_B.pr(idx2plot),...
        "Marker","o","Color",colors(i,:),"LineWidth",2,"LineStyle","none")
    %
    hold on
    %
    plot(data_B.mdot(idx2plot),pr_tt_model_B(idx2plot),...
    "Color",colors(i,:),"LineWidth",2,"LineStyle","-")
    %
end
%
grid on
xlim([2 8])
xticks([2 4 6 8])
ylim([1 2.8])
yticks([1 1.5 2 2.5])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$\dot{m}\;(kg/s)$","Interpreter","latex")
ylabel("$PR_{tt}\;(-)$","Interpreter","latex")
text(0.05,0.95,"(c)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
%
nexttile(t1,7)
%
for i = 1 : 4
    %
    idx2plot = 1 + (5 * (i-1)):5 * (1 + i-1);
    %
    plot(data_B.mdot(idx2plot),data_B.eta_is_tt(idx2plot),...
        "Marker","o","Color",colors(i,:),"LineWidth",2,"LineStyle","none")
    %
    hold on
    %
    plot(data_B.mdot(idx2plot),eta_is_tt_model_B(idx2plot),...
    "Color",colors(i,:),"LineWidth",2,"LineStyle","-")
    %
end
%
grid on
ylim([0.6 0.9])
yticks([0.5 0.6 0.7 0.8 0.9])
xticks([2 4 6 8])
xlim([2 8])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$\dot{m}\;(kg/s)$","Interpreter","latex")
ylabel("$\eta_{tt}\;(-)$","Interpreter","latex")
text(0.05,0.1,"(g)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
%
ax = nexttile(t1,4,[1 1]);
%
x = [data_O.pr; data_A.pr; data_B.pr]; 
y = [pr_tt_model_O pr_tt_model_A pr_tt_model_B];
%
plot(sort(x),sort(x),...
    "Color","k","LineWidth",2)
%
hold on
%
plot( ...
    sort(x),0.95 * sort(x),...
    sort(x),1.05 * sort(x),...
    "Color","k","LineWidth",2,"LineStyle","--")
%
plot( ...
    sort(x),0.90 * sort(x),...
    sort(x),1.10 * sort(x),...
    "Color","k","LineWidth",2,"LineStyle","-.")
%
p4 = plot(data_O.pr,pr_tt_model_O,...
    "Marker","s","Color",red,"LineWidth",3,"LineStyle","none","MarkerSize",9);
%
p5 = plot(data_A.pr,pr_tt_model_A,...
    "Marker","s","Color",blue,"LineWidth",3,"LineStyle","none","MarkerSize",9);
%
p6 = plot(data_B.pr,pr_tt_model_B,...
    "Marker","s","Color",green,"LineWidth",3,"LineStyle","none","MarkerSize",9);
%
grid on
xticks([1 1.5 2 2.5])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1)
xlabel("$PR_{tt,data}\;(-)$","Interpreter","latex")
ylabel("$PR_{tt,model}\;(-)$","Interpreter","latex")
% title("(d)")
text(0.85,0.95,"(d)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
legend(...
    [p4 p5 p6],...
    ["Impeller O","Impeller A","Impeller B"],...
    'location','northwest','Interpreter','latex','NumColumns',1,'FontSize',18)

%
nexttile(t1,8,[1 1])
%
x = [data_O.eta_is_tt; data_A.eta_is_tt; data_B.eta_is_tt]; 
y = [eta_is_tt_model_O eta_is_tt_model_A eta_is_tt_model_B];
%
p1 = plot(sort(x),sort(x),...
    "Color","k","LineWidth",2);
%
hold on
%
p2 = plot(sort(x),0.95 * sort(x), ...
    "Color","k","LineWidth",2,"LineStyle","--");
plot(sort(x),1.05 * sort(x),...
    "Color","k","LineWidth",2,"LineStyle","--");
p3 = plot(sort(x),0.90 * sort(x),...
    "Color","k","LineWidth",2,"LineStyle","-.");
plot(sort(x),1.10 * sort(x),...
    "Color","k","LineWidth",2,"LineStyle","-.");
%
p4 = plot(data_O.eta_is_tt,eta_is_tt_model_O,...
    "Marker","s","Color",red,"LineWidth",3,"LineStyle","none","MarkerSize",9);
%
p5 = plot(data_A.eta_is_tt,eta_is_tt_model_A,...
    "Marker","s","Color",blue,"LineWidth",3,"LineStyle","none","MarkerSize",9);
%
p6 = plot(data_B.eta_is_tt,eta_is_tt_model_B,...
    "Marker","s","Color",green,"LineWidth",3,"LineStyle","none","MarkerSize",9);
%
grid on
%
ylim([0.5 1])
yticks([0.5 0.6 0.7 0.8 0.9 1])
xticks([0.5 0.6 0.7 0.8 0.9])
xlim([0.6 0.9])
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",1,"XTickLabelRotation",0)
xlabel("$\eta_{tt,data}\;(-)$","Interpreter","latex")
ylabel("$\eta_{tt,model}\;(-)$","Interpreter","latex")
text(0.05,0.95,"(h)","Units",'normalized',"FontSize",28,"FontWeight","bold","FontName","Times New Roman")
legend(...
    [p1 p2 p3],...
    ["Ideal","$\pm 5\%$ error","$\pm 10\%$ error"],...
    'location','southeast','Interpreter','latex','NumColumns',1,'FontSize',18)
%
%% Error characterisation
% It calculates the model error bands compared to experimental data
%
figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
%
nexttile
%
err_pr = ([pr_tt_model_O'; pr_tt_model_A'; pr_tt_model_B'] - [data_O.pr; data_A.pr; data_B.pr]) ./ [data_O.pr; data_A.pr; data_B.pr];
err_band_pr = prctile(err_pr,[5 95]);% The error band covers from the 5th to the 95th percentiles of the experimental data set
%
histogram(err_pr)
hold on
xline(err_band_pr(1),"LineWidth",2,"LineStyle","--")
xline(err_band_pr(2),"LineWidth",2,"LineStyle","--")
%
grid on
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",2)
xlabel("PR_{tt} Relative error (-)")
ylabel("Bin count (-)")
legend("Error Histogram","5^{th} and 95^{th} error percentiles")
%
nexttile
%
err_eta = ([eta_is_tt_model_O'; eta_is_tt_model_A'; eta_is_tt_model_B'] - [data_O.eta_is_tt; data_A.eta_is_tt; data_B.eta_is_tt]) ./ [data_O.eta_is_tt; data_A.eta_is_tt; data_B.eta_is_tt];
err_band_eta = prctile(err_eta,[5 95]);
%
histogram(err_eta)
hold on
xline(err_band_eta(1),"LineWidth",2,"LineStyle","--")
xline(err_band_eta(2),"LineWidth",2,"LineStyle","--")%
%
grid on
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",2)
xlabel("\eta_{is,tt} Relative error (-)")
ylabel("Bin count (-)")
% legend("Error Histogram","5^{th} and 95^{th} error percentiles")
%
%