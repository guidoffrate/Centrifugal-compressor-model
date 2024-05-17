%%
%
clc
close all
clearvars
%
%% Geomtry specs
%
geometry.N_blades = 18;                         %(-)
geometry.blade_thickness_inlet = 0.3 * 1e-3;    %(m)
geometry.blade_thickness_outlet = 0.3 * 1e-3;   %(m)
geometry.roughness = 2 * 1e-5;                  %(m)
geometry.r2 = 0.15;                             %(m)
geometry.ratio_r1shroud_r2 = 0.5;               %(m) -  ratio between r1_hub and r2
geometry.ratio_r1hub_r1shroud = 0.53;           %(m) - ratio between r1_hub and r1_shroud
geometry.alpha1 = 0;                            %(deg)
geometry.beta1_blade_rms = -33;                 %(deg)
geometry.ratio_b2_r2 = 0.08;                    %(-) - ratio between b2 and r2
geometry.beta2_blade = -37;                     %(deg)
geometry.ratio_r3_r2 = 1.4;                     %(-) - ratio between r3 and r2
%
%% Inlet conditions
%
p_01 = 2 * 1e5;
p_env = 1 * 1e5;
T_01 = 150 + 273.15;
T_env = 15 + 273.15;
%
fluid = 'air';
%
%% Operating conditions
%
mdot = [1.5:0.1:3.5]; %They vary to plot a map
omega_rpm = [19:28] * 1e3;
%
%%
%
cmp = centrifugal_compressor();                     %it creates the object
cmp = cmp.set_geometry(geometry);                   %it sets the geometry
cmp = cmp.set_inlet_conditions(p_01,T_01,fluid);    %it sets the inlet conditions
%
for i = 1 : length(mdot)
    for j = 1 : length(omega_rpm)
        %
        cmp = cmp.set_operating_conditions(mdot(i),omega_rpm(j));
        cmp = cmp.simulation();
        %
        PR_ts(i,j) = cmp.PR_ts;
        eta_ts(i,j) = cmp.eta_is_ts;
        %
    end
end
%
%% Figures
%
mdot_c = mdot * sqrt(T_01/T_env) / (p_01/p_env);     %corrected mass flowrate
N_c = omega_rpm / sqrt(T_01/T_env);                 %corrected rotating speed
%
[mdot_c,N_c] = meshgrid(mdot_c,N_c);    %transformed into matrices for the plot
%
figure
hold on
[c,h] = contourf(mdot_c,PR_ts',eta_ts',...
    [0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.77 0.78 0.79 0.8 0.85],...
    "ShowText","on","LineWidth",2,"LineStyle","--","EdgeColor","k","FaceAlpha",0.5);
clabel(c,h,"FontSize",24,"FontName","Times New Roman","LabelSpacing",1200)
clim([min(eta_ts,[],"all") max(eta_ts,[],"all")])
%
[c,h] = contour(mdot_c',PR_ts,N_c',...
    [16 18 20 22 23 24]*1e3,...
    "ShowText","on","LineWidth",2,"LineStyle","-","EdgeColor","k");
clabel(c,h,"FontSize",24,"FontName","Times New Roman","LabelSpacing",1200)
%
grid on
box on
ylabel("$PR_{ts} (-)$","Interpreter","latex")
xlabel("$\dot{m}_{c}=\dot{m}\sqrt{T/T_{ref}}/(p/p_{ref})\;(-)$","Interpreter","latex")
set(gca,"FontName","Times New Roman","FontSize",28,"LineWidth",2)
legend("$\eta_{ts}$","$RPM/\sqrt{T/T_{ref}}$",'Interpreter','latex','Location','southwest')
%
%