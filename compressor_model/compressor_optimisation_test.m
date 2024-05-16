%%
%
clc
close all
clearvars
%
%%
%
% Inlet conditions, Operating conditions & Fluid
%
fluid = 'Air';
%
p_01 = 101325;  %(Pa)
T_01 = 95 + 273.15; %(K)
%
mdot = 0.72;    %(kg/s)
omega_rpm = 40 * 1e3;   %(RPM)
%
% Optimisation targets
%
PR_ts_target = 2.2; %(-) PR_ts target for single and multi objective optimisation 
eta_ts_target = 0.85;   %(-) eta_ts target for multi objective optimisation 
%
%% Case 1 - Single Objective Optimisation
%
tic
cmp_soo = centrifugal_compressor();  %It creates the compressor object
t_1 = toc;
%
tic
cmp_soo = cmp_soo.set_inlet_conditions(p_01,T_01,fluid);  %It sets the inlet conditions
%
cmp_soo = cmp_soo.set_operating_conditions(mdot,omega_rpm);   %It sets the operating conditions
t_2 = toc;
%
tic
% It optimizes the geometry to maximise the total-to-static isentropic
% efficiency with the specified operating conditions and total-to-static PR
% equal to the specified target
cmp_soo = cmp_soo.optimisation_geometry_singleObjective(PR_ts_target);    
t_3 = toc;
%
eta_is_tt_so = cmp_soo.eta_is_tt;
eta_is_ts_so = cmp_soo.eta_is_ts;
%
PR_tt_so = cmp_soo.PR_tt;
PR_ts_so = cmp_soo.PR_ts;
PR_ss_so = cmp_soo.PR_ss;
%
%% Case 2 - Multiobjective optimisation
%
tic
cmp_moo = centrifugal_compressor();
t_1 = toc;
%
tic
cmp_moo = cmp_moo.set_inlet_conditions(p_01,T_01,fluid);
%
cmp_moo = cmp_moo.set_operating_conditions(mdot,omega_rpm);
t_2 = toc;
%
tic
% It optimizes the geometry to simultaneosly maximise the total-to-static isentropic
% efficiency and the total-to-static PR with the specified operating conditions
cmp_moo = cmp_moo.optimisation_geometry_multiObjective(PR_ts_target,eta_ts_target);
t_3 = toc;
%
eta_is_tt_moo = cmp_moo.eta_is_tt;
eta_is_ts_moo = cmp_moo.eta_is_ts;
%
PR_tt_moo = cmp_moo.PR_tt;
PR_ts_moo = cmp_moo.PR_ts;
PR_ss_moo = cmp_moo.PR_ss;
%
%% Figures
%
[impeller_profile_soo,vaneless_diffuser_profile_soo] = cmp_soo.plot_geometry();  %Rather than plotting the geometry profiles it generates them
[impeller_profile_moo,vaneless_diffuser_profile_moo] = cmp_moo.plot_geometry();
%
figure
tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
%
nexttile
cmp_soo.plot_geometry() %Rather than generating the geometry profiles it plots them
%
nexttile
cmp_moo.plot_geometry() 
%





