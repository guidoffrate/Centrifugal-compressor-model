classdef centrifugal_compressor
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        %
        % General properties
        %
        N_blades                %(-) - Number of impeller blades
        % blade_thickness       %(m) - Blade thickness
        blade_thickness_inlet   %(m) - Blade thickness at the impeller inlet
        blade_thickness_outlet  %(m) - Blade thickness at the impeller outlet
        roughness               %(m) - surface roughness
        disk_clearance          %(m) - Disk clearance
        radial_clearance        %(m) - Radial clearance
        axial_length            %(m) - Impeller axial length
        %
        mdot                    %(kg/s) - Compressed mass flow rate
        %
        omega_rpm               %(RPM) - Rotating speed
        omega_rads              %(rad/s) - Rotating speed
        %
        convergence             % structure containing convergence residuals for each compressor section
        optimisation            % structure containing optimisation results
        fsolve_options          % structure containing the options for fsolve
        
        %
        % Fluid properties
        %
        fluid                   % operating fluid
        fluid_as                % fluid abstract state generated with coolprop
        %
        cp_fluid                %(kJ/kg/K) - Gas cp (considered constant along the compressor)
        k_fluid                 %(-) - Gas cp/cv (considered constant along the compressor)
        R_fluid                 %(J/K) - Gas constant
        mu_ref                  %(Pa s) reference viscosity @(0 + 273.15 K, 1*1e5 Pa) for Sutherland's law
        %
        % Section 1 (impeller inlet)
        %
        r1_shroud               %(m) - Radius at the shroud
        r1_hub                  %(m) - Radius at the hub
        r1_rms                  %(m) - Radius at the root mean square radius
        alpha1                  %(deg) - Angle of c
        beta1_blade_rms         %(deg) - Blade angle at the rms
        beta1_blade_hub         %(deg) - Blade angle at the hub
        beta1_blade_shroud      %(deg) - Blade angle at the shroud
        beta1_shroud            %(deg) - Relative speed angle at the shroud
        beta1_hub               %(deg) - Relative speed angle at the hub
        beta1_rms               %(deg) - Relative speed angle at the root mean square radius
        %
        A1                      %(m2) - Flow area
        %
        T_1                     %(K) - Static temperature
        p_1                     %(bar) - Static pressure
        h_1                     %(kJ/kg) - Static enthalpy
        s_1                     %(kJ/kg/K) - Entropy
        rho_1                   %(kg/m3) - Density
        %
        T_01                    %(K) - Total temperature
        p_01                    %(bar) - Total pressure
        h_01                    %(bar) - Total enthalpy
        s_01                    %(kJ/kg/K) - Entropy at p_01 and T_01                
        rho_01                  %(kg/m3) - Density at p_01 and T_01
        %
        Ma_1                    %(-) - Mach number
        %
        u1_rms                  %(m/s) - tangential speed at the root mean square radius                  
        u1_hub                  %(m/s) - tangential speed at the hub
        u1_shroud               %(m/s) - tangential speed at the shroud 
        c1_meridional           %(m/s) - absolute speed meridional component
        c1_tangential           %(m/s) - absolute speed tangential component
        c1                      %(m/s) - absolute speed
        w1_rms                  %(m/s) - relative speed at the root mean square radius
        w1_shroud               %(m/s) - relative speed at shroud
        w1_hub                  %(m/s) - relative speed at the hub
        %
        % Section 1 throat (minimum flow area near to Section 1)
        %
        A1_throat               %(m2) - Flow area
        %
        T_1throat               %(K) - Static temperature
        p_1throat               %(bar) - Static pressure
        rho_1throat             %(kg/m3) - Density
        %
        Ma_1throat              %(-) - Mach number
        %
        w1_throat               %(m/s) - relative speed
        %
        % Section 2 (impeller outlet)
        %
        r2                      %(m) - Radius at the impeller outlet
        b2                      %(m) - Blade height at the impeller outlet
        alpha2                  %(deg) - Absolute speed angle
        beta2                   %(deg) - Relative speed angle
        beta2_blade             %(deg) - Blade angle             
        %
        A2                      %(m2) - Flow area
        %
        T_2                     %(K) - Static temperature
        T_2is                   %(K) - Static temperature at the impeller's isentropic outlet conditions
        p_2                     %(bar) - Static pressure
        h_2                     %(kJ/kg) - Static enthalpy
        h_2is                   %(kJ/kg) - Static enthalpy at the impeller's isentropic outlet conditions
        s_2                     %(kJ/kg/K) - Entropy
        rho_2                   %(kg/m3) - Density
        %
        T_02                    %(K) - Total temperature
        T_02is                  %(K) - Total temperature at the impeller's isentropic outlet conditions
        p_02                    %(bar) - Total pressure
        p_02is                  %(bar) - Total pressure at the impeller's isentropic outlet conditions
        h_02                    %(kJ/kg) - Total enthalpy
        h_02is                  %(kJ/kg) - Total enthalpy at the impeller's isentropic outlet conditions
        %
        Ma_2                    %(-) - Mach number
        %
        u2                      %(m/s) - tangential speed at impeller outlet
        c2_tangential           %(m/s) - absolute speed tangential component
        c2_meridional           %(m/s) - absolute speed meridional component
        c2                      %(m/s) - absolute speed
        w2_tangential           %(m/s) - relative speed tangential component
        w2                      %(m/s) - relative speed
        %
        % Section 3 (vaneless diffuser outlet)
        %
        r3                      %(m) - Vaneless diffuser outlet radius
        b3                      %(m) - Vaneless diffuser height
        %
        A3                      %(m2) - Flow area
        %
        T_3                     %(K) - Static temperature
        T_3is                   %(K) - Static temperature at the vaneless diffuser's isentropic outlet conditions
        p_3                     %(bar) - Static pressure
        h_3                     %(kJ/kg) - Static enthalpy
        h_3is                   %(kJ/kg) - Static enthalpy at the vaneless diffuser's isentropic outlet conditions
        s_3                     %(kJ/kg/K) - Entropy
        rho_3                   %(kg/m3) - Density
        %
        T_03                    %(K) - Total temperature
        T_03is                  %(K) - Total temperature at the vaneless diffuser's isentropic outlet conditions
        p_03                    %(bar) - Total pressure
        p_03is                  %(bar) - Total pressure at the vaneless diffuser's isentropic outlet conditions
        h_03                    %(kJ/kg) - Total enthalpy
        h_03is                  %(kJ/kg) - Total enthalpy at the vaneless diffuser's isentropic outlet conditions
        rho_3is                 %(kJ/kg) - Density at the vaneless diffuser's isentropic outlet conditions
        %
        Ma_3                    %(-) - Mach number
        Re_3                    %(-) - Reynolds number
        %
        c3_tangential           %(m/s) - absolute speed tangential component
        c3_meridional           %(m/s) - absolute speed meridional component
        c3                      %(m/s) - absolute speed
        %
        % Section 4 (volute)
        %
        r4                      %(m) - Volute external radius
        %
        A4                      %(m2) - Flow area at the volute outlet (assumed as a circle or diameter r4 - r3)                      
        %
        T_4                     %(K) - Static temperature
        p_4                     %(bar) - Static pressure 
        h_4                     %(kJ/kg) - Static enthalpy
        s_4                     %(kJ/kg/K) - Entropy
        rho_4                   %(kg/m3) - Density
        %
        T_04                    %(K) - Total temperature
        p_04                    %(bar) - Total pressure
        h_04                    %(kJ/kg) - Total enthalpy
        %
        Ma_4                    %(-) - Mach number
        %
        c4                      %(m/s) - absolute speed
        %
        % Section 5 (Exhaust cone diffuser)
        %
        cone_length             %(m) - Exhaust cone's length
        cone_diameter_inlet     %(m) - Exhaust cone's inlet diameter (same as volute outlet diameter)
        cone_diameter_outlet    %(m) - Exhaust cone's outlet diameter
        cone_divergence_angle   %(deg) - Exhaust cone's wall angle with respect of the horizonal direction
        %
        A5                      %(m2) - Flow area at the cone outlet (assumed as a circle or diameter cone_diameter_outlet)
        %
        T_5                     %(K) - Static temperature
        p_5                     %(bar) - Static pressure 
        h_5                     %(kJ/kg) - Static enthalpy
        s_5                     %(kJ/kg/K) - Entropy
        rho_5                   %(kg/m3) - Density
        %
        T_05                    %(K) - Total temperature
        p_05                    %(bar) - Total pressure
        h_05                    %(kJ/kg) - Total enthalpy
        %
        Ma_5                    %(-) - Mach number
        %
        c5                      %(m/s) - absolute speed
        %
        % Stage quantities ("Stage" = "impeller" + "vaneless diffuser")
        %
        eta_is_tt_stage         %(-) - total-to-total isentropic efficiency
        PR_tt_stage             %(-) - total-to-total pressure ratio
        PR_ts_stage             %(-) - total-to-static pressure ratio
        PR_ss_stage             %(-) - static-to-static pressure ratio
        %
        % Global quantities ("global" = "impeller" + "vaneless diffuser" +
        % "Volute" + "cone")
        %
        eta_is_tt               %(-) - total-to-total isentropic efficiency
        eta_is_ts               %(-) - total-to-static isentropic efficiency
        PR_tt                   %(-) - total-to-total pressure ratio
        PR_ts                   %(-) - total-to-static pressure ratio
        PR_ss                   %(-) - static-to-static pressure ratio
        %
        L_eulero                %(kJ/kg) - Euler's work
        Dh_internal_loss        %(kJ/kg) - Internal losses
        Dh_external_loss        %(kJ/kg) - External losses
        RR                      %(kJ/kg) - Rothalpy calculated between impeller's inlet and outlet
        %
    end

    methods
        function cmp = centrifugal_compressor()
            % centrifugal_compressor constructs an instance of this class
            % Nothing is done here: it just creates the compressor object
            %
        end

        function [cmp] = set_geometry(cmp,geometry_specs)
            % set_geometry Sets the cmpressor geometrical specifications
            %
            cmp.N_blades = geometry_specs.N_blades;
            % cmp.blade_thickness = geometry_specs.blade_thickness;
            cmp.blade_thickness_inlet = geometry_specs.blade_thickness_inlet;
            cmp.blade_thickness_outlet = geometry_specs.blade_thickness_outlet;
            %
            cmp.disk_clearance = 1 * 1e-3; % (m) - imposed here. It can be moved outside.
            cmp.radial_clearance = 0.15 * 1e-3; % (m) - imposed here. It can be moved outside.

            %
            cmp.roughness = geometry_specs.roughness;
            cmp.cone_divergence_angle = 5;  % (deg) - imposed here. It can be moved outside, but it is probably always constant
            %
            cmp.r2 = geometry_specs.r2;
            %
            cmp.r1_shroud = cmp.r2 * geometry_specs.ratio_r1shroud_r2;
            cmp.r1_hub = cmp.r1_shroud * geometry_specs.ratio_r1hub_r1shroud;
            cmp.r1_rms = sqrt((cmp.r1_shroud^2 + cmp.r1_hub^2) / 2);
            cmp.A1 = pi*(cmp.r1_shroud^2 - cmp.r1_hub^2);
            cmp.alpha1 = geometry_specs.alpha1;
            cmp.beta1_blade_rms = geometry_specs.beta1_blade_rms;
            %
            cmp.b2 = cmp.r2 * geometry_specs.ratio_b2_r2;
            cmp.A2 = (2 * pi * cmp.r2 - cmp.N_blades * cmp.blade_thickness_outlet) * cmp.b2;
            cmp.beta2_blade = geometry_specs.beta2_blade;
            %
            cmp.r3 = cmp.r2 * geometry_specs.ratio_r3_r2;
            cmp.b3 = cmp.b2;
            cmp.A3 = 2 * pi * cmp.r3 * cmp.b3;
            %
            cmp.r4 = 1.4 * cmp.r3; % (m) - imposed here. It can be moved outside.
            cmp.A4 = pi * (cmp.r4 - cmp.r3)^2;
            %
            cmp.cone_length = cmp.r3 * 1.2; % (m) - imposed here. It can be moved outside.
            cmp.cone_diameter_inlet = 2 * (cmp.r4 - cmp.r3);
            cmp.cone_diameter_outlet = cmp.cone_diameter_inlet + 2 * cmp.cone_length * sind(cmp.cone_divergence_angle / 2);
            cmp.A5 = pi * cmp.cone_diameter_outlet^2 / 4;
            %
            cmp.axial_length = 0.4 * (2 * cmp.r2 - (2 * cmp.r1_shroud + 2 * cmp.r1_hub) / 2); 
            %
        end

        function [cmp] = set_inlet_conditions(cmp,p_01,T_01,fluid)
            % set_inlet_conditions sets the thermodynamic properties at the
            % compressor suction
            %
            cmp.fluid = fluid;
            %
            cmp.T_01 = T_01;
            cmp.p_01 = p_01;
            %
            cmp.cp_fluid = py.CoolProp.CoolProp.PropsSI('CPMASS','T',cmp.T_01,'P',cmp.p_01,cmp.fluid);  % hypothesis that it is constant and equal to the inlet value
            cmp.R_fluid =  py.CoolProp.CoolProp.PropsSI('gas_constant',cmp.fluid) / py.CoolProp.CoolProp.PropsSI('molarmass',cmp.fluid);
            cmp.k_fluid = cmp.cp_fluid / py.CoolProp.CoolProp.PropsSI('CVMASS','T',cmp.T_01,'P',cmp.p_01,cmp.fluid);
            %
            cmp.mu_ref = py.CoolProp.CoolProp.PropsSI('viscosity','T',273.15,'P',101325,cmp.fluid);
            cmp.fluid_as = py.CoolProp.CoolProp.AbstractState("BICUBIC&HEOS", fluid);
            %
            cmp.T_01 = T_01;
            cmp.p_01 = p_01;
            cmp.h_01 = cmp.T_01 * cmp.cp_fluid;
            cmp.s_01 = py.CoolProp.CoolProp.PropsSI('S','T',cmp.T_01,'P',cmp.p_01,cmp.fluid);
            cmp.rho_01 = py.CoolProp.CoolProp.PropsSI('DMASS','T',cmp.T_01,'P',cmp.p_01,cmp.fluid);
            %
        end

        function [cmp] = set_operating_conditions(cmp,mdot,omega_rpm)
            % set_operating_conditions sets the compressor operating
            % conditions
            %
            cmp.mdot = mdot;
            %
            cmp.omega_rpm = omega_rpm;
            cmp.omega_rads = 2 * pi * cmp.omega_rpm / 60;
            %
        end

        function [cmp] = simulation_stage(cmp)
            % It runs the model of the compressor stage (impeller + diffuser) for the
            % specified geometry, inlet conditions, fluid, rotating speed
            % and mass flow rate
            %
            % calculation of the inlet velocities & conditions
            %
            cmp.u1_rms = cmp.omega_rads * cmp.r1_rms;
            cmp.u1_hub = cmp.omega_rads * cmp.r1_hub;
            cmp.u1_shroud = cmp.omega_rads * cmp.r1_shroud;
            %
            f0 = @(x)(cmp.implicit_calculation_1(x));
            x0 = cmp.rho_01;
            cmp.fsolve_options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6);
            % options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-4,'FunctionTolerance',1e-4);
            [rho__1,cmp.convergence.residual_rho1,cmp.convergence.exitflag_rho1] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_rho1 == 1)  %Security check on the convergence of the inlet density calculation
                rho__1 = x0;
            end
            %
            cmp = cmp.calculation_of_section_1(rho__1);
            %
            cmp.c1_meridional = cmp.c1 * cosd(cmp.alpha1);
            cmp.c1_tangential = cmp.c1 * sind(cmp.alpha1);
            %
            cmp.w1_rms = sqrt(cmp.c1^2 + cmp.u1_rms^2);
            cmp.beta1_rms = acosd(cmp.c1 / cmp.w1_rms);
            cmp.w1_shroud = sqrt(cmp.u1_shroud^2 + cmp.c1^2);
            cmp.beta1_shroud = acosd(cmp.c1 / cmp.w1_shroud);
            %
            cmp.w1_hub = sqrt(cmp.u1_hub^2 + cmp.c1^2);
            cmp.beta1_hub = acosd(cmp.c1 / cmp.w1_hub);
            %
            cmp.h_1 = cmp.h_01 - cmp.c1^2 / 2;
            cmp.fluid_as.update(py.CoolProp.CoolProp.PT_INPUTS,cmp.p_1,cmp.T_1)
            cmp.s_1 = cmp.fluid_as.smass;
            %
            % calculation of the throat velocities & conditions
            %
            c_blade_rms = cmp.u1_rms / tand(abs(cmp.beta1_blade_rms));
            c_blade_hub = c_blade_rms;
            c_blade_shroud = c_blade_rms;
            %
            cmp.beta1_blade_hub = -atand(cmp.u1_hub / c_blade_hub);
            cmp.beta1_blade_shroud = -atand(cmp.u1_shroud / c_blade_shroud);
            %
            s_shroud = 2 * pi * cmp.r1_shroud / cmp.N_blades - cmp.blade_thickness_inlet;
            s_hub = 2 * pi * cmp.r1_hub / cmp.N_blades - cmp.blade_thickness_inlet;
            s_rms = 2 * pi * cmp.r1_rms / cmp.N_blades - cmp.blade_thickness_inlet;
            %
            o_shroud = s_shroud * cosd(cmp.beta1_blade_shroud);
            o_hub = s_hub * cosd(cmp.beta1_blade_hub);
            o_rms = s_rms * cosd(cmp.beta1_blade_rms);
            %
            cmp.A1_throat = cmp.N_blades / 2 * (o_hub * (cmp.r1_rms - cmp.r1_hub) +...
                o_rms * (cmp.r1_shroud - cmp.r1_hub) +...
                o_shroud * (cmp.r1_shroud - cmp.r1_rms));
            %
            f0 = @(x)(cmp.implicit_calculation_1throat(x));
            x0 = cmp.rho_1;
            [rho__1throat,cmp.convergence.residual_rho1throat,cmp.convergence.exitflag_rho1throat] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_rho1throat == 1)  %Security check on the convergence of the inlet density calculation
                rho__1throat = x0;
            end
            %
            cmp = calculation_of_section_1throat(cmp,rho__1throat);
            %
            % Impeller outlet
            %
            cmp.u2 = cmp.omega_rads * cmp.r2;
            %
            f0 = @(x)(cmp.implicit_calculation_2(x));
            x0 = cmp.rho_1;
            [rho__2,cmp.convergence.residual_rho2,cmp.convergence.exitflag_rho2] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_rho2 == 1)  %Security check on the convergence of the inlet density calculation
                rho__2 = x0;
            end
            %
            cmp = calculation_of_section_2(cmp,rho__2);
            %
            a2 = sqrt(cmp.k_fluid * cmp.R_fluid * cmp.T_2);
            cmp.Ma_2 = cmp.c2 / a2;
            cmp.p_02 = cmp.p_2 * (1 + ((cmp.k_fluid - 1) / 2) * cmp.Ma_2^2)^(cmp.k_fluid / (cmp.k_fluid - 1));
            % cmp.s_2 = py.CoolProp.CoolProp.PropsSI('S','T',cmp.T_2,'P',cmp.p_2,cmp.fluid);
            cmp.fluid_as.update(py.CoolProp.CoolProp.PT_INPUTS,cmp.p_2,cmp.T_2)
            cmp.s_2 = cmp.fluid_as.smass;
            %
            cmp.RR = (cmp.h_2 + cmp.w2^2 / 2 - cmp.u2^2 / 2) - (cmp.h_1 + cmp.w1_rms^2 / 2 - cmp.u1_rms^2 / 2);   % Rothalpy conservation (RR should be zero)
            %
            % Vaneless Diffuser 
            %
            cmp.h_03 = cmp.h_02;    %total enthalpy conservation
            cmp.T_03 = cmp.T_02;
            %
            f0 = @(x)(cmp.implicit_calculation_3(x));
            x0 = cmp.rho_2;
            [rho__3,cmp.convergence.residual_rho3,cmp.convergence.exitflag_rho3] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_rho3 == 1)  %Security check on the convergence of the inlet density calculation
                rho__3 = x0;
            end
            %
            cmp = calculation_of_section_3(cmp,rho__3);
            %
            cmp.fluid_as.update(py.CoolProp.CoolProp.PT_INPUTS,cmp.p_3,cmp.T_3)
            cmp.s_3 = cmp.fluid_as.smass;
            %
            a3 = sqrt(cmp.k_fluid * cmp.R_fluid * cmp.T_3);
            cmp.Ma_3 = cmp.c3 / a3;
            cmp.p_03 = cmp.p_3* (1 + (cmp.k_fluid - 1) / 2 * cmp.Ma_3^2)^(cmp.k_fluid / (cmp.k_fluid - 1));
            cmp.h_3is = cmp.h_03is - cmp.c3^2 / 2;
            %
            % Stage performance
            %
            cmp.eta_is_tt_stage = (cmp.L_eulero - cmp.Dh_internal_loss - cmp.loss_vaneless_diffuser()) / ...
                (cmp.L_eulero + cmp.Dh_external_loss);
            %
            cmp.PR_tt_stage = cmp.p_03 / cmp.p_01;
            cmp.PR_ts_stage = cmp.p_3 / cmp.p_01;
            cmp.PR_ss_stage = cmp.p_3 / cmp.p_1;
            %
        end
        
        function [cmp] = simulation(cmp)
            % It runs the complete model of the compressor adding to the 
            % impeller and diffuser the volute and the exit cone for the
            % specified geometry, inlet conditions, fluid, rotating speed
            % and mass flow rate
            %
            % Stage (impeller + diffuser)
            %
            [cmp] = simulation_stage(cmp);
            %
            % Volute
            %
            cmp.h_04 = cmp.h_03;    %total enthalpy conservation
            cmp.T_04 = cmp.T_03;
            %
            f0 = @(x)(cmp.implicit_calculation_4(x));
            x0 = cmp.rho_3;
            [rho__4,cmp.convergence.residual_rho4,cmp.convergence.exitflag_rho4] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_rho4 == 1)  %Security check on the convergence of the inlet density calculation
                rho__4 = x0;
            end
            %
            cmp = calculation_of_section_4(cmp,rho__4);
            %
            cmp.h_4 = cmp.h_04 - cmp.c4^2 / 2;
            cmp.fluid_as.update(py.CoolProp.CoolProp.PT_INPUTS,cmp.p_4,cmp.T_4)
            cmp.s_4 = cmp.fluid_as.smass;
            %
            % Exit cone
            %
            cmp.h_05 = cmp.h_04;    %total enthalpy conservation
            cmp.T_05 = cmp.T_04;
            %
            rho__5 = cmp.rho_4; %It is assumed that rho_4 = rho_5 due to the low Mach number in the cone
            %
            cmp = calculation_of_section_5(cmp,rho__5);
            %
            cmp.h_5 = cmp.h_05 - cmp.c5^2 / 2;
            cmp.fluid_as.update(py.CoolProp.CoolProp.PT_INPUTS,cmp.p_5,cmp.T_5)
            cmp.s_5 = cmp.fluid_as.smass;
            %
            % Global performance
            %
            rho_volute_ave = (cmp.rho_3 + cmp.rho_4) / 2;
            Dh_loss_volute = (cmp.p_03 - cmp.p_04) / rho_volute_ave;    % [J/kg] - Total pressure losses converted into hentalpy losses with average density
            %
            rho_cone_ave = (cmp.rho_4 + cmp.rho_5) / 2;
            Dh_loss_cone = (cmp.p_04 - cmp.p_05) / rho_cone_ave;    % [J/kg] - Total pressure losses converted into hentalpy losses with average density
            %
            Dh0_id = cmp.L_eulero - cmp.Dh_internal_loss - cmp.loss_vaneless_diffuser() - Dh_loss_volute - Dh_loss_cone;    %ideal Dh input
            Dh0_act = cmp.L_eulero + cmp.Dh_external_loss;  %actual Dh input
            cmp.eta_is_tt =  Dh0_id / Dh0_act;
            %
            cmp.eta_is_ts = (Dh0_id - cmp.c5^2 / 2) / Dh0_act;    % Total to static is total to total where the exit dynamic term is considered as lost
            %
            cmp.PR_tt = cmp.p_05 / cmp.p_01;
            cmp.PR_ts = cmp.p_5 / cmp.p_01;
            cmp.PR_ss = cmp.p_5 / cmp.p_1;
            %
        end

        function [cmp] = optimisation_geometry_singleObjective(cmp,PR_ts_target)
            % It optimise the compressor geometry to achieve maximum
            % total-to-static isentropic efficiency while also achieving
            % the specified target  total-to-static Pr
            %
            %
            u2_max = 400; % (m/s) - max u2 value to limit Mach
            r2_max = u2_max / cmp.omega_rads; % (m) - consequential max r2 value, given the imposed omega_rpm
            %
            % opt vars      r2              r1_h/r1_s   r1_s/r2     b2/r2   beta2_b     beta1_b_rms     r3/r2                                                   
            problem.lb = [  0.005           0.1         0.4         0.02    -45         -70             1.4];
            problem.ub = [  r2_max          0.7         0.7         0.3     -5          -10             2];
            problem.x0 = [  0.95 * r2_max   0.4         0.6         0.08    -40         -50             1.6];
            problem.solver = 'fmincon';
            problem.options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
                'ScaleProblem',true,'MaxFunctionEvaluations',1e3,'FiniteDifferenceType','forward',...
                'OptimalityTolerance',1e-4,'ConstraintTolerance',1e-3);
            %
            problem.objective = @(x)(cmp.optimisation_geometry_objective(x,PR_ts_target));
            problem.nonlcon = @(x)(cmp.optimisation_geometry_constraints(x,PR_ts_target));
            %
            [cmp.optimisation.x_opt,~,cmp.optimisation.exitflag,cmp.optimisation.output] = fmincon(problem);
            % ms = MultiStart("Display","off");
            % [cmp.optimisation.x_opt,~,cmp.optimisation.exitflag,cmp.optimisation.output] = run(ms,problem,10);
            %
            cmp = cmp.optimisation_geometry_assign(cmp.optimisation.x_opt,PR_ts_target);
            %
            cmp = cmp.simulation();
            %
        end

        function [cmp] = optimisation_geometry_multiObjective(cmp,PR_ts_target,eta_ts_target)
            % It optimise the compressor geometry in a multi-objective sense
            % to simultaneously achieve maximum total-to-static isentropic 
            % efficiency and maximum total-to-static PR
            %
            u2_max = 400; % (m/s) - max u2 value to limit Mach
            r2_max = u2_max / cmp.omega_rads; % (m) - consequential max r2 value, given the imposed omega_rpm
            %
            % opt vars      r2              r1_h/r1_s   r1_s/r2     b2/r2   beta2_b     beta1_b_rms     r3/r2
            problem.lb = [  0.005           0.1         0.4         0.02    -45         -70             1.4];
            problem.ub = [  r2_max          0.7         0.7         0.3     -5          -10             2];
            problem.x0 = [  0.95 * r2_max   0.4         0.6         0.08    -40         -50             1.6];
            problem.solver = 'fmincon';
            problem.options = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
                'ScaleProblem',true,'MaxFunctionEvaluations',1e3,'FiniteDifferenceType','forward',...
                'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-4);
            %
            problem.objective = @(x)(cmp.optimisation_geometry_scalarised_objective(x,PR_ts_target,eta_ts_target));
            problem.nonlcon = [];
            %
            [cmp.optimisation.x_opt,cmp.optimisation.fval,cmp.optimisation.exitflag,cmp.optimisation.output] = fmincon(problem);
            % ms = MultiStart("Display","off");
            % [cmp.optimisation.x_opt,~,cmp.optimisation.exitflag,cmp.optimisation.output] = run(ms,problem,10);
            %
            cmp = cmp.optimisation_geometry_assign(cmp.optimisation.x_opt,PR_ts_target);
            %
            cmp = cmp.simulation();
            %
        end

        function [cmp] = optimisation_geometry_multiObjective_multistart(cmp,PR_ts_target,eta_ts_target)
            %
            u2_max = 400; % (m/s) - max u2 value to limit Mach
            r2_max = u2_max / cmp.omega_rads; % (m) - consequential max r2 value, given the imposed omega_rpm
            %
            problem.lb = [0.005         0.1     0.4     0.02    -45     -70     1.4];
            problem.ub = [r2_max        0.7     0.7     0.3     -5      -10     2];
            problem.x0 = [0.95 * r2_max 0.4     0.6     0.08    -40     -50     1.6];
            problem.solver = 'fmincon';
            problem.options = optimoptions('fmincon','Display','off','Algorithm','sqp',...
                'ScaleProblem',true,'MaxFunctionEvaluations',1e3,'FiniteDifferenceType','forward',...
                'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-4);
            %
            problem.objective = @(x)(cmp.optimisation_geometry_scalarised_objective(x,PR_ts_target,eta_ts_target));
            problem.nonlcon = [];
            %
            % [cmp.optimisation.x_opt,cmp.optimisation.fval,cmp.optimisation.exitflag,cmp.optimisation.output] = fmincon(problem);
            ms = MultiStart("Display","off");
            [cmp.optimisation.x_opt,~,cmp.optimisation.exitflag,cmp.optimisation.output] = run(ms,problem,10);
            %
            cmp = cmp.optimisation_geometry_assign(cmp.optimisation.x_opt,PR_ts_target);
            %
            cmp = cmp.simulation();
            %
        end

        %
        % Validation
        %
        
        function [cmp] = set_geometry_eckardt_impeller_O(cmp)
            % It sets the compressor geometrical specifications
            % according to the data related the Impeller "O" of Eckardt
            % experiments for validation purposes - (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            %
            cmp.N_blades = 20;
            cmp.blade_thickness_inlet = 2.11 * 1e-3;
            cmp.blade_thickness_outlet = 1.08 * 1e-3;
            %
            cmp.disk_clearance = 0.372 * 1e-3;   %(m)
            cmp.radial_clearance = 0.372 * 1e-3;    %(m)
            %
            cmp.roughness = 0.002 * 1e-3;
            cmp.cone_divergence_angle = 5;  % (deg) - imposed here. It can be moved outside, but it is probably always constant
            %
            cmp.r2 = 0.2;
            %
            cmp.r1_shroud = 0.14;
            cmp.r1_hub = 0.045;
            cmp.r1_rms = sqrt((cmp.r1_shroud^2 + cmp.r1_hub^2) / 2);
            cmp.A1 = pi*(cmp.r1_shroud^2 - cmp.r1_hub^2);
            cmp.alpha1 = 0;
            %
            % Eckardt specifies beta1_blade_hub and beta1_blade_shroud.
            % Here we use one of the two to calculate beta1_blade_rms (need
            % by the model) - the other could be used as well.
            beta1_blade__hub = -32;
            beta1_blade__shroud = -63;
            cmp.beta1_blade_rms = mean(...
                [- atand(cmp.r1_rms / cmp.r1_hub * tand(abs(beta1_blade__hub))), ...
                - atand(cmp.r1_rms / cmp.r1_shroud * tand(abs(beta1_blade__shroud)))]);
            %
            cmp.b2 = 0.026;
            cmp.A2 = (2 * pi * cmp.r2 - cmp.N_blades * cmp.blade_thickness_outlet) * cmp.b2;
            cmp.beta2_blade = 0;
            %
            cmp.r3 = cmp.r2 * 1.69;
            cmp.b3 = cmp.b2 * 0.51;
            cmp.A3 = 2 * pi * cmp.r3 * cmp.b3;
            %
            cmp.r4 = 1.3 * cmp.r3; % (m)
            cmp.A4 = pi * (cmp.r4 - cmp.r3)^2;
            %
            cmp.cone_length = cmp.r3 * 1.1;
            cmp.cone_diameter_inlet = 2 * (cmp.r4 - cmp.r3);
            cmp.cone_diameter_outlet = cmp.cone_diameter_inlet + 2 * cmp.cone_length * sind(cmp.cone_divergence_angle / 2);
            cmp.A5 = pi * cmp.cone_diameter_outlet^2 / 4;
            %
            cmp.axial_length = 0.13;
            %
        end

        function [cmp] = set_geometry_eckardt_impeller_A(cmp)
            %
            % It sets the compressor geometrical specifications
            % according to the data related the Impeller "A" of Eckardt
            % experiments for validation purposes - (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            %
            cmp.N_blades = 20;
            cmp.blade_thickness_inlet = 2.11 * 1e-3;
            cmp.blade_thickness_outlet = 1.08 * 1e-3;
            %
            cmp.disk_clearance = 0.235 * 1e-3;   %(m)
            cmp.radial_clearance = 0.19 * 1e-3;    %(m)
            %
            cmp.roughness = 0.002 * 1e-3;
            cmp.cone_divergence_angle = 5;  % (deg) - imposed here. It can be moved outside, but it is probably always constant
            %
            cmp.r2 = 0.2;
            %
            cmp.r1_shroud = 0.14;
            cmp.r1_hub = 0.06;
            cmp.r1_rms = sqrt((cmp.r1_shroud^2 + cmp.r1_hub^2) / 2);
            cmp.A1 = pi*(cmp.r1_shroud^2 - cmp.r1_hub^2);
            cmp.alpha1 = 0;
            %
            % Eckardt specifies beta1_blade_hub and beta1_blade_shroud.
            % Here we use one of the two to calculate beta1_blade_rms (need
            % by the model) - the other could be used as well.
            beta1_blade__hub = -40;
            beta1_blade__shroud = -63;
            cmp.beta1_blade_rms = mean(...
                [- atand(cmp.r1_rms / cmp.r1_hub * tand(abs(beta1_blade__hub))), ...
                - atand(cmp.r1_rms / cmp.r1_shroud * tand(abs(beta1_blade__shroud)))]);
            %
            cmp.b2 = 0.026;
            cmp.A2 = (2 * pi * cmp.r2 - cmp.N_blades * cmp.blade_thickness_outlet) * cmp.b2;
            cmp.beta2_blade = -30;
            %
            cmp.r3 = cmp.r2 * 1.69;
            cmp.b3 = cmp.b2 * 0.51;
            cmp.A3 = 2 * pi * cmp.r3 * cmp.b3;
            %
            cmp.r4 = 1.3 * cmp.r3; % (m) - It comes from the mass conservation by imposing rho_4 = rho_3 and c4 = c3_meridional
            cmp.A4 = pi * (cmp.r4 - cmp.r3)^2;
            %
            cmp.cone_length = cmp.r3 * 1.1;
            cmp.cone_diameter_inlet = 2 * (cmp.r4 - cmp.r3);
            cmp.cone_diameter_outlet = cmp.cone_diameter_inlet + 2 * cmp.cone_length * sind(cmp.cone_divergence_angle / 2);
            cmp.A5 = pi * cmp.cone_diameter_outlet^2 / 4;
            %
            cmp.axial_length = 0.13;
            %
        end
        
        function [cmp] = set_geometry_eckardt_impeller_B(cmp)
            %
            % It sets the compressor geometrical specifications
            % according to the data related the Impeller "A" of Eckardt
            % experiments for validation purposes - (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            %
            cmp.N_blades = 20;
            cmp.blade_thickness_inlet = 2.11 * 1e-3;
            cmp.blade_thickness_outlet = 1.08 * 1e-3;
            %
            cmp.disk_clearance = 0.372 * 1e-3;   %(m)
            cmp.radial_clearance = 0.372 * 1e-3;    %(m)
            %
            cmp.roughness = 0.002 * 1e-3;
            cmp.cone_divergence_angle = 5;  % (deg) - imposed here. It can be moved outside, but it is probably always constant
            %
            cmp.r2 = 0.2;
            %
            cmp.r1_shroud = 0.14;
            cmp.r1_hub = 0.0959;
            cmp.r1_rms = sqrt((cmp.r1_shroud^2 + cmp.r1_hub^2) / 2);
            cmp.A1 = pi*(cmp.r1_shroud^2 - cmp.r1_hub^2);
            cmp.alpha1 = 0;
            %
            % Eckardt specifies beta1_blade_hub and beta1_blade_shroud.
            % Here we use one of the two to calculate beta1_blade_rms (need
            % by the model) - the other could be used as well.
            beta1_blade__hub = -45;
            beta1_blade__shroud = -60;
            cmp.beta1_blade_rms = mean(...
                [- atand(cmp.r1_rms / cmp.r1_hub * tand(abs(beta1_blade__hub))), ...
                - atand(cmp.r1_rms / cmp.r1_shroud * tand(abs(beta1_blade__shroud)))]);
            %
            cmp.b2 = 0.026;
            cmp.A2 = (2 * pi * cmp.r2 - cmp.N_blades * cmp.blade_thickness_outlet) * cmp.b2;
            cmp.beta2_blade = -40;
            %
            cmp.r3 = cmp.r2 * 1.69;
            cmp.b3 = cmp.b2 * 0.51;
            cmp.A3 = 2 * pi * cmp.r3 * cmp.b3;
            %
            cmp.r4 = 1.3 * cmp.r3; % (m) - It comes from the mass conservation by imposing rho_4 = rho_3 and c4 = c3_meridional
            cmp.A4 = pi * (cmp.r4 - cmp.r3)^2;
            %
            cmp.cone_length = cmp.r3 * 1.1;
            cmp.cone_diameter_inlet = 2 * (cmp.r4 - cmp.r3);
            cmp.cone_diameter_outlet = cmp.cone_diameter_inlet + 2 * cmp.cone_length * sind(cmp.cone_divergence_angle / 2);
            cmp.A5 = pi * cmp.cone_diameter_outlet^2 / 4;
            %
            cmp.axial_length = 0.13;
            %
        end

        %
        % Loss correlations
        %

        function [loss_tab] = loss_analysis(cmp)
            % It generates a table summarising the contribution of each 
            % loss mechanism
            %
            rho_volute_ave = (cmp.rho_3 + cmp.rho_4) / 2;
            Dh_loss_volute = (cmp.p_03 - cmp.p_04) / rho_volute_ave;    % [J/kg] - Total pressure losses converted into hentalpy losses with average density
            %
            rho_cone_ave = (cmp.rho_4 + cmp.rho_5) / 2;
            Dh_loss_cone = (cmp.p_04 - cmp.p_05) / rho_cone_ave;    % [J/kg] - Total pressure losses converted into hentalpy losses with average density
            %
            loss_tab = table(...
                cmp.loss_inducer_incidence(), cmp.loss_impeller_blade_loading(), cmp.loss_impeller_clearance(), ...
                cmp.loss_impeller_disk_friction(), cmp.loss_impeller_mixing(), cmp.loss_impeller_recirculation(), ...
                cmp.loss_impeller_skin_friction(), cmp.loss_vaneless_diffuser(),...
                Dh_loss_volute, Dh_loss_cone,...
                'VariableNames',[...
                "Inducer_Incidence","Impeller_Blade_Loading","Impeller_Clearance",...
                "Impeller_Disk_Friction","Impeller_Mixing","Impeller_Recirculation",...
                "Impeller_Skin_Friction","Vaneless_diffuser","Volute","Exhaust Cone"],...
                'RowNames',"Losses (J/kg)");
            %
        end

        %
        % Figures
        %

        function [varargout] = plot_geometry(cmp)
            %
            x1 = linspace(0,cmp.axial_length - cmp.b2,200);
            y1 = nan(size(x1));
            % options = optimoptions('fsolve','Display','off','OptimalityTolerance',1e-4,'FunctionTolerance',1e-4);
            for i = 1 : length(x1)
                fun = @(y)(cmp.ellipse_internal(x1(i),y));
                y1(i) = fsolve(fun,0.04,cmp.fsolve_options);
            end
            %
            x2 = linspace(0,cmp.axial_length,200);
            y2 = nan(size(x2));
            for i = 1 : length(x2)
                fun = @(y)cmp.ellipse_external(x2(i),y);
                y2(i) = fsolve(fun,0.04,cmp.fsolve_options);
            end
            %
            impeller_profile(:,1) = [...
                0 0,...
                x1,...
                cmp.axial_length - cmp.b2 cmp.axial_length,...
                fliplr(x2)]';
            %
            impeller_profile(:,2) = [...
                cmp.r1_hub cmp.r1_shroud,...
                y1,...
                cmp.r2 cmp.r2,...
                fliplr(y2)]';
            %
            vaneless_profile(:,1) = [...
                cmp.axial_length - cmp.b2 cmp.axial_length - cmp.b2,...
                cmp.axial_length - cmp.b2 cmp.axial_length,...
                cmp.axial_length cmp.axial_length,...
                cmp.axial_length cmp.axial_length - cmp.b2]';

            vaneless_profile(:,2) = [...
                cmp.r2 cmp.r3,...
                cmp.r3 cmp.r3,...
                cmp.r3 cmp.r2,...
                cmp.r2 cmp.r2]';
            %
            r_v = cmp.cone_diameter_inlet / 2;
            x_center_v = cmp.axial_length - cmp.b3 / 2;
            y_center_v = cmp.r3 + r_v;
            %
            a = -2 * x_center_v;    %parameters of the general equation of circle - x^2 + y^2 + ax + by + c = 0
            b = -2 * y_center_v;
            c = x_center_v^2 + y_center_v^2 - r_v^2;
            %
            x_volute = linspace(x_center_v - r_v,x_center_v + r_v);
            alpha = 1;  %parameters of the equation of second degree to find y - alpha*y^2 + beta*y + gamma = 0
            beta = b;
            gamma = x_volute.^2 + a * x_volute + c;
            %
            y_plus_volute = (-beta + sqrt(beta^2 - 4 * alpha * gamma)) / (2 * alpha);
            y_minus_volute = (-beta - sqrt(beta^2 - 4 * alpha * gamma)) / (2 * alpha);
            %
            x_volute = [x_volute fliplr(x_volute)];
            y_volute = [y_plus_volute fliplr(y_minus_volute)]';
            %
            volute_profile(:,1) = x_volute;
            volute_profile(:,2) = y_volute;
            %
            if nargout == 0 %just the plot
                %
                gcf;
                %
                k1 = plot(impeller_profile(:,1),impeller_profile(:,2),...
                    "color",[0.7961 0.0941 0.1137],"LineWidth",2);
                hold on
                k2 = plot(vaneless_profile(:,1),vaneless_profile(:,2),...
                    "color",[0.1294 0.4431 0.7098],"LineWidth",2);
                k3 = plot(volute_profile(:,1),volute_profile(:,2),...
                    "color",[0.1294 0.4431 0.7098],"LineWidth",2);
                axis equal
                xlabel('Radial dimension (m)')
                ylabel('Axial dimension (m)')
                grid on
                box on
                set(gca,"FontSize",24)
                legend([k1 k2 k3],{'Impeller','Vaneless Diffuser','Volute outlet'},'Location','southeast')
                %
            elseif nargout == 2
                varargout{1} = impeller_profile;
                varargout{2} = vaneless_profile;
            elseif nargout == 3
                varargout{1} = impeller_profile;
                varargout{2} = vaneless_profile;
                varargout{3} = volute_profile;
            else
                error("Only 0 or 2 output arguments are supported")
            end


        end
    end

    methods (Access = private, Hidden = true)

        function [cmp] = calculation_of_section_1(cmp,rho_1)
            %
            cmp.rho_1 = rho_1;
            %
            cmp.c1 = cmp.mdot / (cmp.A1 * cmp.rho_1);
            cmp.T_1 = cmp.T_01 - cmp.c1^2 / (2 * cmp.cp_fluid);
            a1 = sqrt(cmp.k_fluid * cmp.R_fluid * cmp.T_1);
            cmp.Ma_1 = cmp.c1 / a1;
            cmp.p_1 = cmp.p_01 /...
                ((1 + ((cmp.k_fluid - 1) / 2) * cmp.Ma_1^2)^(cmp.k_fluid / (cmp.k_fluid - 1))); % [Pa]
            %
        end

        function [err] = implicit_calculation_1(cmp,rho_1guess)
            %
            cmp = cmp.calculation_of_section_1(rho_1guess);
            %
            rho__1 = cmp.p_1 / (cmp.R_fluid * cmp.T_1);
            err = rho_1guess - rho__1;
            %
        end

        function [cmp] = calculation_of_section_1throat(cmp,rho_1throat)
            %
            cmp.rho_1throat = rho_1throat;
            %
            cmp.w1_throat = cmp.mdot / (cmp.rho_1throat * cmp.A1_throat);
            w1_throat_tangential = cmp.w1_throat * sind(abs(cmp.beta1_blade_rms));
            w1_throat_meridional = cmp.w1_throat * cosd(cmp.beta1_blade_rms);
            c1_throat_tangential = cmp.u1_rms - w1_throat_tangential;
            c1_throat = sqrt(c1_throat_tangential^2 + w1_throat_meridional^2);
            %
            h1_throat = cmp.h_01 - c1_throat^2 / 2;  % we neglect u, as the variation of radius is limited
            cmp.T_1throat = h1_throat / cmp.cp_fluid;
            a_1throat = sqrt(cmp.k_fluid * cmp.R_fluid * cmp.T_1throat);
            cmp.Ma_1throat = c1_throat / a_1throat;
            cmp.p_1throat = cmp.p_01 / ((1 + (cmp.k_fluid - 1) / 2 * cmp.Ma_1throat^2)^(cmp.k_fluid / (cmp.k_fluid - 1))); %[Pa]
            %
        end

        function [err] = implicit_calculation_1throat(cmp,rho_1throatguess)
            %
            cmp = cmp.calculation_of_section_1throat(rho_1throatguess);
            %
            rho__1throat = cmp.p_1throat / (cmp.R_fluid * cmp.T_1throat);
            err = rho_1throatguess - rho__1throat;
            %
        end

        function [cmp] = calculation_of_section_2(cmp,rho_2)
            %
            cmp.rho_2 = rho_2;
            %
            cmp.c2_meridional = cmp.mdot / (cmp.rho_2 * cmp.A2);
            cmp.c2_tangential = cmp.slip_factor() * cmp.u2 - cmp.c2_meridional * tand(abs(cmp.beta2_blade));
            cmp.c2 = sqrt(cmp.c2_meridional^2 + cmp.c2_tangential^2);
            %
            cmp.w2_tangential = cmp.u2 - cmp.c2_tangential;
            cmp.w2 = sqrt(cmp.w2_tangential^2 + cmp.c2_meridional^2);
            cmp.alpha2 = acosd(cmp.c2_meridional / cmp.c2);
            cmp.beta2 = acosd(cmp.c2_meridional / cmp.w2);
            %
            cmp.L_eulero = cmp.u2 * cmp.c2_tangential - cmp.u1_rms * cmp.c1_tangential;
            %
            cmp.h_02 = cmp.h_01 + cmp.L_eulero;
            cmp.T_02 = cmp.h_02 / cmp.cp_fluid;
            %
            cmp.h_2 = cmp.h_02 - cmp.c2^2 / 2;
            cmp.T_2 = cmp.h_2 / cmp.cp_fluid;
            %
            cmp.Dh_internal_loss = cmp.loss_inducer_incidence() + cmp.loss_impeller_skin_friction() + cmp.loss_impeller_blade_loading() + cmp.loss_impeller_clearance();
            cmp.Dh_external_loss = cmp.loss_impeller_mixing() + cmp.loss_impeller_disk_friction() + cmp.loss_impeller_recirculation();
            %
            cmp.h_02is = cmp.h_02 - cmp.Dh_internal_loss;
            cmp.h_2is = cmp.h_02is - cmp.c2^2 / 2;
            %
            cmp.T_02is = cmp.h_02is / cmp.cp_fluid;
            cmp.T_2is = cmp.h_2is / cmp.cp_fluid;
            %
            cmp.p_02is = cmp.p_01 * (cmp.T_02is / cmp.T_01)^(cmp.k_fluid / (cmp.k_fluid - 1));
            cmp.p_2 = cmp.p_02is * (cmp.T_02is / cmp.T_2is)^(cmp.k_fluid / (1 - cmp.k_fluid));
            %
        end

        function [err] = implicit_calculation_2(cmp,rho_2guess)
            %
            cmp = cmp.calculation_of_section_2(rho_2guess);
            %
            rho__2 = cmp.p_2 / (cmp.R_fluid * cmp.T_2);
            err = rho_2guess - rho__2;
            %
        end

        function [cmp] = calculation_of_section_3(cmp,rho_3)
            %
            cmp.rho_3 = rho_3;
            %
            cmp.c3_meridional = cmp.mdot / (cmp.rho_3 * cmp.A3);
            %
            % Implicit calculation of the Reynolds (to determine T_3, required to calculate rho_3 and check the guess)
            %
            f0 = @(x)(cmp.implicit_calculation_3_Re(x));
            D_hydr = 2 * cmp.b3;
            %
            cmp.fluid_as.update(py.CoolProp.CoolProp.DmassP_INPUTS,cmp.rho_3,cmp.p_2)
            mu_3guess = cmp.fluid_as.viscosity;
            Re_3guess = cmp.rho_3 * cmp.c2 * D_hydr / mu_3guess;
            x0 = Re_3guess;
            [Re__3,cmp.convergence.residual_Re3,cmp.convergence.exitflag_Re3] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_Re3 == 1)  %Security check on the convergence of the Reynolds number calculation on Section 3
                Re__3 = x0;
            end
            %
            cmp = cmp.calculation_of_section_3_Re(Re__3);
            %
            % Implicit calculation of the isentropic vaneless diffuser (to determine p_3, required to calculate rho_3 and check the guess)
            %
            cmp.h_03is = cmp.h_03 - cmp.loss_vaneless_diffuser();
            cmp.T_03is = cmp.h_03is / cmp.cp_fluid;
            cmp.p_03is = cmp.p_02 * (cmp.T_02 / cmp.T_03is)^(cmp.k_fluid  / (1 - cmp.k_fluid));
            %
            f0 = @(x)(cmp.implicit_calculation_3_rhois(x));
            x0 = cmp.rho_3;
            [rho__3is,cmp.convergence.residual_rho3is,cmp.convergence.exitflag_rho3is] = fsolve(f0,x0,cmp.fsolve_options);
            %
            if not(cmp.convergence.exitflag_rho3is == 1)  %Security check on the convergence of the Reynolds number calculation on Section 3
                rho__3is = x0;
            end
            %
            cmp = cmp.calculation_of_section_3_rhois(rho__3is);
            %
        end

        function [err] = implicit_calculation_3(cmp,rho_3guess)
            %
            cmp = cmp.calculation_of_section_3(rho_3guess);
            %
            rho__3 = cmp.p_3 / (cmp.R_fluid * cmp.T_3);
            %
            err = rho_3guess - rho__3;
            %
        end

        function [cmp] = calculation_of_section_3_Re(cmp,Re_3)
            %
            cmp.Re_3 = Re_3;
            %
            cf = 0.0058 * (1.8 * 1e5 / cmp.Re_3)^0.2; %(Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            %
            cmp.c3_tangential = cmp.c2_tangential / (cmp.r3 / cmp.r2 + ....
                (2 * pi * cf * cmp.rho_2 * cmp.c2_tangential * (cmp.r3^2 - cmp.r2 * cmp.r3)) / cmp.mdot);
            cmp.c3 = sqrt(cmp.c3_meridional^2 + cmp.c3_tangential^2);
            %
            cmp.h_3 = cmp.h_03 - cmp.c3^2 / 2;
            cmp.T_3 = cmp.T_03 - (cmp.c3^2 / (2 * cmp.cp_fluid));
            %
        end

        function [err] = implicit_calculation_3_Re(cmp,Re_3guess)
            %
            cmp = cmp.calculation_of_section_3_Re(Re_3guess);
            %
            cmp.fluid_as.update(py.CoolProp.CoolProp.DmassT_INPUTS,cmp.rho_3,cmp.T_3)
            mu_3 = cmp.fluid_as.viscosity;
            %
            D_hydr = 2 * cmp.b3;
            Re__3 = (cmp.rho_3 * cmp.c3 * D_hydr) / mu_3;
            %
            err = Re_3guess - Re__3;
            %
        end

        function [cmp] = calculation_of_section_3_rhois(cmp,rho_3is)
            %
            cmp.rho_3is = rho_3is;
            %
            c3_meridional_is = cmp.mdot / (cmp.rho_3is * cmp.A3);
            c3_tangential_is = cmp.c2_tangential * (cmp.r2 / cmp.r3);   %free vortex without losses
            c3_is = sqrt(c3_meridional_is^2 + c3_tangential_is^2);
            cmp.T_3is = cmp.T_03is - c3_is^2 / (2 * cmp.cp_fluid);
            cmp.p_3 = cmp.p_03is * (cmp.T_03is / cmp.T_3is)^(cmp.k_fluid / (1 - cmp.k_fluid));
            %
        end

        function [err] = implicit_calculation_3_rhois(cmp,rho_3isguess)
            %
            cmp = cmp.calculation_of_section_3_rhois(rho_3isguess);
            %
            rho__3is = cmp.p_3 / (cmp.R_fluid * cmp.T_3is);
            %
            err = rho_3isguess - rho__3is;
            %
        end

        function [cmp] = calculation_of_section_4(cmp,rho_4)
            %
            cmp.rho_4 = rho_4;
            %
            cmp.c4 = cmp.mdot / (cmp.rho_4 * cmp.A4);
            cmp.T_4 = cmp.T_04 - cmp.c4^2 / (2 * cmp.cp_fluid);
            %
            cmp.p_04 = cmp.p_03 - (cmp.p_03 - cmp.p_3) * cmp.loss_factor_volute();
            %
            a4 = sqrt(cmp.k_fluid * cmp.R_fluid * cmp.T_4);
            cmp.Ma_4 = cmp.c4 / a4;
            cmp.p_4 = cmp.p_04 /...
                ((1 + ((cmp.k_fluid - 1) / 2) * cmp.Ma_4^2)^(cmp.k_fluid / (cmp.k_fluid - 1))); % [Pa]
            %
        end

        function [err] = implicit_calculation_4(cmp,rho_4guess)
            %
            cmp = cmp.calculation_of_section_4(rho_4guess);
            %
            rho__4 = cmp.p_4 / (cmp.R_fluid * cmp.T_4);
            err = rho_4guess - rho__4;
            %
        end

        function [cmp] = calculation_of_section_5(cmp,rho_5)
            %
            cmp.rho_5 = rho_5;
            %
            cmp.c5 = cmp.mdot / (cmp.rho_5 * cmp.A5);
            cmp.T_5 = cmp.T_05 - cmp.c5^2 / (2 * cmp.cp_fluid);
            %
            cmp.p_05 = cmp.p_04 - (cmp.p_03 - cmp.p_3) * cmp.loss_factor_cone();
            %
            a5 = sqrt(cmp.k_fluid * cmp.R_fluid * cmp.T_5);
            cmp.Ma_5 = cmp.c5 / a5;
            cmp.p_5 = cmp.p_05 /...
                ((1 + ((cmp.k_fluid - 1) / 2) * cmp.Ma_5^2)^(cmp.k_fluid / (cmp.k_fluid - 1))); % [Pa]
            %
        end

        %
        % Loss correlations
        %

        function [slip_factor] = slip_factor(cmp)
            %
            % Reference 
            % (Wiesner, 1967) - doi:10.1115/1.3616734 
            % (Meroni, 2018) - doi:10.1016/j.apenergy.2018.09.210
            % (Aungier, 2000) - doi:10.1115/1.800938
            %
            sigma = 1 - sqrt(cosd(cmp.beta2_blade)) / cmp.N_blades^0.7;
            %
            sigma_star = sind(19 + 0.2 * (90 - abs(cmp.beta2_blade)));
            %
            A = (sigma - sigma_star) / (1 - sigma_star);    % (r1_rms/r2)_lim
            B = (cmp.r1_rms / cmp.r2 - A) / (1 - A);
            %
            if cmp.r1_rms / cmp.r2 < A
                slip_factor = sigma;
            else
                % sigma_cor = sigma  *( 1 - B^sqrt((90-abd(beta2_blade)/10)));
                slip_factor = sigma * ( 1 - B^sqrt((90 - abs(cmp.beta2_blade) / 10)));
            end
            %
        end

        function [Dh_loss] = loss_inducer_incidence(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            beta1_optimal = atand(cmp.A1 / cmp.A1_throat * tand(abs(cmp.beta1_blade_rms)));
            E = sind(abs(beta1_optimal - cmp.beta1_rms));
            Dh_loss = cmp.w1_rms^2 * E^2 / 2;
            %
        end

        function [Dh_loss] = loss_impeller_skin_friction(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            cmp.fluid_as.update(py.CoolProp.CoolProp.DmassT_INPUTS,cmp.rho_2,cmp.T_2)
            mu_2 = cmp.fluid_as.viscosity;
            %
            E = (cosd(cmp.beta1_blade_shroud) + cosd(cmp.beta1_blade_hub)) / 2;
            Lb = (pi / 8) * (2 * cmp.r2 - (cmp.r1_shroud + cmp.r1_hub) - cmp.b2 + 2 * cmp.axial_length) * (2 / (E + cosd(cmp.beta2_blade)));
            F = 2 * cmp.r2 / ((cmp.N_blades / pi * cosd(cmp.beta2_blade)) + 2 * cmp.r2 / cmp.b2);
            G = tand(cmp.beta1_blade_shroud);
            lamda = cmp.r1_hub / cmp.r1_shroud;
            H = 2 * cmp.r1_shroud / (2 / (1 - lamda) + (2 * cmp.N_blades / (pi * (1 + lamda))) * sqrt(1 + G^2 * (1 + lamda^2 / 2)));
            D_hydr = F + H; % (m) hydraulic diameter
            %
            Re = (cmp.rho_2 * cmp.c2 * D_hydr) / mu_2;
            Re_e = (Re - 2000) * (cmp.roughness / D_hydr);
            %
            Cf = cmp.friction_factor(Re,Re_e,D_hydr);
            %
            K_1 = sqrt((cmp.w1_rms^2 + cmp.w2^2) / 2);
            K_2 = sqrt((cmp.w1_throat^2 + cmp.w2^2) / 2);
            %
            W = max(K_1,K_2);
            %
            Dh_loss = 4 * Cf * (Lb / D_hydr) * W^2;
            %
        end

        function [Dh_loss] = loss_impeller_blade_loading(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            E = (cmp.N_blades / pi * (1 - cmp.r1_shroud / cmp.r2) + 2 * (cmp.r1_shroud / cmp.r2));
            F = 0.75 * cmp.u2 * cmp.c2_tangential / cmp.u2^2 * cmp.w2 / cmp.w1_shroud;
            Df = 1 - (cmp.w2 / cmp.w1_shroud) + F / E;
            %
            Dh_loss = 0.05 * Df^2 * cmp.u2^2;
            %
        end

        function [Dh_loss] = loss_impeller_clearance(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            E = 4 * pi / (cmp.b2 * cmp.N_blades);
            F = (cmp.r1_shroud^2 - cmp.r1_hub^2) / ((cmp.r2 - cmp.r1_shroud) * (1 + cmp.rho_2 / cmp.rho_1));
            Dh_loss = 0.6 * (cmp.radial_clearance / cmp.b2) * abs(cmp.c2_tangential) * ...
                sqrt(E * F * abs(cmp.c2_tangential) * cmp.c1);
            %
        end

        function [Dh_loss] = loss_impeller_mixing(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            b_star = 1;
            csi = 0.15;
            epsilon_w = (-0.07 + sqrt(0.07^2 + 4 * 0.93 * csi)) / (2 * 0.93);
            E = tand(cmp.alpha2);
            Dh_loss = (1 / (1 + E^2))*((1 - epsilon_w - b_star) / (1 - epsilon_w))^2 * cmp.c2^2 / 2;
            %
        end

        function [Dh_loss] = loss_impeller_disk_friction(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            cmp.fluid_as.update(py.CoolProp.CoolProp.DmassT_INPUTS,cmp.rho_2,cmp.T_2)
            mu_2 = cmp.fluid_as.viscosity;
            %
            rho_ave = (cmp.rho_1 + cmp.rho_2) / 2;
            Re = cmp.rho_2 * cmp.u2 * cmp.r2 / mu_2;
            %
            if Re  > 3 * 1e5
                Kf = 0.102 * (cmp.disk_clearance / cmp.b2)^0.1 / Re^0.2;
            else
                Kf = 3.7 * (cmp.disk_clearance / cmp.b2)^0.1 / Re^0.5;
            end
            %
            Dh_loss = 0.25 * rho_ave * cmp.u2^3 * cmp.r2^2 * Kf / cmp.mdot;
            %
        end

        function [Dh_loss] = loss_impeller_recirculation(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % (Oh et al., 1997) - 10.1243/0957650971537231
            %
            E = (cmp.N_blades / pi * (1 - cmp.r1_shroud / cmp.r2) + 2 * (cmp.r1_shroud / cmp.r2));
            F = 0.75 * cmp.u2 * cmp.c2_tangential / cmp.u2^2 * cmp.w2 / cmp.w1_shroud;
            Df = 1 - (cmp.w2 / cmp.w1_shroud) + F / E;
            %
            Dh_loss = 8 * 1e-5 * sinh(3.5 * deg2rad(cmp.alpha2)^3)* Df^2 * cmp.u2^2;
            %
        end

        function [Dh_loss] = loss_vaneless_diffuser(cmp)
            %
            % (Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            %
            cf = 0.005 * (1.8 * 1e5 / cmp.Re_3)^0.2;
            Dh_loss = (cf * cmp.r2 * (1 - (cmp.r2 / cmp.r3)^1.5) * cmp.c2^2)/(1.5 * cmp.b2 * cosd(cmp.alpha2));
            %
        end

        function [loss_factor] = loss_factor_volute(cmp)
            %
            % (Aungier, 2000) - doi:10.1115/1.800938
            %
            % Meridional velocity loss
            %
            loss_1 = (cmp.c3_meridional / cmp.c3)^2;
            %
            % Tangential velocity loss
            %
            sp = cmp.r3 * cmp.c3_tangential / (cmp.r4 * cmp.c4);
            %
            if sp >= 1
                loss_2 = 0.5 * cmp.r3 / cmp.r4 * (cmp.c3_tangential / cmp.c3)^2 * (1 - 1 / sp^2);
            else
                loss_2 = cmp.r3 / cmp.r4 * (cmp.c3_tangential / cmp.c3)^2 * (1 - 1 / sp)^2;
            end
            %
            % Skin friction losses
            %
            L_ave = pi * (cmp.r3 + cmp.r4) / 2; % (m) - Average length travelled by the fluid in the volute
            D_hydr = sqrt(4 * cmp.A4 / pi);
            cmp.fluid_as.update(py.CoolProp.CoolProp.DmassT_INPUTS,cmp.rho_4,cmp.T_4)
            mu_4 = cmp.fluid_as.viscosity;
            Re = (cmp.rho_4 * cmp.c4 * D_hydr) / mu_4;
            Re_e = (Re - 2000) * (cmp.roughness / D_hydr);
            loss_3 = 4 * cmp.friction_factor(Re,Re_e,D_hydr) * (cmp.c4 / cmp.c3)^2 * L_ave / D_hydr;
            %
            loss_factor = loss_1 + loss_2 + loss_3; %Total pressure loss factor
            %
        end

        function [loss_factor] = loss_factor_cone(cmp)
            %
            % (Aungier, 2000) - doi:10.1115/1.800938
            %
            loss_factor = ((cmp.c4 - cmp.c5) / cmp.c3)^2; %Total pressure loss factor
            %
        end

        function [Cf] = friction_factor(cmp,Re,Re_e,D_hydr)
            %
            % Calculation of the friction factor - (Aungier, 2000) -
            % 10.1115/1.800938 - pp. 74 - 75
            %
            f0 = @(x)(2.51 * 10^(1 / (2 * x)) - Re * x);    %colebrook function - CFS
            x0 = 1;
            cfs = fsolve(f0,x0,cmp.fsolve_options)^2 / 4;
            %
            f0 = @(x)(10^(1 / (2 * x)) * cmp.roughness - 3.71 * D_hydr);    %colebrook function - CFR
            x0 = 0.5;
            cfr = fsolve(f0,x0,cmp.fsolve_options)^2 / 4;
            %
            if Re_e < 60
                Cf = cfs;
            else
                Cf = cfs + (cfr - cfs) * (1 - 60 / Re_e);
            end
            %
        end

        %
        % Geometry optimisation
        %

        function [obj] = optimisation_geometry_objective(cmp,x,PR_ts_target)
            % It calculates the objective function
            %
            cmp = cmp.optimisation_geometry_assign(x,PR_ts_target);
            %
            cmp = cmp.simulation();
            %
            obj = -cmp.eta_is_ts;
            %
        end

        function [c,ceq] = optimisation_geometry_constraints(cmp,x,PR_ts_target)
            % It calculates the constraints
            %
            cmp = cmp.optimisation_geometry_assign(x,PR_ts_target);
            %
            cmp = cmp.simulation();
            %
            c = [];
            %
            ceq(1) = cmp.PR_ts - PR_ts_target;    % Designed to achieve a target PR_ts
            %
        end

        function [obj] = optimisation_geometry_scalarised_objective(cmp,x,PR_ts_target,eta_ts_target)
            %It calculates of the scalarised objective function to solve
            % a multi-objective problem according to (Wierzbicki, 1980) - 10.1007/978-3-642-48782-8_32
            %
            cmp = cmp.optimisation_geometry_assign(x,PR_ts_target);
            %
            cmp = cmp.simulation();
            %
            obj_1 = -cmp.eta_is_ts;
            obj_2 = -cmp.PR_ts;
            %
            obj_1_target = -eta_ts_target;
            obj_2_target = -PR_ts_target;
            %
            aug_factor = 1e-4;
            obj = max(obj_1 - obj_1_target, obj_2 - obj_2_target) + aug_factor * sum([obj_1 - obj_1_target, obj_2 - obj_2_target]);
            %
        end

        function [cmp] = optimisation_geometry_assign(cmp,x,PR_ts_target)
            %
            % Here is where the x components are assigned to the related geometry_specs entries
            %
            % Non optimised parameters
            %
            geometry_specs.N_blades = round(12.03 + 2.544 * PR_ts_target); %(Meroni et al., 2018) - 10.1016/j.apenergy.2018.09.210
            % geometry_specs.blade_thickness = 0.2 * 1e-3;
            geometry_specs.blade_thickness_inlet = 0.15 * 1e-3;
            geometry_specs.blade_thickness_outlet = 0.15 * 1e-3;
            geometry_specs.alpha1 = 0;
            geometry_specs.roughness = 2 * 1e-6;
            %
            % Optimised parameters
            %
            geometry_specs.r2 = x(1);
            geometry_specs.ratio_r1hub_r1shroud = x(2);
            geometry_specs.ratio_r1shroud_r2 = x(3);
            geometry_specs.ratio_b2_r2 = x(4);
            geometry_specs.beta2_blade = x(5);
            geometry_specs.beta1_blade_rms = x(6);
            geometry_specs.ratio_r3_r2 = x(7);
            %
            cmp = cmp.set_geometry(geometry_specs);
            %
        end

        %
        % Figures
        %

        function [f] = ellipse_internal(cmp,x,y)
            % It generates the points to draw the "shroud" line
            % (approximated as an ellipses)
            %
            f = x^2 / (cmp.axial_length - cmp.b2)^2 + (y - cmp.r2)^2 / (cmp.r1_shroud - cmp.r2)^2 - 1;
            %
        end

        function [f] = ellipse_external(cmp,x,y)
            % It generates the points to draw the "hub" line
            % (approximated as an ellipses)
            %
            f = x^2 / cmp.axial_length^2 + (y - cmp.r2)^2 / (cmp.r1_hub - cmp.r2)^2 - 1;
            %
        end

    end
end
%
%