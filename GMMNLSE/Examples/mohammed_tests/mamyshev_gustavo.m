% This code runs a ring Mamyshev oscillator with a rate-equation-gain model.
% 
% DUAL-ARM CONFIGURATION:
% - Arm 1: 6 µm core diameter, 4 m active + 2 m passive fiber
% - Arm 2: 20 µm core diameter, 3 m active + 1 m passive fiber
%
% Pulse path sequence:
%   1 m 6µm passive → 4 m 6µm active → 1 m 6µm passive → 1 m 20µm passive → 3 m 20µm active
%
% Total cavity length: 10 meters
% Active segments: 2 (4m 6µm + 3m 20µm active fiber)
% Passive segments: 3 (1m 6µm + 1m 6µm + 1m 20µm)
%
% The fiber_segments array defines each segment: {arm_id, core_diameter_um, length_m, is_active_boolean}

close all; clearvars;

addpath('../../GMMNLSE algorithm/','../../user_helpers/');

%% Gain info
% Please find details of all the parameters in "gain_info.m" if not specified here.
% Note that the use of single spatial mode is different from multi-spatial modes.
% "reuse_data" and "linear_oscillator_model" are activated and some related parameters are set.
gain_rate_eqn.cross_section_filename = 'Liekki Yb_AV_20160530.txt';
gain_rate_eqn.core_diameter = 6; % um
gain_rate_eqn.cladding_diameter = 125; % um
gain_rate_eqn.core_NA = 0.12;
% Set independent gain model parameters for each cavity arm
arm_gain_params = struct( ...
    'arm1', struct('absorption_wavelength_to_get_N_total', 920, 'absorption_to_get_N_total', 0.55, 'pump_wavelength', 976), ...
    'arm2', struct('absorption_wavelength_to_get_N_total', 920, 'absorption_to_get_N_total', 5.1, 'pump_wavelength', 976));
gain_rate_eqn.absorption_wavelength_to_get_N_total = arm_gain_params.arm2.absorption_wavelength_to_get_N_total; % default; overwritten per active segment by arm id
gain_rate_eqn.absorption_to_get_N_total = arm_gain_params.arm2.absorption_to_get_N_total; % default; overwritten per active segment by arm id
gain_rate_eqn.pump_wavelength = arm_gain_params.arm2.pump_wavelength; % default; overwritten per active segment by arm id
% Set independent copump powers for each cavity arm (W)
arm_copump_power = struct('arm1', 2, 'arm2', 3);
gain_rate_eqn.copump_power = arm_copump_power.arm2; % default; overwritten per active segment by arm id
gain_rate_eqn.counterpump_power = 0; % W
gain_rate_eqn.reuse_data = false; % For a ring or linear cavity, the pulse will enter a steady state eventually.
                                 % If reusing the pump and ASE data from the previous roundtrip, the convergence can be much faster, especially for counterpumping.
gain_rate_eqn.linear_oscillator = false; % For a linear oscillator, there are pulses from both directions simultaneously, which will deplete the gain;
                                         % therefore, the backward-propagating pulses need to be taken into account.
gain_rate_eqn.tau = 840e-6; % lifetime of Yb in F_(5/2) state (Paschotta et al., "Lifetme quenching in Yb-doped fibers"); in "s"
gain_rate_eqn.export_N2 = true; % whether to export N2, the ion density in the upper state or not
gain_rate_eqn.ignore_ASE = true;
gain_rate_eqn.max_iterations = 5; % If there is ASE, iterations are required.
gain_rate_eqn.tol = 1e-5; % the tolerance for the above iterations
gain_rate_eqn.verbose = true; % show the information(final pulse energy) during iterations of computing the gain

%% Setup fiber parameters
% load_default_GMMNLSE_propagate loads default parameters like 
%
% loading fiber.betas and fiber.SR based on your multimode folder above
% fiber.MFD = 6.2; % um
% sim.Raman_model = 1; Use isotropic Raman model
% sim.gain_model = 0; Don't use gain model = passive propagation
% sim.gpu_yes = true; Use GPU (default to true)
% ......

% General parameters
sim.lambda0 = 1030e-9;
%sim.progress_bar = false;
sim.save_period = 0.1;
sim.gpu_yes = false; % For only the fundamental mode, running with CPU is faster if the number of points is lower than 2^(~18).
sim.gain_model = 2;
sim.progress_bar_name = 'Gain';

%% Define fiber segments for dual-arm configuration
% Arm 1: 6 um core diameter (4 m active + 2 m passive = 6 m total)
% Arm 2: 20 um core diameter (3 m active + 1 m passive = 4 m total)
% Pulse path: 1 m 6um passive → 4 m 6um active → 1 m 6um passive → 1 m 20um passive → 3 m 20um active

% Define fiber segments: {arm_id, core_diameter (um), length (m), is_active}
fiber_segments = {
    % Arm 1: Start with 1m 6um passive (part of 2m passive section)
    {1, 6, 1, false},
    % Arm 1: 4m 6um active fiber
    {1, 6, 4, true},
    % Arm 1: end with 1m 6um passive
    {1, 6, 1, false},
    % Arm 2: 1m 20um passive, 3m 20um active
    {2, 20, 1, false},
    {2, 20, 3, true},
};

% Create fiber structures for each segment
num_segments = length(fiber_segments);
fiber_array = cell(num_segments, 1);
gain_array = cell(num_segments, 1);

for seg_idx = 1:num_segments
    arm_id = fiber_segments{seg_idx}{1};
    core_diam = fiber_segments{seg_idx}{2};
    fiber_length = fiber_segments{seg_idx}{3};
    is_active = fiber_segments{seg_idx}{4};
    
    % Create fiber structure for this segment
    fiber_temp = struct();
    fiber_temp.L0 = fiber_length;
    
    % Load fiber properties based on core diameter
    if core_diam == 6
        [fiber_temp, sim_temp] = load_default_GMMNLSE_propagate(fiber_temp, sim, 'single_mode');
    elseif core_diam == 20
        [fiber_temp, sim_temp] = load_default_GMMNLSE_propagate(fiber_temp, sim, 'single_mode');
    else
        error('Unsupported core diameter: %d um', core_diam);
    end

    % Keep global simulation settings consistent with loader-populated fields (e.g., sim.f0).
    sim = sim_temp;
    
    fiber_array{seg_idx} = fiber_temp;
    
    % Setup gain parameters for this segment
    if is_active
        gain_temp = gain_rate_eqn;
        gain_temp.core_diameter = core_diam;
        if arm_id == 1
            gain_temp.absorption_wavelength_to_get_N_total = arm_gain_params.arm1.absorption_wavelength_to_get_N_total;
            gain_temp.absorption_to_get_N_total = arm_gain_params.arm1.absorption_to_get_N_total;
            gain_temp.pump_wavelength = arm_gain_params.arm1.pump_wavelength;
            gain_temp.copump_power = arm_copump_power.arm1;
        elseif arm_id == 2
            gain_temp.absorption_wavelength_to_get_N_total = arm_gain_params.arm2.absorption_wavelength_to_get_N_total;
            gain_temp.absorption_to_get_N_total = arm_gain_params.arm2.absorption_to_get_N_total;
            gain_temp.pump_wavelength = arm_gain_params.arm2.pump_wavelength;
            gain_temp.copump_power = arm_copump_power.arm2;
        else
            error('Unsupported arm id: %d', arm_id);
        end
        % Store per-segment gain settings here; gain_info is initialized later after lambda is defined.
        gain_array{seg_idx} = gain_temp;
    else
        % For passive segments, create a structure with no gain
        gain_array{seg_idx} = struct('is_passive', true);
    end
end

fiber = fiber_array{1}; % Keep original structure for compatibility

%% Setup general cavity parameters
max_rt = 100;
N = 2^13; % the number of time/freq points
time_window = 50; % ps
dt = time_window/N;
f = sim.f0+(-N/2:N/2-1)'/(N*dt); % THz
t = (-N/2:N/2-1)'*dt; % ps
c = 299792458; % m/s
lambda = c./(f*1e12)*1e9; % nm
OC = 0.9; % output coupler
loss = 0.1; % coupling loss
tol_convergence = 1e-5;

% Calculate total cavity length
total_cavity_length = 0;
for i = 1:num_segments
    total_cavity_length = total_cavity_length + fiber_array{i}.L0;
end

%% Filter parameters
spectral_filter = struct('bw', 4, ...    % bandwidth (nm)
                         'cw', {1025, 1035}); % center wavelength (nm)

%% Setup initial conditions
tfwhm = 1; % ps
total_energy = 1; % nJ
pedestal_energy = 0.01; % nJ

phis = 0; % no extra spectral phase on the seed pulse
prop_output = build_noisy_MMgaussian(tfwhm, inf, time_window, total_energy, pedestal_energy, 1, N, sim.lambda0, phis, 0.01);

%% Saved field information
% Allocate larger save arrays to account for multiple segments
save_num = int64(total_cavity_length / sim.save_period + 10);
save_num = double(save_num);
saved_z = zeros(1, save_num);
splice_z = zeros(1, num_segments + 1);
splice_z(1) = 0;
for i = 1:num_segments
    splice_z(i+1) = splice_z(i) + fiber_array{i}.L0;
end
% Note: We'll add filter displacements after the 6um and 20um active sections.
filter_displacement = sim.save_period / 25;
output_field = zeros(N, 1, max_rt);
output_field2 = zeros(N, 1, max_rt);

%% Load gain parameters
L_air = 1; % 1 is the free-space length
c = 299792458; % m/s
v = 1/fiber_array{1}.betas(2)*1e12; % velocity in the fiber (use first segment as reference)
gain_rate_eqn.t_rep = total_cavity_length/v + L_air/c; % s (in s)

% Setup gain structures for all active segments
% Note: passive segments will use the is_passive flag instead of gain model
for i = 1:num_segments
    if fiber_segments{i}{4} % if is_active
        gain_temp = gain_rate_eqn;
        gain_temp.core_diameter = fiber_segments{i}{2};
        if fiber_segments{i}{1} == 1
            gain_temp.absorption_wavelength_to_get_N_total = arm_gain_params.arm1.absorption_wavelength_to_get_N_total;
            gain_temp.absorption_to_get_N_total = arm_gain_params.arm1.absorption_to_get_N_total;
            gain_temp.pump_wavelength = arm_gain_params.arm1.pump_wavelength;
            gain_temp.copump_power = arm_copump_power.arm1;
        elseif fiber_segments{i}{1} == 2
            gain_temp.absorption_wavelength_to_get_N_total = arm_gain_params.arm2.absorption_wavelength_to_get_N_total;
            gain_temp.absorption_to_get_N_total = arm_gain_params.arm2.absorption_to_get_N_total;
            gain_temp.pump_wavelength = arm_gain_params.arm2.pump_wavelength;
            gain_temp.copump_power = arm_copump_power.arm2;
        else
            error('Unsupported arm id: %d', fiber_segments{i}{1});
        end
        gain_array{i} = gain_info(fiber_array{i}, sim, gain_temp, ifftshift(lambda, 1));
    end
end

%% Run the cavity simulation
func = analyze_sim;

% Initialize tracking arrays for field evolution
field = zeros(N, 1, save_num, 'single');
N2 = zeros(1, 1, save_num, 'single');
pump = zeros(1, 1, save_num, 'single');

% Initialize some parameters
output_energy = zeros(max_rt,1);
rt_num = 0;
pulse_survives = true;
while rt_num < max_rt
    current_z = 0;
    time_delay = 0;
    zn = 1;
    saved_z_rt = zeros(1, save_num);
    rt_num = rt_num +1;
    
    t_iteration_start = tic;
    cprintf('*[1 0.5 0.31]','Iteration %d', rt_num);
    % -----------------------------------------------------------------
    for seg_idx = 1:num_segments
        current_fiber = fiber_array{seg_idx};
        current_segment = fiber_segments{seg_idx};
        core_diam = current_segment{2};
        is_active = current_segment{4};
        
        % Propagation inside fiber segments
        if is_active
            % Propagate through active fiber with gain
            prop_output = GMMNLSE_propagate(current_fiber, prop_output, sim, gain_array{seg_idx});
        else
            % Propagate through passive fiber (no gain)
            % Temporarily disable gain model
            sim_passive = sim;
            sim_passive.gain_model = 0; % No gain for passive fiber
            prop_output = GMMNLSE_propagate(current_fiber, prop_output, sim_passive);
        end
        
        time_delay = time_delay + prop_output.t_delay(end);

        % Save the field information
        num_saved = length(prop_output.z);
        field(:,:,zn:zn+num_saved-1) = prop_output.fields;
        saved_z_rt(zn:zn+num_saved-1) = current_z + prop_output.z;

        % Save the gain info if available
        if is_active && isfield(prop_output, 'N2')
            N2(:,:,zn:zn+num_saved-1) = prop_output.N2;
            pump(:,:,zn:zn+num_saved-1) = prop_output.Power.pump.forward;
        else
            % For passive sections, save zeros
            if ~exist('N2', 'var')
                N2 = zeros(1,1,save_num);
            end
            if ~exist('pump', 'var')
                pump = zeros(1,1,save_num);
            end
            N2(:,:,zn:zn+num_saved-1) = 0;
            pump(:,:,zn:zn+num_saved-1) = 0;
        end

        current_z = current_z + prop_output.z(end);
        zn = zn + num_saved - 1;
        
        % Spectral filter after the 6um active section (segment 2)
        if seg_idx == 2
            if rt_num ~= 1
                close(fig_filter);
            end
            [prop_output, fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter(1).cw, spectral_filter(1).bw, 3, true);

            prop_output.fields = prop_output.fields(:,:,end)*sqrt(1-loss);

            % Save the filtered field as its own cavity point.
            current_z = current_z + filter_displacement;
            zn = zn + 1;
            field(:,:,zn) = prop_output.fields;
            saved_z_rt(zn) = current_z;
            N2(:,:,zn) = N2(:,:,zn-1);
            pump(:,:,zn) = pump(:,:,zn-1);
        end

        % Output coupler and spectral filter after the 20um active section (segment 5)
        if seg_idx == 5  % After 3m 20um active fiber
            end_active_field = prop_output.fields(:,:,end);

            % Output coupler is placed immediately after the 20 um active fiber.
            output_field(:,:,rt_num) = sqrt(OC)*end_active_field;
            prop_output.fields = sqrt(1-OC)*end_active_field;

            output_field2(:,:,rt_num) = end_active_field;
            if rt_num ~= 1
                close(fig_filter);
            end
            [prop_output, fig_filter] = gaussian_spectral_filter(prop_output, sim.f0, spectral_filter(2).cw, spectral_filter(2).bw, 3, true);
            
            prop_output.fields = prop_output.fields(:,:,end)*sqrt(1-loss);
            
            % Save the filtered field as its own cavity point.
            current_z = current_z + filter_displacement;
            zn = zn + 1;
            field(:,:,zn) = prop_output.fields;
            saved_z_rt(zn) = current_z;
            N2(:,:,zn) = N2(:,:,zn-1);
            pump(:,:,zn) = pump(:,:,zn-1);
        end
    end
    saved_z = saved_z_rt(1:zn);
    
    % -----------------------------------------------------------------
    % Energy of the output field
    output_energy(rt_num) = sum(trapz(abs(output_field(:,:,rt_num)).^2))*prop_output.dt/1e3; % energy in nJ

    % If the energies stop changing, then we're done!
    if rt_num ~= 1
        close(fig);
    end
    warning('off')
    output_field(:,:,rt_num) = pulse_tracker(output_field(:,:,rt_num));
    warning('on');
    [converged_yes,fig] = check_convergence( output_energy,output_field(:,:,rt_num),f,t,tol_convergence,true );
    
    % ---------------------------------------------------------------------
    % Display running time
    t_iteration_end = toc(t_iteration_start);
    t_iteration_spent = datevec(t_iteration_end/3600/24);
    fprintf(': %1u:%2u:%3.1f\n',t_iteration_spent(4),t_iteration_spent(5),t_iteration_spent(6));
    
    % ---------------------------------------------------------------------
    % Update the repetition rate based on "time_delay"
    for seg_idx = 1:num_segments
        if fiber_segments{seg_idx}{4}  % if is_active
            gain_array{seg_idx}.t_rep = gain_rate_eqn.t_rep + time_delay*1e-12;
        end
    end
    % ---------------------------------------------------------------------
    % Plot
    if rt_num ~= 1
        close(fig_evolution);
    end
    fig_evolution = func.analyze_fields(t,f,field(:,:,1:zn),saved_z,splice_z);
    
    if rt_num ~= 1
        close(fig_gain);
    end
    pump_plot.forward = pump(:,:,1:zn);
    pump_plot.backward = zeros(1,1,length(saved_z));
    fig_gain = func.analyze_gain(saved_z,splice_z,pump_plot,N2(:,:,1:zn));
    
    % ---------------------------------------------------------------------
    % Break if converged
    if converged_yes
        cprintf('blue','The field has converged!\n');
        break;
    end
    % Break if pulse dies
    if output_energy(rt_num) < 0.01
        disp('The pulse dies.');
        pulse_survives = false;
        break;
    end
end

%% Finish the simulation and save the data
% Clear reducdant parts of the data
output_field = output_field(:,:,1:rt_num);
output_field2 = output_field2(:,:,1:rt_num);
energy = output_energy(arrayfun(@any,output_energy)); % clear zero

% -------------------------------------------------------------------------
% Save the final output field
save('ring_Mamyshev_oscillator.mat', 't','f','output_field','output_field2','time_delay','energy',...
                         'saved_z','splice_z','field',...
                         'N2','pump','fiber_array','gain_array',...
                         'fiber','sim','fiber_segments',... % cavity parameters
                         '-v7.3'); % saved mat file version
% -------------------------------------------------------------------------

close(fig,fig_filter,fig_evolution,fig_gain);