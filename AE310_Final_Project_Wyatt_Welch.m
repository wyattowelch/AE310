%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AE310 Final Project - Wyatt Welch, Nicolas Ngo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all

%% --------------------------------------------------------------------- %% Initilize
clc 

dat = load("NACA747A415_dat.txt");
xdat = dat(:,1);
ydat = dat(:,2);
xover1 = xdat(26:51);

ym1 = [-75, 0, 75, 75, 0, -75] / 100; 
zm1 = [16.775, 25.825, 14.625, -8.825, -11.05, -7.725] / 100;

AB = 84.51/100;
BC = 76.93/100;
CD = 23.45/100;
DE = 76.79/100;
EF = 75.27/100;
AF = 24.50/100;
BE = 36.88/100;
t1 = 00.25/100;
t2 = 00.60/100;
S1 = 05.50/100;
S2 = 04.50/100;
St = 01.00/100;

span = 30;
root = 2.5;
fuselage = 2;
engineW = 8e3;
w_TO = 800e3;
rho_fuel = 750;
g = 9.81;
a = 1.25 * g;
chord = 2.5;
E_skin = 71.7e9;
E_web = 68.9e9;
G_skin = 26e9;

%% --------------------------------------------------------------------- %% Wing Box Geometry
clc 

A_cs = [t1*AB, t1*BC, t1*EF, t1*DE, t2*AF, t2*BE, t2*CD, 8*St*S1, 8*St*S2, 8*St*S1, ...
    8*St*S2, 8*St*S1, 8*St*S2, 8*St*S1, 8*St*S2, 8*St*S2, 8*St*S1, 8*St*S2, 8*St*S1, ...
    8*St*S2, 8*St*S1, 8*St*S2, 8*St*S1];

y_bar = [AB/2, AB+BC/2, t2/2, EF, EF+DE-t2/2, EF/2, EF+DE/2, t2+S1/2, t2+St/2, ...
    AB-t2/2-S1/2, AB-t2/2-St/2, AB+t2/2+S1/2, AB+t2/2+St/2, AB+BC-t2/2-S1/2, ...
    AB+BC-t2/2-St/2, t2+St/2, t2+S1/2, EF-t2/2-St/2, EF-t2/2-S1/2, EF+t2/2+St/2, ...
    EF+t2/2+S1/2, EF+DE-t2/2-St/2, EF+DE-t2/2-S1/2];

z_bar = [t1+AF+t1/2, t1+BE+t1/2, t1+AF/2, t1+BE/2, t1+CD/2, t1/2, t1/2, t1+AF-St/2, ...
    t1+AF-St-S2/2, t1+BE-St/2, t1+BE-St-S2/2, t1+BE-St/2, t1+BE-St-S2/2, t1+CD-St/2, ...
    t1+CD-St-S2/2, t1+St+S2/2, t1+St/2, t1+St+S2/2, t1+St/2, t1+St+S2/2, t1+St/2, ...
    t1+St+S2/2, t1+St/2];

y_cen = sum(A_cs .* y_bar) / sum(A_cs);
z_cen = sum(A_cs .* z_bar) / sum(A_cs);

dy = y_bar - y_cen;
dz = z_bar - z_cen;

by = [AB, BC, t2, t2, t2, EF, DE, S1, St, S1, St, S1, St, S1, St, St, S1, St, ...
    S1, St, S1, St, S1];
hy = [t1, t1, AF, BE, CD, t1, t1, St, S2, St, S2, St, S2, St, S2, S2, St, S2, ...
    St, S2, St, S2, St];
hz = [AB, BC, t2, t2, t2, EF, DE, S1, St, S1, St, S1, St, S1, St, St, S1, St, ...
    S1, St, S1, St, S1];
bz = [t1, t1, AF, BE, CD, t1, t1, St, S2, St, S2, St, S2, St, S2, S2, St, S2, ...
    St, S2, St, S2, St];

Iy = sum((1/12) .* by .* (hy .^ 3) + A_cs .* (dz .^ 2));
Iz = sum((1/12) .* bz .* (hz .^ 3) + A_cs .* (dy .^ 2));
Iyz = sum(A_cs .* dy .* dz);

fprintf("Cross Sectional Area = %+.4E \nCentroid (y_bar, z_bar) = (%+.4E, %+.4E) " + ...
    "    \nArea Moments of Inertia (I_y, I_z, I_yz) = (%+.4E, %+.4E, %+.4E)\n\n", ...
    sum(A_cs), y_cen, z_cen, Iy, Iz, Iyz)

% PART A COMPLETE %
%% --------------------------------------------------------------------- %% Wing Loading Calculations
clc 

% Shrenk's
A_tot = 4563.25/10000;
A_fuel = A_tot - sum(A_cs);
x = linspace(0,span, 201);

FL_to = w_TO * 2.25;
FL_cr = w_TO;
Lift_half_tot = FL_cr;
% Internal Forces on Wn

w_to = w_TO; % total plane weight
w_fuel = rho_fuel * A_fuel * g * (span * .7);
rho_al = 2800;
w_wing = sum(A_cs) * span * rho_al;
w_eng = 8e3;

w_tot = w_fuel + w_wing + w_eng;
Pz = w_tot .* (sqrt(1 - .9 * (x ./ span) .^ 8)) .* (1 + .1 * (x ./ span));
Vz = trapz(x,Pz) - cumtrapz(x,Pz);
My = -(trapz(x, Vz) - cumtrapz(x, Vz));

Py = Pz * .007;
Vy = trapz(x,Py) - cumtrapz(x,Py);
Mz = trapz(x, Vy) - cumtrapz(x, Vy);

e_offset = (.25 * chord) - y_cen;
T = (trapz(x,Pz) - cumtrapz(x, Pz)) * e_offset;
Mx = Vz * e_offset;

V_max = max(abs(sqrt(Vy .^ 2 + Vz .^ 2)));
M_max = max(abs(sqrt(My .^ 2 + Mz .^ 2)));
T_max = max(abs(T));

fprintf("Maximum Bending Moment = %+.4E \nMaximum Shear Force = %+.4E " + ...
    "\nMaximum Torque = %+.4E", M_max, V_max, T_max)

% PART B COMPLETE%
%% --------------------------------------------------------------------- %% Deflection Analysis
clc 

alp1 = (1/E_skin) * (Iyz / (Iy * Iz - Iyz^2));
alp2 = (1/E_skin) * (Iy / (Iy * Iz - Iyz^2));

dvdx = cumtrapz(x, alp1 * My + alp2 * Mz);
v = cumtrapz(x, dvdx);
dwdx = cumtrapz(x, alp1 * Mz + alp2 * My);
w = cumtrapz(x, dwdx);
LovrT = [AB/t1, BC/t1, AF/t2, BE/t2, CD/t2, EF/t1, DE/t1];
J_inner = 4 * (A_fuel ^ 2) / sum(LovrT);
dThetaDx = Mx / (G_skin * J_inner);
theta = cumtrapz(x, dThetaDx);

% PART C COMPLETE, POSSIBLY NO AXIAL LOAD %
%% --------------------------------------------------------------------- %% Stress Analysis
clc 

Mx0 = Mx(1);
My0 = My(1);
Mz0 = Mz(1);
K1 = (My0 * Iyz + Mz0 * Iy) / (Iy * Iz - Iyz^2);
K2 = (Mz0 * Iyz + My0 * Iz) / (Iy * Iz - Iyz^2);

Ymid = (ym1 + ym1([2:6,1]))/2;
Zmid = (zm1 + zm1([2:6,1]))/2;
checkpoints = [ym1, zm1; Ymid, Zmid]';

sigma = -(checkpoints(:,1) - y_cen) .* K1 + K2 .* (checkpoints(:,2) - z_cen);
[sigma_max_ten, idx_t] = max(sigma);
[sigma_max_com, idx_c] = min(sigma);
max_ten_loc = checkpoints(idx_t);
max_com_loc = checkpoints(idx_c);

Vy0 = Vy(1);
Vz0 = Vz(1);

L_skin = [AB BC AF CD EF DE];
t_skin = [t1 t1 t2 t2 t1 t1];
A_skin = L_skin .* t_skin;
A_stringer = .001;

A_corner = A_stringer + (A_skin + circshift(A_skin, 1)) / 6;
A_mid = (2/3) * A_skin;
A_lump = [A_corner, A_mid];

Qz = A_lump .* (checkpoints(:,1)' - y_cen);
Qy = A_lump .* (checkpoints(:,2)' - z_cen);

q_set = Vz0 / Iy * Qz + Vy0 / Iz * Qy;
q0 = -sum(q_set(7:12) .* (L_skin ./ t_skin)) / sum(L_skin ./ t_skin);
q_shear = q_set + q0;

q_T = Mx0 / (2 * A_fuel);

q_tot = q_shear + q_T;
tau = q_tot ./ [t1 t1 t1 t1 t1 t1 t1 t2 t2 t2 t1 t1];
[tau_max, idx_max] = max(abs(tau));
max_shear = q_tot(idx_max);
max_shear_loc = (idx_max);

fprintf("Maximum Bending Stresses: Tension = %+.4E, Compression = %+.4E " + ...
    "\nMaximum Bending Locations: Tension = %+.4E, Compression = %+.4E " + ...
    "\nMax Shear = %+.4E, Location = %+.4E\n",sigma_max_ten,sigma_max_com,max_ten_loc,max_com_loc,max_shear,max_shear_loc) 
disp("A_lump:")
disp(A_lump')
disp("q_shear:")
disp(q_shear')

% PART D COMPLETE %
%% --------------------------------------------------------------------- %% Failure Analysis
clc 

Sy = 505e6;
sigma_v = sqrt(sigma' .^ 2 + 3 * (tau .^ 2));
[sigma_v_max, idx_vm] = max(sigma_v);
sigma_loc = checkpoints(idx_vm);
FS = (.8 * Sy) ./ sigma_v_max;

nu = .33;
sigma_cr = (pi ^ 2 * E_skin * ((t1/AB)^2)) / (12 * (1 - nu ^ 2));
Ks = 5.34 + 4 * (chord/AB) ^ 2;
tau_cr = sigma_cr * Ks;

sigma_bend = sigma(7);
tau_shear = tau(7);
lambda = sqrt((sigma_cr / sigma_bend) ^ 2 + (tau_cr / tau_shear) ^ 2);
M_buckle = (lambda - 1) * 100;

fprintf('Max Von Mises = %+.4E, Located at (y, z) = (%0.3f, %0.3f), FS = %.2f\n', ...
        sigma_v_max, sigma_loc, sigma_loc, FS);
fprintf("Buckling Margin = %4.2f, on panel AB's midpoint\n", M_buckle);

%% --------------------------------------------------------------------- %% Plots
clc 

figure(1) % Lift Distribution per Span
hold on, grid on
plot(x, Pz)
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Lift (N/m)', 'Interpreter','latex')
title('Lift Distribution per Span', 'Interpreter','latex')


figure(2) % Shear Forces in y and z
subplot(2,1,2)
plot(x, Vy)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Transverse Shear (N)', 'Interpreter','latex')
title('Transverse Shear in y Axis', 'Interpreter','latex')

subplot(2,1,1)
plot(x, Vz)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Transverse Shear (N)', 'Interpreter','latex')
title('Transverse Shear in z Axis', 'Interpreter','latex')


figure(3) % Bending Moments in y and z
subplot(3,1,3)
plot(x, Mx)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Bending Moment (N/m)', 'Interpreter','latex')
title('Bending Moment in x Axis', 'Interpreter','latex')

subplot(3,1,2)
plot(x, My)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Bending Moment (N/m)', 'Interpreter','latex')
title('Bending Moment in y Axis', 'Interpreter','latex')

subplot(3,1,1)
plot(x, Mz)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Bending Moment (N/m)', 'Interpreter','latex')
title('Bending Moment in z Axis', 'Interpreter','latex')


figure(4) % Torque
plot(x, T)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('Torque', 'Interpreter','latex')
title('Torue Distribution per Span', 'Interpreter','latex')


figure(5) % Bending Deflections
subplot(2,1,1)
plot(x,v)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('v(x) (m)', 'Interpreter','latex')
title('Deflection v in y Direction', 'Interpreter','latex')

subplot(2,1,2)
plot(x,w)
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('w(x) (m)', 'Interpreter','latex')
title('Deflection w in z Direction', 'Interpreter','latex')


figure(6)
plot(x, theta * (180/pi))
hold on, grid on
xlabel("Span (m)", 'Interpreter','latex')
ylabel('$\theta$(x) (m)', 'Interpreter','latex')
title('Twist Angle Distribution', 'Interpreter','latex')


%% --------------------------------------------------------------------- %% End
clc, clear all