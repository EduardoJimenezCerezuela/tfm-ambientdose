% =========================================================================
%  Dosis.m 
%  Autor: Eduardo Jiménez Cerezuela
% =========================================================================
%
%  1. Genera una función del espectro parametrizada con valores arbitrarios
%  y se visualiza con el espectro real. Saca el chi^2 entre ambas funciones.
%
%  2. A partir de esos valores, se pregunta al usuario si se desea aumentar 
%  porcentualmente el valor integral de cada subfunción del espectro. De 
%  esta forma, podemos analizar qué tan bueno es nuestro espectro modelo* si
%  lo asumimos como real*. Visualizamos, por tanto, el espectro real; el 
%  modelo (parametrizado); y, el modelo*/real* (parametrizado con nuevos pesos).
%
%  3. Calcula la dosis estiamda (H*) como la c.l. de las funciones de respuesta 
%  (tasas de conteo) con el espectro real*. Además, saca la dosis calculada 
%  con los coeficientes h* y las compara con el error relativo.
%
%  4. Genera casos sobre el espectro parametrizado y analiza la correlación
%  entre la dosis estimada y calculada. Aplica correciones a los
%  coeficientes del ICRP 74 (e ICRU si se desea) para bajar el error
%  relativo
%
% =========================================================================

clear all; close all; clc; 
addpath('C:\Users\eduji\Desktop\Master\TFM\Data');

%% LECTURA DE DATOS
[datos.vEphi,  datos.vPhi ]   = read_2cols('EXPACSJavalambre_phi.txt'); % Espectro
[datos.vEEphi, datos.vEPhi]  = read_2cols('EXPACSJavalambre_Ephi.txt'); % Espectro de letargia
[datos.vER, datos.vRbare, datos.vRcore, datos.vRhepb] = read_4cols('Response_Functions_HENSApp_OAJ_Phase0.dat'); % Funciones de respuesta
[datos.vEh95,   datos.vh95]       = read_2cols('h_Emax_ICRU95'); % Coeficientes h*
[datos.vEh1074,   datos.vh1074]       = read_2cols('h10_ICRP74'); % Coeficientes h* (10)
[datos.vEh1074ferrari,   datos.vh1074ferrari]       = read_2cols('h10_extension_ferrarietal1997'); % Extensión coeficientes h* (10) ferrari
[datos.vEh1074sannikov,   datos.vh1074sannikov]       = read_2cols('h10_extension_sannikovetal1997'); % Extensión coeficientes h* (10) sannikov

% Unimos ICRP 74 a cada extensión, manteniendo solo las filas donde la 
% energía >= 500 MeV o 300 (depende de la extension)
% Lectura Ferrari
datos.vh1074ferrari = datos.vh1074ferrari(datos.vEh1074ferrari(:,1) >= 500, :);
datos.vEh1074ferrari = datos.vEh1074ferrari(datos.vEh1074ferrari(:,1) >= 500, :);
% concatenar Ferrari
datos.vh1074ferrari = [datos.vh1074; datos.vh1074ferrari];
datos.vEh1074ferrari = [datos.vEh1074; datos.vEh1074ferrari];

% Lectura Sannikov
datos.vh1074sannikov = datos.vh1074sannikov(datos.vEh1074sannikov(:,1) >= 300, :);
datos.vEh1074sannikov = datos.vEh1074sannikov(datos.vEh1074sannikov(:,1) >= 300, :);
% concatenar Sannikov
datos.vh1074sannikov = [datos.vh1074; datos.vh1074sannikov];
datos.vEh1074sannikov = [datos.vEh1074; datos.vEh1074sannikov];


% Visualización de los 3 tipos de coeficientes de fluencia-dosis
figure;
loglog(datos.vEh95,datos.vh95,'r', 'DisplayName', 'ICRU 95', 'LineWidth', 2); hold on;
loglog(datos.vEh1074ferrari,datos.vh1074ferrari,'g', 'DisplayName', 'ICRP 74 F', 'LineWidth', 2);
loglog(datos.vEh1074sannikov,datos.vh1074sannikov,'b', 'DisplayName', 'ICRP 74 S', 'LineWidth', 2);
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
legend;
xlabel('Energía (MeV)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('h^* / h^*(10) (pSv cm²)', 'FontSize', 18, 'FontWeight', 'bold');
title('Coeficientes de conversión fluencia-dosis', 'FontSize', 18, 'FontWeight', 'bold');

fprintf(['[LECTURA] EXPACS: %d pts | (ICRU 95) h*(E): %d pts | (ICRP 74 - Ferrari) h*(10): %d pts |' ...
    ' \n         | (ICRP 74 - Sannikov) h*(10): %d pts| Respuesta: %d pts | \n\n'], ...
        length(datos.vEphi), length(datos.vEh95), length(datos.vEh1074ferrari), length(datos.vEh1074sannikov), length(datos.vER));

%% MALLA DE ENERGÍA COMÚN PARA CÁLCULOS E INTEGRALES
Emin   = max([datos.vEphi(1), datos.vEEphi(1), datos.vEh95(1), datos.vEh1074ferrari(1), datos.vEh1074sannikov(1) datos.vER(1)]);
Emax   = min([datos.vEphi(end), datos.vEEphi(end), datos.vEh95(end), datos.vER(end)]);
E_plot = logspace(log10(Emin), log10(100000), 700);

fprintf('[LECTURA] %d puntos de malla E activos en [%.2e, %.2e] MeV\n\n', ...
        length(E_plot), E_plot(1), E_plot(end));

%% INTERPOLACIÓN LOG-LOG DE LOS DATOS
% Interpolación de tipo PCHIP (Piecewise Cubic Hermite Interpolating Polynomial).
%
%   Ev: Energía (x)
%   Fv: valores correpondientes a Ev (y)
%   Eq: Los puntos que queremos consultar (Los bines) - Los guardamos en
%       una estructura
%
% =========================================================================

interp_ll = @(Ev, Fv, Eq) ...
    10.^interp1(log10(Ev), log10(Fv), log10(Eq), 'pchip'); % Representación log-log 

datos.vPhi_interp = interp_ll(datos.vEphi, datos.vPhi, E_plot);
datos.vEPhi_interp = interp_ll(datos.vEEphi, datos.vEPhi, E_plot);
datos.vh95_interp = interp_ll(datos.vEh95, datos.vh95, E_plot); % Coeficientes h*
datos.vh1074ferrari_interp = interp_ll(datos.vEh1074ferrari, datos.vh1074ferrari, E_plot); % Coeficientes h*(10) ferrari
datos.vh1074sannikov_interp = interp_ll(datos.vEh1074sannikov, datos.vh1074sannikov, E_plot); % Coeficientes h*(10) sannikov
datos.vRbare_interp = interp_ll(datos.vER, datos.vRbare, E_plot); % Respuesta 1
datos.vRcore_interp = interp_ll(datos.vER, datos.vRcore, E_plot); % Respuesta 2
datos.vRhepb_interp = interp_ll(datos.vER, datos.vRhepb, E_plot); % Respuesta 3

%% FUNCIÓN PARAMETRIZADA
% Coeficientes: AJUSTAR A MANO
coef.k  = 4.209e-02;  % factor k  – escala global
coef.p1 = 0.1944;     % p1– peso térmico (%)
coef.To = 2.513e-8;   % To – temp. térmica (MeV)
coef.p2 = 0.2366;     % p2 – peso epitérmico (%)
coef.b  = 0.0503;     % b - pendiente epitérmica
coef.bp = 2.0;        % bp – corte epitérmico (MeV)
coef.Ed = 1.723e-7;   % Ed – corte Cd (MeV)
coef.p3 = 0.3162;     % p3 – peso evaporación (%)
coef.bf = 0.2389;     % bf – beta
coef.gf = -0.5919;    % gf – gamma
coef.p4 = 0.2528;     % p4 – peso alta energía (%)
coef.bh = 0.6261;     % bh – beta
coef.gh = 5.0;        % gh – gamma

%% 1. REPRESENTACIÓN PHI EXPACS CON PHI PARAMÉTRICA + CHI^2

fprintf('------ INICIO DE LA SECCIÓN 1 ------ \n\n');

% Obtenemos el espectro del modelo
% Aseguramos que los vectores tengan la misma orientación con (:)
[phiparam.tot, phiparam.t, phiparam.e, phiparam.f, phiparam.h] = compute_phi(E_plot, coef);
phiparam.E_phiparam_tot = E_plot(:)' .* phiparam.tot(:)';
phiparam.E_phiparam_t   = E_plot(:)' .* phiparam.t(:)';
phiparam.E_phiparam_e   = E_plot(:)' .* phiparam.e(:)';
phiparam.E_phiparam_f   = E_plot(:)' .* phiparam.f(:)';
phiparam.E_phiparam_h   = E_plot(:)' .* phiparam.h(:)';

fprintf('[ANÁLISIS] Componente Térmica: %.1f %% \n', coef.p1.*100);
fprintf('[ANÁLISIS] Componente Epitérmica: %.1f %% \n', coef.p2.*100);
fprintf('[ANÁLISIS] Componente Evaporación: %.1f %% \n', coef.p3.*100);
fprintf('[ANÁLISIS] Componente Alta Energía: %.1f %% \n', coef.p4.*100);

% Cálculo de Chi^2 + protección para evitar división por cero
chi_cuadrado = sum( ((datos.vEPhi_interp - phiparam.E_phiparam_tot) ./ (0.01 .* datos.vEPhi_interp)).^2 );
fprintf('[ANÁLISIS] Chi^2: %.4e\n', chi_cuadrado);
n_params = numel(fieldnames(coef));
fprintf('[ANÁLISIS] Chi^2_r (reducido): %.4e\n\n', chi_cuadrado/(length(E_plot)-n_params)); % Chi^2 reducido

% Representación Gráfica
fprintf('[FIGURA] Comparativa de espectros de Letargia\n\n');
figure('Name', 'Espectro de Letargia', 'Color', 'w', 'Position', [80 200 950 550]);

% Espectro EXPACS
semilogx(E_plot, datos.vEPhi_interp, 'r.', 'MarkerSize', 8, 'DisplayName', 'Espectro EXPACS');
hold on;

% Espectro Paramétrico Total
semilogx(E_plot, phiparam.E_phiparam_tot, 'k-', 'LineWidth', 2.5, 'DisplayName', ...
    sprintf('Modelo Parametrizado'));

% Subfunciones Paramétricas (Líneas discontinuas)
semilogx(E_plot, phiparam.E_phiparam_t, 'b--', 'LineWidth', 1.2, ...
'DisplayName', sprintf('Componente Térmica'));
semilogx(E_plot, phiparam.E_phiparam_e, 'g--', 'LineWidth', 1.2, ...
'DisplayName', sprintf('Componente Epitérmica'));
semilogx(E_plot, phiparam.E_phiparam_f, 'm--', 'LineWidth', 1.2, ...
'DisplayName', sprintf('Componente Evaporación'));
semilogx(E_plot, phiparam.E_phiparam_h, 'c--', 'LineWidth', 1.2, ...
'DisplayName', sprintf('Componente Alta Energía'));

% Líneas para los rangos de integración
xline(3.5e-7, ':', '0.35 eV', 'Color', '#80B3FF', 'LineWidth', 2, 'HandleVisibility', 'off', 'FontWeight', 'bold', 'FontSize', 13);
xline(0.1, ':', '0.1 MeV', 'Color', '#80B3FF', 'LineWidth', 2, 'HandleVisibility', 'off', 'FontWeight', 'bold', 'FontSize', 13);
xline(20, ':', '20 MeV', 'Color', '#80B3FF', 'LineWidth', 2, 'HandleVisibility', 'off', 'FontWeight', 'bold','FontSize', 13);

% Configurar propiedades de los ejes
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
grid on; grid minor;

xlabel('Energía (MeV)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('E · φ̇ (E) (cm⁻² s⁻¹)', 'FontSize', 18, 'FontWeight', 'bold');
title('Espectro de Letargia: EXPACS vs. Modelo Parametrizado', 'FontSize', 18, 'FontWeight', 'bold');

legend('Location', 'best', 'FontSize', 13);
xlim([E_plot(1) E_plot(end)]);

%% REPRESENTACIÓN DEL ERROR RELATIVO EXPACS vs. MODELO PARAMETRIZADO
fprintf('[FIGURA] Error relativo por zonas energéticas\n\n');

% Error relativo (%)
err_rel = 100 .* (datos.vEPhi_interp - phiparam.E_phiparam_tot) ./ datos.vEPhi_interp;

% ── Línea de error cero ────────────────────────────────────────────────────
yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
figure;
% ── Error relativo por zona ───────────────────────────────────────────────
semilogx(E_plot, err_rel, 'b-',  'LineWidth', 1.8, ...
    'DisplayName', sprintf('µ (media) = %+.1f%%, σ (std) = %.1f%%', mean(err_rel), std(err_rel)));

% ── Formato ───────────────────────────────────────────────────────────────
set(gca, 'XScale', 'log');
xlabel('Energía (MeV)', 'FontWeight', 'bold');
ylabel('Error relativo (%)', 'FontWeight', 'bold');
title('Error relativo  E · φ̇ (E): EXPACS vs. Modelo Parametrizado', 'FontSize', 12);
grid on; grid minor;
legend('Location', 'best', 'FontSize', 9);
xlim([E_plot(1) E_plot(end)]);
ylim([-100 100]);

%% 2. ANÁLISIS DE INTENSIDAD (ESPECTRO PARAMETRIZADO)
%  En esta sección "inyectamos" neutrones adicionales en cada zona, de tal 
%  forma que al aumentar la fluencia en una zona, mantenemos el área 
%  absoluta de las otras, a cambio de variar el factor k.
%
%  Matemáticamente:
%
%  k'   = k · ( 1 + sum(p_i · Delta_i) )
%  p_i' = p_i · ( 1 + Delta_i / 1 + sum(p_i · Delta_i) )
%
% =========================================================================
fprintf('------ INICIO DE LA SECCIÓN 2 ------ \n\n');

fprintf('[LECTURA] VARIACIÓN DEL ESPECTRO PARAMETIZADO (Perturbación de Intensidad) ---\n');
fprintf('[LECTURA] Introduce el incremento PORCENTUAL de cada zona\n');
fprintf('          Pulsar ENTER asume 0%%. Valores a partir de -99.9%%.\n\n');

% Captura de Inputs con valores por defecto (0) ---
prompt_analisis = '> Incremento %% en zona %s: ';

val = input(sprintf(prompt_analisis, 'TÉRMICA (p1)'));
if isempty(val), val = 0; end
inc_p1 = max(val, -99.9) / 100; % Seguridad: no permitimos <= -100%

val = input(sprintf(prompt_analisis, 'EPITÉRMICA (p2)'));
if isempty(val), val = 0; end
inc_p2 = max(val, -99.9) / 100;

val = input(sprintf(prompt_analisis, 'RÁPIDA (p3)'));
if isempty(val), val = 0; end
inc_p3 = max(val, -99.9) / 100;

val = input(sprintf(prompt_analisis, 'ALTA ENERGÍA (p4)'));
if isempty(val), val = 0; end
inc_p4 = max(val, -99.9) / 100;

% Definimos los vectores y guardamos los incremento
p_viejos = [coef.p1, coef.p2, coef.p3, coef.p4];
vector_incremento = [inc_p1, inc_p2, inc_p3, inc_p4];

% El factor de crecimiento total (1 + sum(p_i · Delta_i))
factor_crecimiento = 1 + sum(p_viejos .* vector_incremento);

coef_star = coef; % Copiamos los parámetros para sobreescribir los p viejos

% Calculamos p_i' = p_i · (1 + Delta_i / 1 + sum(p_i · Delta_i)) =
%                 = (Area_vieja * (1 + incremento)) / Area_Total_Nueva
coef_star.p1 = (coef.p1 * (1 + inc_p1)) / factor_crecimiento;
coef_star.p2 = (coef.p2 * (1 + inc_p2)) / factor_crecimiento;
coef_star.p3 = (coef.p3 * (1 + inc_p3)) / factor_crecimiento;
coef_star.p4 = (coef.p4 * (1 + inc_p4)) / factor_crecimiento;

coef_star.k  = coef.k * factor_crecimiento; % Nuevo factor k 

% Nuevo Espectro Parametrizado
[phiparam.phistar.tot, phiparam.phistar.t, phiparam.phistar.e, phiparam.phistar.f, phiparam.phistar.h] = compute_phi(E_plot, coef_star);
phiparam.phistar.E_tot = E_plot(:)' .* phiparam.phistar.tot(:)';

%% 3. DOSIS
fprintf('\n\n------ INICIO DE LA SECCIÓN 3 ------ \n\n');

% Dosis h95 calculada
    h95_interp = datos.vh95_interp(:)'; 
        dosis_95_EXPACS = trapz(E_plot, (datos.vPhi_interp(:)') .* h95_interp);
        dosis_95_PARAM  = trapz(E_plot, phiparam.tot(:)' .* h95_interp);
        dosis_95_PARAM2 = trapz(E_plot, phiparam.phistar.tot(:)' .* h95_interp);

% Dosis h95 estimada
    coef.a1 = 0;
    coef.a2 = 0.983010270763716;
    coef.a3 = 13.332602841589240;
        dosis_CL95_EXPACS = trapz(E_plot, datos.vPhi_interp.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));
        dosis_CL95_PARAM  = trapz(E_plot, phiparam.tot.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));
        dosis_CL95_PARAM2 = trapz(E_plot, phiparam.phistar.tot.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));

% Dosis h74 Ferrari calculada
    h74F_interp = datos.vh1074ferrari_interp(:)'; 
        dosis_74F_EXPACS = trapz(E_plot, (datos.vPhi_interp(:)') .* h74F_interp);
        dosis_74F_PARAM  = trapz(E_plot, phiparam.tot(:)' .* h74F_interp);
        dosis_74F_PARAM2 = trapz(E_plot, phiparam.phistar.tot(:)' .* h74F_interp);

% Dosis h74 Ferrari estimada
    coef.a1 = 0;
    coef.a2 = 9.592616973362901;
    coef.a3 = 8.118038494845630;
        dosis_CL74F_EXPACS = trapz(E_plot, datos.vPhi_interp.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));
        dosis_CL74F_PARAM  = trapz(E_plot, phiparam.tot.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));
        dosis_CL74F_PARAM2 = trapz(E_plot, phiparam.phistar.tot.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));

% Dosis h74 Sannikov calculada
    h74S_interp = datos.vh1074sannikov_interp(:)'; 
        dosis_74S_EXPACS = trapz(E_plot, (datos.vPhi_interp(:)') .* h74S_interp);
        dosis_74S_PARAM  = trapz(E_plot, phiparam.tot(:)' .* h74S_interp);
        dosis_74S_PARAM2 = trapz(E_plot, phiparam.phistar.tot(:)' .* h74S_interp);

% Dosis CL h74 Sannikov estimada
    coef.a1 = 0;
    coef.a2 = 9.239092852385420;
    coef.a3 = 8.603464735268430;
        dosis_CL74S_EXPACS = trapz(E_plot, datos.vPhi_interp.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));
        dosis_CL74S_PARAM  = trapz(E_plot, phiparam.tot.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));
        dosis_CL74S_PARAM2 = trapz(E_plot, phiparam.phistar.tot.* ...
            (datos.vRbare_interp*coef.a1 + datos.vRcore_interp*coef.a2 + datos.vRhepb_interp*coef.a3));

fprintf('\n==================================================\n');
fprintf('        TABLA DE DOSIS Y ERRORES RELATIVOS  \n');
fprintf('==================================================\n');
fprintf(' H* 95 EXPACS : %.4e pSv/s\n', dosis_95_EXPACS);
fprintf(' H* 95 PARAM: %.4e pSv/s\n', dosis_95_PARAM);
fprintf(' H* 95 PARAM mod.: %.4e pSv/s\n', dosis_95_PARAM2);
fprintf(' H* CL EXPACS : %.4e pSv/s\n', dosis_CL95_EXPACS);
fprintf(' H* CL PARAM : %.4e pSv/s\n', dosis_CL95_PARAM);
fprintf(' H* CL PARAM mod.: %.4e pSv/s\n', dosis_CL95_PARAM2);
fprintf('==================================================\n');

fprintf(' H* 74F EXPACS : %.4e pSv/s\n', dosis_74F_EXPACS);
fprintf(' H* 74F PARAM: %.4e pSv/s\n', dosis_74F_PARAM);
fprintf(' H* 74F PARAM mod.: %.4e pSv/s\n', dosis_74F_PARAM2);
fprintf(' H* CL EXPACS : %.4e pSv/s\n', dosis_CL74F_EXPACS);
fprintf(' H* CL PARAM : %.4e pSv/s\n', dosis_CL74F_PARAM);
fprintf(' H* CL PARAM mod.: %.4e pSv/s\n', dosis_CL74F_PARAM2);
fprintf('==================================================\n');

fprintf(' H* 74S EXPACS : %.4e pSv/s\n', dosis_74S_EXPACS);
fprintf(' H* 74S PARAM: %.4e pSv/s\n', dosis_74S_PARAM);
fprintf(' H* 74S PARAM mod.: %.4e pSv/s\n', dosis_74S_PARAM2);
fprintf(' H* CL EXPACS : %.4e pSv/s\n', dosis_CL74S_EXPACS);
fprintf(' H* CL PARAM : %.4e pSv/s\n', dosis_CL74S_PARAM);
fprintf(' H* CL PARAM mod.: %.4e pSv/s\n', dosis_CL74S_PARAM2);
fprintf('==================================================\n');

fprintf('[FIGURA] Comparativa de espectros \n\n');
figure('Name', 'Análisis de Intensidad', 'Color', 'w', 'Position', [100 100 1000 600]);

semilogx(E_plot, phiparam.E_phiparam_tot, 'k-', 'LineWidth', 2, ...
    'DisplayName', sprintf('Modelo Parametrizado'));
hold on;
semilogx(E_plot, phiparam.phistar.E_tot, 'b--', 'LineWidth', 2, ...
    'DisplayName', sprintf('Modelo Parametrizado modificado'));

set(gca, 'FontSize', 13, 'FontWeight', 'bold');

% Líneas para los rangos de integración
xline(3.5e-7, ':', '0.35 eV', 'Color', '#80B3FF', 'LineWidth', 2, 'HandleVisibility', 'off', 'FontSize', 13, 'FontWeight', 'bold');
xline(0.1, ':', '0.1 MeV', 'Color', '#80B3FF', 'LineWidth', 2, 'HandleVisibility', 'off', 'FontSize', 13, 'FontWeight', 'bold');
xline(20, ':', '20 MeV', 'Color', '#80B3FF', 'LineWidth', 2, 'HandleVisibility', 'off', 'FontSize', 13, 'FontWeight', 'bold');

xlabel('Energía (MeV)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('E · φ̇ (E) (cm⁻² s⁻¹)', 'FontSize', 18, 'FontWeight', 'bold');
title('Análisis de Sensibilidad de Intensidad', 'FontSize', 18, 'FontWeight', 'bold');


grid on; 
legend('show', 'Position', [0.3, 0.45, 0.18, 0.15]);
xlim([E_plot(1) E_plot(end)]);

% Etiquetas de variación de Intensidad en cada zona
% Calculamos el % 
var_p1 = inc_p1 * 100;
var_p2 = inc_p2 * 100;
var_p3 = inc_p3 * 100;
var_p4 = inc_p4 * 100;

% Límites del eje X e Y de la gráfica actual
lim_x = xlim; 
x_inicio = lim_x(1); 
x_fin = lim_x(2);
lim_y = ylim;
y_pos = lim_y(2) * 0.92; % Altura: al 92% del tope de la gráfica

% Fronteras de las zonas 
z1 = 3.5e-7;
z2 = 0.1;
z3 = 20;

% Centros de cada zona (raíz cuadrada porque el eje X es logarítmico)
centro_z1 = sqrt(x_inicio * z1);
centro_z2 = sqrt(z1 * z2);
centro_z3 = sqrt(z2 * z3);
centro_z4 = sqrt(z3 * x_fin);

% Colocamos los textos centrados con formato de recuadro
props = {'HorizontalAlignment', 'center', 'BackgroundColor', 'w', ...
         'EdgeColor', '#80B3FF', 'Margin', 3, 'FontSize', 12, 'FontWeight', 'bold'};
text(centro_z1, y_pos, sprintf('\\delta: %+.2f%%', var_p1), props{:});
text(centro_z2, y_pos, sprintf('\\delta: %+.2f%%', var_p2), props{:});
text(centro_z3, y_pos, sprintf('\\delta: %+.2f%%', var_p3), props{:});
text(centro_z4, y_pos, sprintf('\\delta: %+.2f%%', var_p4), props{:});

%% 4. ESTUDIO DE CORRELACIÓN Y CALIBRACIÓN INDEPENDIENTE
fprintf('\n\n------ INICIO DE LA SECCIÓN 4 ------ \n\n');

% Definición de los 21 casos de la tabla [dT, dE, dF, dH]
casos_variacion = [
     0,   0,   0,   0;   
     10,   0,   0,   0;    
     0,  10,   0,   0;    
     0,   0,  10,   0;    
     0,   0,   0,  10;
     30,   0,   0,   0;    
     0,  30,   0,   0;    
     0,   0,  30,   0;    
     0,   0,   0,  30;   
     50,  10, -20,   0;
    -20,   0,  40, -10;  
    -10,   0,   0,   0;    
      0, -10,   0,   0;    
      0,   0, -10,   0;    
      0,   0,   0, -10;
    -30,   0,   0,   0;    
      0, -30,   0,   0;    
      0,   0, -30,   0;    
      0,   0,   0, -30;    
      0,   0,  30, -10;
      0,   0, -30,  20
];
num_casos = size(casos_variacion, 1);


%% 4A. ICRP 74F (FERRARI)
fprintf('\n\n------ INICIO DEL ANÁLISIS AISLADO: ICRP 74F (FERRARI) ------ \n\n');


% Inicializar vectores de dosis para Ferrari
H_calc_F = zeros(num_casos, 1); 
H_est_F  = zeros(num_casos, 1);

% Combinación lineal original para Ferrari
Respuesta_CL_F = datos.vRbare_interp*0 + datos.vRcore_interp*9.592616973362901 + datos.vRhepb_interp*8.118038494845630;
h74F_interp = datos.vh1074ferrari_interp(:)'; 
p_viejos = [coef.p1, coef.p2, coef.p3, coef.p4];

% Bucle de cálculo para Ferrari
for i = 1:num_casos
    inc = casos_variacion(i, :) / 100;
    factor_crecimiento = 1 + sum(p_viejos .* inc);
    c_tmp = coef; 
    c_tmp.p1 = (coef.p1 * (1 + inc(1))) / factor_crecimiento;
    c_tmp.p2 = (coef.p2 * (1 + inc(2))) / factor_crecimiento;
    c_tmp.p3 = (coef.p3 * (1 + inc(3))) / factor_crecimiento;
    c_tmp.p4 = (coef.p4 * (1 + inc(4))) / factor_crecimiento;
    c_tmp.k  = coef.k * factor_crecimiento;
    
    [Phi_tmp_tot, ~, ~, ~, ~] = compute_phi(E_plot, c_tmp);
    
    H_calc_F(i) = trapz(E_plot, Phi_tmp_tot(:)' .* h74F_interp);
    H_est_F(i)  = trapz(E_plot, Phi_tmp_tot(:)' .* Respuesta_CL_F);
end

% Ajuste por mínimos cuadrados para encontrar el factor de corrección C_F
C_F = H_est_F \ H_calc_F; 
H_est_F_corr = H_est_F * C_F;

% Cálculo de errores relativos (%)
err_F_orig = 100 * (H_est_F - H_calc_F) ./ H_calc_F;
err_F_corr = 100 * (H_est_F_corr - H_calc_F) ./ H_calc_F;

% Impresión de resultados en consola
fprintf('====================================================================\n');
fprintf('         RESULTADOS DE CALIBRACIÓN - MODELO ICRP 74F (FERRARI)\n');
fprintf('====================================================================\n');
fprintf(' * Factor de Corrección Global (C_F) = %.4f\n', C_F);
fprintf(' * Ecuación de calibración: H*_corregida = %.4f * H*_estimada_bruta\n', C_F);
fprintf(' * Error Absoluto Medio Inicial: %.2f%%\n', mean(abs(err_F_orig)));
fprintf(' * Error Absoluto Medio Post-Corrección: %.2f%%\n', mean(abs(err_F_corr)));
fprintf('====================================================================\n');

% Gráficas exclusivas para Ferrari - Distribución optimizada 2x2
figure('Name', 'Calibración Exclusiva - ICRP 74F (Ferrari)', 'Color', 'w', 'Position', [100 100 1200 700]);

% --- SUBPLOT 1 (Arriba a la Izquierda): Evolución del Error Relativo ---
subplot(2,2,1);
plot(1:num_casos, err_F_orig, 'r--o', 'LineWidth', 1.2, 'DisplayName', 'Error Bruto'); hold on;
plot(1:num_casos, err_F_corr, 'g-s', 'LineWidth', 1.6, 'MarkerFaceColor', 'g', 'DisplayName', 'Error Corregido');
yline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Línea ideal');
xlabel('Caso (ID)', 'FontWeight', 'bold');
ylabel('Error Relativo (%)', 'FontWeight', 'bold');
title('Evolución del Error Relativo - Ferrari');
legend('Location', 'best'); grid on; 
xticks(1:num_casos); xlim([0.5 num_casos+0.5]);

% --- SUBPLOT 3 (Abajo a la Izquierda): NUEVO - Comparativa de Dosis Absolutas por ID ---
subplot(2,2,3);
plot(1:num_casos, H_calc_F, 'k-d', 'LineWidth', 1.8, 'MarkerFaceColor', 'k', 'DisplayName', 'H* Calculada'); hold on;
plot(1:num_casos, H_est_F, 'r--o', 'LineWidth', 1.2, 'DisplayName', 'H* Estimada');
% plot(1:num_casos, H_est_F_corr, 'g-.s', 'LineWidth', 1.4, 'MarkerFaceColor', 'g', 'DisplayName', 'H* Estimada');
xlabel('Caso (ID)', 'FontWeight', 'bold');
ylabel('H*(10) [pSv/s]', 'FontWeight', 'bold');
title('Normativa ICRP 74 con la extensión de Ferrari');
legend('Location', 'best'); grid on;
xticks(1:num_casos); xlim([0.5 num_casos+0.5]);

% --- SUBPLOT 2 y 4 (Columna Derecha Completa): Correlación de Dosis ---
subplot(2,2,[2,4]);
scatter(H_est_F, H_calc_F, 50, 'r', 'filled', 'DisplayName', 'Datos Brutos'); hold on;
scatter(H_est_F_corr, H_calc_F, 50, 'g', 'filled', 'DisplayName', 'Datos Corregidos');
x_fit_F = linspace(min(H_est_F)*0.9, max([H_est_F; H_est_F_corr])*1.1, 100);
plot(x_fit_F, C_F*x_fit_F, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Ajuste Lineal (C_F = %.3f)', C_F));
plot([0 max([H_est_F; H_est_F_corr])], [0 max([H_est_F; H_est_F_corr])], 'k:', 'LineWidth', 1.2, 'DisplayName', 'Línea Ideal (1:1)');
xlabel('H*(10) Estimada (pSv/s)', 'FontWeight', 'bold');
ylabel('H*(10) Calculada (pSv/s)', 'FontWeight', 'bold');
title('Correlación Dosis Estimada vs Calculada - Ferrari');
legend('Location', 'northwest'); grid on; axis tight;

%% 4B. ICRP 74S (SANNIKOV)
fprintf('\n\n------ INICIO DEL ANÁLISIS AISLADO: ICRP 74S (SANNIKOV) ------ \n\n');

% Inicializar vectores de dosis para Sannikov
H_calc_S = zeros(num_casos, 1); 
H_est_S  = zeros(num_casos, 1);

% Combinación lineal original para Sannikov
Respuesta_CL_S = datos.vRbare_interp*0 + datos.vRcore_interp*9.239092852385420 + datos.vRhepb_interp*8.603464735268430;
h74S_interp = datos.vh1074sannikov_interp(:)'; 

% Bucle de cálculo para Sannikov
for i = 1:num_casos
    inc = casos_variacion(i, :) / 100;
    factor_crecimiento = 1 + sum(p_viejos .* inc);
    c_tmp = coef; 
    c_tmp.p1 = (coef.p1 * (1 + inc(1))) / factor_crecimiento;
    c_tmp.p2 = (coef.p2 * (1 + inc(2))) / factor_crecimiento;
    c_tmp.p3 = (coef.p3 * (1 + inc(3))) / factor_crecimiento;
    c_tmp.p4 = (coef.p4 * (1 + inc(4))) / factor_crecimiento;
    c_tmp.k  = coef.k * factor_crecimiento;
    
    [Phi_tmp_tot, ~, ~, ~, ~] = compute_phi(E_plot, c_tmp);
    
    H_calc_S(i) = trapz(E_plot, Phi_tmp_tot(:)' .* h74S_interp);
    H_est_S(i)  = trapz(E_plot, Phi_tmp_tot(:)' .* Respuesta_CL_S);
end

% Ajuste por mínimos cuadrados para encontrar el factor de corrección C_S
C_S = H_est_S \ H_calc_S; 
H_est_S_corr = H_est_S * C_S;

% Cálculo de errores relativos (%)
err_S_orig = 100 * (H_est_S - H_calc_S) ./ H_calc_S;
err_S_corr = 100 * (H_est_S_corr - H_calc_S) ./ H_calc_S;

% Impresión de resultados en consola
fprintf('====================================================================\n');
fprintf('         RESULTADOS DE CALIBRACIÓN - MODELO ICRP 74S (SANNIKOV)\n');
fprintf('====================================================================\n');
fprintf(' * Factor de Corrección Global (C_S) = %.4f\n', C_S);
fprintf(' * Ecuación de calibración: H*_corregida = %.4f * H*_estimada_bruta\n', C_S);
fprintf(' * Error Absoluto Medio Inicial: %.2f%%\n', mean(abs(err_S_orig)));
fprintf(' * Error Absoluto Medio Post-Corrección: %.2f%%\n', mean(abs(err_S_corr)));
fprintf('====================================================================\n');

% Gráficas exclusivas para Sannikov - Distribución optimizada 2x2
figure('Name', 'Calibración Exclusiva - ICRP 74S (Sannikov)', 'Color', 'w', 'Position', [150 150 1200 700]);

% --- SUBPLOT 1 (Arriba a la Izquierda): Evolución del Error Relativo ---
subplot(2,2,1);
plot(1:num_casos, err_S_orig, 'b--o', 'LineWidth', 1.2, 'DisplayName', 'Error Bruto'); hold on;
plot(1:num_casos, err_S_corr, 'g-s', 'LineWidth', 1.6, 'MarkerFaceColor', 'g', 'DisplayName', 'Error Corregido');
yline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Línea ideal');
xlabel('Caso (ID)', 'FontWeight', 'bold');
ylabel('Error Relativo (%)', 'FontWeight', 'bold');
title('Evolución del Error Relativo - Sannikov');
legend('Location', 'best'); grid on; 
xticks(1:num_casos); xlim([0.5 num_casos+0.5]);

% --- SUBPLOT 3 (Abajo a la Izquierda): Comparativa de Dosis Absolutas por ID ---
subplot(2,2,3);
plot(1:num_casos, H_calc_S, 'k-d', 'LineWidth', 1.8, 'MarkerFaceColor', 'k', 'DisplayName', 'H* Calculada'); hold on;
plot(1:num_casos, H_est_S, 'b--o', 'LineWidth', 1.2, 'DisplayName', 'H* Estimada');
% plot(1:num_casos, H_est_S_corr, 'g-.s', 'LineWidth', 1.4, 'MarkerFaceColor', 'g', 'DisplayName', 'H* Estimada Corregida');
xlabel('Caso (ID)', 'FontWeight', 'bold');
ylabel('H*(10) [pSv/s]', 'FontWeight', 'bold');
title('Normativa ICRP 74 con la extensión de Sannikov');
legend('Location', 'best'); grid on;
xticks(1:num_casos); xlim([0.5 num_casos+0.5]);

% --- SUBPLOT 2 y 4 (Columna Derecha Completa): Correlación de Dosis ---
subplot(2,2,[2,4]);
scatter(H_est_S, H_calc_S, 50, 'b', 'filled', 'DisplayName', 'Datos Brutos'); hold on;
scatter(H_est_S_corr, H_calc_S, 50, 'g', 'filled', 'DisplayName', 'Datos Corregidos');
x_fit_S = linspace(min(H_est_S)*0.9, max([H_est_S; H_est_S_corr])*1.1, 100);
plot(x_fit_S, C_S*x_fit_S, 'm-', 'LineWidth', 2, 'DisplayName', sprintf('Ajuste Lineal (C_S = %.3f)', C_S));
plot([0 max([H_est_S; H_est_S_corr])], [0 max([H_est_S; H_est_S_corr])], 'k:', 'LineWidth', 1.2, 'DisplayName', 'Línea Ideal (1:1)');
xlabel('H*(10) Estimada (pSv/s)', 'FontWeight', 'bold');
ylabel('H*(10) Calculada (pSv/s)', 'FontWeight', 'bold');
title('Correlación Dosis Estimada vs Calculada - Sannikov');
legend('Location', 'northwest'); grid on; axis tight;

%% 4C. ICRU 95
fprintf('\n\n------ INICIO DEL ANÁLISIS AISLADO: ICRU 95 ------ \n\n');

% Inicializar vectores de dosis para ICRU 95
H_calc_95 = zeros(num_casos, 1); 
H_est_95  = zeros(num_casos, 1);

% =========================================================================
% REVISAR: Sustituye estos coeficientes por los valores reales que obtuviste
% en tu optimización específica para la curva de la ICRU 95.
% =========================================================================
Respuesta_CL_95 = datos.vRbare_interp*0 + datos.vRcore_interp*0.915052439549803 + datos.vRhepb_interp*14.269786202324562; 

h95_interp = datos.vh95_interp(:)'; 

% Bucle de cálculo para ICRU 95
for i = 1:num_casos
    inc = casos_variacion(i, :) / 100;
    factor_crecimiento = 1 + sum(p_viejos .* inc);
    c_tmp = coef; 
    c_tmp.p1 = (coef.p1 * (1 + inc(1))) / factor_crecimiento;
    c_tmp.p2 = (coef.p2 * (1 + inc(2))) / factor_crecimiento;
    c_tmp.p3 = (coef.p3 * (1 + inc(3))) / factor_crecimiento;
    c_tmp.p4 = (coef.p4 * (1 + inc(4))) / factor_crecimiento;
    c_tmp.k  = coef.k * factor_crecimiento;
    
    [Phi_tmp_tot, ~, ~, ~, ~] = compute_phi(E_plot, c_tmp);
    
    H_calc_95(i) = trapz(E_plot, Phi_tmp_tot(:)' .* h95_interp);
    H_est_95(i)  = trapz(E_plot, Phi_tmp_tot(:)' .* Respuesta_CL_95);
end

% Ajuste por mínimos cuadrados para encontrar el factor de corrección C_95
C_95 = H_est_95 \ H_calc_95; 
H_est_95_corr = H_est_95 * C_95;

% Cálculo de errores relativos (%)
err_95_orig = 100 * (H_est_95 - H_calc_95) ./ H_calc_95;
err_95_corr = 100 * (H_est_95_corr - H_calc_95) ./ H_calc_95;

% Impresión de resultados en consola
fprintf('====================================================================\n');
fprintf('         RESULTADOS DE CALIBRACIÓN - MODELO ICRU 95\n');
fprintf('====================================================================\n');
fprintf(' * Factor de Corrección Global (C_95) = %.4f\n', C_95);
fprintf(' * Ecuación de calibración: H*_corregida = %.4f * H*_estimada_bruta\n', C_95);
fprintf(' * Error Absoluto Medio Inicial: %.2f%%\n', mean(abs(err_95_orig)));
fprintf(' * Error Absoluto Medio Post-Corrección: %.2f%%\n', mean(abs(err_95_corr)));
fprintf('====================================================================\n');

% Definición del color Naranja personalizado [R G B]
color_naranja = [1.0 0.5 0.0];

% Gráficas exclusivas para ICRU 95 - Distribución optimizada 2x2
figure('Name', 'Calibración Exclusiva - ICRU 95', 'Color', 'w', 'Position', [200 200 1200 700]);

% --- SUBPLOT 1 (Arriba a la Izquierda): Evolución del Error Relativo ---
subplot(2,2,1);
plot(1:num_casos, err_95_orig, '--o', 'Color', color_naranja, 'LineWidth', 1.2, 'DisplayName', 'Error Bruto'); hold on;
plot(1:num_casos, err_95_corr, 'g-s', 'LineWidth', 1.6, 'MarkerFaceColor', 'g', 'DisplayName', 'Error Corregido');
yline(0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Línea ideal');
xlabel('Caso (ID)', 'FontWeight', 'bold');
ylabel('Error Relativo (%)', 'FontWeight', 'bold');
title('Evolución del Error Relativo - ICRU 95');
legend('Location', 'best'); grid on; 
xticks(1:num_casos); xlim([0.5 num_casos+0.5]);

% --- SUBPLOT 3 (Abajo a la Izquierda): Comparativa de Dosis Absolutas por ID ---
subplot(2,2,3);
plot(1:num_casos, H_calc_95, 'k-d', 'LineWidth', 1.8, 'MarkerFaceColor', 'k', 'DisplayName', 'H* Calculada'); hold on;
plot(1:num_casos, H_est_95, '--o', 'Color', color_naranja, 'LineWidth', 1.2, 'DisplayName', 'H* Estimada');
% plot(1:num_casos, H_est_95_corr, 'g-.s', 'LineWidth', 1.4, 'MarkerFaceColor', 'g', 'DisplayName', 'H* Estimada Corregida');
xlabel('Caso (ID)', 'FontWeight', 'bold');
ylabel('H*(10) [pSv/s]', 'FontWeight', 'bold');
title('Normativa ICRU 95');
legend('Location', 'best'); grid on;
xticks(1:num_casos); xlim([0.5 num_casos+0.5]);

% --- SUBPLOT 2 y 4 (Columna Derecha Completa): Correlación de Dosis ---
subplot(2,2,[2,4]);
scatter(H_est_95, H_calc_95, 50, color_naranja, 'filled', 'DisplayName', 'Datos Brutos'); hold on;
scatter(H_est_95_corr, H_calc_95, 50, 'g', 'filled', 'DisplayName', 'Datos Corregidos');
x_fit_95 = linspace(min(H_est_95)*0.9, max([H_est_95; H_est_95_corr])*1.1, 100);
plot(x_fit_95, C_95*x_fit_95, 'm-', 'LineWidth', 2, 'DisplayName', sprintf('Ajuste Lineal (C_{95} = %.3f)', C_95));
plot([0 max([H_est_95; H_est_95_corr])], [0 max([H_est_95; H_est_95_corr])], 'k:', 'LineWidth', 1.2, 'DisplayName', 'Línea Ideal (1:1)');
xlabel('H*(10) Estimada (pSv/s)', 'FontWeight', 'bold');
ylabel('H*(10) Calculada (pSv/s)', 'FontWeight', 'bold');
title('Correlación Dosis Estimada vs Calculada - ICRU 95');
legend('Location', 'northwest'); grid on; axis tight;

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================

function [Phi_tot, phi_t, phi_e, phi_f, phi_h] = compute_phi(E, P)
    
    E = E(:)';

    % 1. Definición de funciones base
    u_t = (E ./ P.To^2) .* exp(-E ./ P.To);
    u_e = (1 - exp(-(E./P.Ed).^2)) .* max(E,1e-40).^(P.b-1) .* exp(-E./P.bp);
    u_f    = exp( -P.bf .* (log(E)).^2 + P.gf .* log(E) );
    u_h    = exp( -P.bh .* (log(E)).^2 + P.gh .* log(E) );
    
    % Elimina NaNs e Infinitos sustituyéndolos por cero
    u_t(~isfinite(u_t)) = 0; u_e(~isfinite(u_e)) = 0;
    u_f(~isfinite(u_f)) = 0; u_h(~isfinite(u_h)) = 0;

    % Normalización inicial de las formas sobre todo el dominio para estabilidad
    phi_t    = u_t    ./ max(trapz(E, u_t),    1e-35);
    phi_e    = u_e    ./ max(trapz(E, u_e),    1e-35);
    phi_f    = u_f    ./ max(trapz(E, u_f),    1e-35);
    phi_h    = u_h    ./ max(trapz(E, u_h),    1e-35);


    % 2. Definir las máscaras para los 4 rangos físicos (E está en MeV)
    m1 = (E <= 3.5e-7);             % R1: [Emin, 0.35 eV]  -> 0.35 eV = 3.5e-7 MeV
    m2 = (E > 3.5e-7) & (E <= 0.1); % R2: (0.35 eV, 0.1 MeV]
    m3 = (E > 0.1)    & (E <= 20);  % R3: (0.1 MeV, 20 MeV]
    m4 = (E > 20);                  % R4: (20 MeV, Emax]

    masks  = {m1, m2, m3, m4};
    shapes = {phi_t, phi_e, phi_f, phi_h};

    % 3. Construir Matriz M (4x4)
    % M(i,j) = Integral en la Región i de la Función de Forma j
    M = zeros(4,4);
    for i = 1:4
        idx = masks{i};
        if sum(idx) > 1
            for j = 1:4
                M(i,j) = trapz(E(idx), shapes{j}(idx));
            end
        end
    end

    % 4. Vector de integrales objetivo en las 4 zonas (sliders del usuario)
    p_target = max([P.p1, P.p2, P.p3, P.p4]', 1e-10);
    p_target = p_target / sum(p_target); % Imponer que la suma de áreas sea 1

    % 5. Resolver el sistema M * w = p_target para encontrar las amplitudes
    % Usamos lsqnonneg para evitar amplitudes negativas debido al solapamiento
    opts = optimset('Display','off');
    w = lsqnonneg(M, p_target, opts);

    % 6. Construir espectro normalizado y forzar integral total = 1 numéricamente
    Phi_norm = w(1)*phi_t + w(2)*phi_e + w(3)*phi_f + w(4)*phi_h;
    Phi_norm = Phi_norm ./ max(trapz(E, Phi_norm), 1e-35);

    % 7. Escalar por el factor global k
    Phi_tot = P.k * Phi_norm;

    % Componentes individuales escaladas para el plot
    phi_t = P.k * w(1) * phi_t;
    phi_e = P.k * w(2) * phi_e;
    phi_f = P.k * w(3) * phi_f;
    phi_h = P.k * w(4) * phi_h;

    % --- CORTE DE ALTA ENERGÍA (100000 MeV) ---
    E_corte = 100000; 
    Phi_tot(E > E_corte) = 0;
    phi_t(E > E_corte) = 0;
    phi_e(E > E_corte) = 0;
    phi_f(E > E_corte) = 0;
    phi_h(E > E_corte) = 0;
end