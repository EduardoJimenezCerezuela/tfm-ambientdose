% =========================================================================
%  coeficientes.m 
%  Autor: Eduardo Jiménez Cerezuela
% =========================================================================
%  HENSA++ – Ajuste de coeficientes según IEC 61005:2014
%  Representación por BINES de energía (Sección 6.4.3 + Tabla 2)
%  Optimización: mínimos cuadrados ponderados (WLS) sobre respuesta relativa
% =========================================================================

addpath('C:\Users\eduji\Desktop\Máster\TFM\Data')

close all; clear all; clc;

%% 1. LECTURA DE DATOS
[vEphi, vPhi]            = read_2cols('EXPACSJavalambre_phi.txt');

% Cambiar según el modo
modo = 'ICRU95';   % Opciones: 'Ferrari' | 'Sannikov' | 'ICRU95'

switch modo
    case 'ICRU95'
        [vEh, vh] = read_2cols('h_Emax_ICRU95');

    case 'Ferrari'
        [vEh, vh] = read_2cols('h10_ICRP74');
        [vEh2, vh2] = read_2cols('h10_extension_ferrarietal1997'); % tu fichero Ferrari
        % ajusta el umbral de empalme según tus datos
        vh2  = vh2(vEh2 >= 500, :);
        vEh2 = vEh2(vEh2 >= 500, :);
        vEh = [vEh; vEh2];
        vh  = [vh;  vh2];

    case 'Sannikov'
        [vEh, vh] = read_2cols('h10_ICRP74');
        [vEh2, vh2] = read_2cols('h10_extension_sannikovetal1997');
        vh2  = vh2(vEh2 >= 300, :);
        vEh2 = vEh2(vEh2 >= 300, :);
        vEh = [vEh; vEh2];
        vh  = [vh;  vh2];
end

[vER, vRbare, vRcore, vRhepb] = read_4cols('Response_Functions_HENSApp_OAJ_Phase0.dat');

fprintf('[LECTURA] EXPACS: %d pts | h*(10): %d pts | Respuesta: %d pts\n', ...
        length(vEphi), length(vEh), length(vER));


%%  2. BINES DE ENERGÍA IEC 61005:2014
%
%  La sección 6.4.3 exige los siguientes puntos de validación:
%    a) ≥2 energías < 50 keV,  una de ellas térmica
%    b) ≥3 energías en 50 keV – 10 MeV
%    c) ≥1 fuente extensa (252Cf / 241Am-Be)
%    d) ≥1 energía > 10 MeV
%
%  La sección 6.5.2 exige al menos un punto calculado por DÉCADA de energía.
%  Los bines siguientes cubren ambas exigencias.
%
%  Tolerancias normativas (Tabla 2 / Sección 6.4.2):
%    Zona térmica–50 keV   → respuesta relativa en [0.2, 8.0]  
%    Zona 50 keV–10 MeV   → respuesta relativa en [0.5, 2.0]  
%    Zona >10 MeV          → respuesta relativa en [0.2, 2.0]  
% =========================================================================

%                    E (MeV)     Zona
E_bins = [
    2.53e-8;   % Térmico (0.025 eV)    — zona a)  [0.2–8.0]
    1.00e-7;   % 0.1 eV                — zona a)
    1.00e-6;   % 1 eV                  — zona a)
    5.00e-6;   % 5 eV                  — zona a)
    1.00e-5;   % 10 eV                 — zona a)
    5.00e-5;   % 50 eV                 — zona a)
    1.00e-4;   % 100 eV                — zona a)
    5.00e-4;   % 500 eV                — zona a)
    1.00e-3;   % 1 keV                 — zona a)
    1.00e-2;   % 10 keV                — zona a)
    5.00e-2;   % 50 keV  ← frontera   — zona a)/b)
    1.00e-1;   % 100 keV               — zona b)  [0.5–2.0]
    2.50e-1;   % 250 keV               — zona b)
    5.00e-1;   % 500 keV               — zona b)
    1.00e+0;   % 1 MeV                 — zona b)
    2.00e+0;   % 2 MeV                 — zona b)
    5.00e+0;   % 5 MeV                 — zona b)
    1.00e+1;   % 10 MeV ← frontera   — zona b)/c)
    1.40e+1;   % 14 MeV (fuente D-T)  — zona c)  [0.2–2.0]
    2.00e+1;   % 20 MeV               — zona c)
    5.00e+1;   % 50 MeV               — zona c)
    7.50e+1;   % 75 MeV               — zona c)
    1.00e+2;   % 100 MeV               — zona c)
    1.50e+2;   % 150 MeV               — zona c)
    2.00e+2;   % 200 MeV               — zona c)
    5.00e+2;   % 500 MeV               — zona c)
    7.50e+2;   % 750 MeV               — zona c)
    1.00e+3;   % 1 GeV               — zona c)
];

% Límites normativos en cada bin [lo, hi]
lim_lo = zeros(length(E_bins), 1);
lim_hi = zeros(length(E_bins), 1);
for k = 1:length(E_bins)
    if E_bins(k) < 0.050
        lim_lo(k) = 0.2;  lim_hi(k) = 8;  % Zona térmica–50 keV
    elseif E_bins(k) <= 10.0
        lim_lo(k) = 0.5;  lim_hi(k) = 2;  % Zona crítica (OBLIGATORIA)
    else
        lim_lo(k) = 0.2;  lim_hi(k) = 2;  % Zona >10 MeV
    end
end

% Sigma para chi² ponderado:  σ_k = (hi_k - lo_k) / 2  (semi-amplitud del intervalo)
% Esto hace que chi²=1 cuando r se aleja exactamente una semi-amplitud de r=1.
sigma_rel = (lim_hi - lim_lo) / 2;

% =========================================================================
%  PONDERACIÓN sigma_rel  (depende del modo)
% =========================================================================
switch modo

    case 'ICRU95'
        idx_evap = (E_bins >= 0.5 & E_bins <= 5.0); % Pico de Evaporación (0.5 a 5 MeV)
        sigma_rel(idx_evap) = sigma_rel(idx_evap) / 7;

        idx_transicion = (E_bins > 10 & E_bins < 50); % Zona de cruce (10 a 40 MeV)
        sigma_rel(idx_transicion) = sigma_rel(idx_transicion) / 3;
        
        idx_cascada = (E_bins >= 50 & E_bins <= 250); % Pico de Cascada (50 a 200 MeV)
        sigma_rel(idx_cascada) = sigma_rel(idx_cascada) / 6;

    case 'Ferrari'
        % Diferencias estructurales clave vs ICRU95:
        %  - Pico de evaporación desplazado a ~1-2 MeV (416-425 pSv·cm²)
        %  - PICO SECUNDARIO pronunciado en 14-20 MeV (520-600 pSv·cm²)
        %    vs ICRU95 que DESCIENDE en esa zona → comportamiento opuesto
        %  - Descenso fuerte 20→150 MeV (600→245 pSv·cm²) mucho más abrupto
        %  - Mínimo extendido y profundo 200-500 MeV (Ferrari: 221→290 pSv·cm²,
        %    frente a ICRU95 que sube suavemente 448→533)
        %  - Ferrari tiene pocos puntos a alta energía → subida lenta y esparcida

        idx_evap = (E_bins >= 0.5 & E_bins <= 5.0); % Pico de Evaporación (0.5 a 5 MeV)
        sigma_rel(idx_evap) = sigma_rel(idx_evap) / 5;

        idx_pico_sec = (E_bins > 10 & E_bins <= 30); % Pico secundario ICRP74 (10 a 30 MeV) - Exige más precisión que la zona de cruce de ICRU95 (que es monótona)
        sigma_rel(idx_pico_sec) = sigma_rel(idx_pico_sec) / 8;

        idx_min_ferr = (E_bins > 30 & E_bins <= 300); % Mínimo extendido Ferrari (30 a 300 MeV) - Cubre desde la caída post-pico (30→50 MeV) hasta el fondo Ferrari (~200-500 MeV)
        sigma_rel(idx_min_ferr) = sigma_rel(idx_min_ferr) / 4;

    case 'Sannikov'
        % Diferencias estructurales clave vs ICRU95:
        %  - Pico de evaporación desplazado a ~1-2 MeV (416-425 pSv·cm²)
        %  - PICO SECUNDARIO pronunciado en 14-20 MeV (520-600 pSv·cm²)
        %    vs ICRU95 que DESCIENDE en esa zona → comportamiento opuesto
        %  - Descenso fuerte 20→150 MeV (600→245 pSv·cm²), mínimo en 150 MeV
        %  - Empalme suave Sannikov en 300 MeV (306 pSv·cm² ≈ ICRP74 en 260)
        %  - Subida rápida Sannikov 300→750 MeV (306→558 pSv·cm²)
        %    llega a ~647 a 1 GeV, casi idéntico a ICRU95 en ese punto

        idx_evap_crit = (E_bins >= 0.5 & E_bins <= 10.0); % Zona IEC crítica (0.5 a 10 MeV) - Se amplía hasta 10 MeV (antes solo 0.5-5) para anclar el pico secundario
        sigma_rel(idx_evap_crit) = sigma_rel(idx_evap_crit) / 5;

        idx_pico_sec = (E_bins > 10 & E_bins <= 30);
        sigma_rel(idx_pico_sec) = sigma_rel(idx_pico_sec) / 6; % 5 Pico secundario ICRP74 (10 a 30 MeV) — Peso moderado: ya está anclado por el bin de 10 MeV del paso anterior

        idx_cascada = (E_bins > 30 & E_bins <= 200); % Descenso cascada ICRP74 (30 a 200 MeV) - Peso reducido: no debe competir con la zona crítica
        sigma_rel(idx_cascada) = sigma_rel(idx_cascada) / 3;

        idx_sannikov = (E_bins > 300 & E_bins <= 750); % Subida Sannikov (300 MeV a 750 MeV)
        sigma_rel(idx_sannikov) = sigma_rel(idx_sannikov) / 2;
end


%% 3. FILTRADO Y RESTRICCIÓN AL RANGO VÁLIDO COMÚN

Emin_valid = max([vEphi(1), vEh(1), vER(1)]);
Emax_valid = min([vEphi(end), vEh(end), vER(end)]);

switch modo
    case 'ICRU95'
        mask      = (E_bins >= Emin_valid) & (E_bins <= Emax_valid);

    case 'Ferrari'
        mask      = (E_bins >= Emin_valid) & (E_bins <= Emax_valid);
        mask = mask & (E_bins >= 1.0); % Priorizamos el ajuste a partir del MeV

    case 'Sannikov'
        mask      = (E_bins >= Emin_valid) & (E_bins <= Emax_valid);
        mask = mask & (E_bins >= 1.0); % Priorizamos el ajuste a partir del MeV
end

E_bins    = E_bins(mask);
sigma_rel = sigma_rel(mask);
lim_lo    = lim_lo(mask);
lim_hi    = lim_hi(mask);
nBins     = sum(mask);

fprintf('[BINES] %d bines activos en [%.2e, %.2e] MeV\n\n', ...
        nBins, E_bins(1), E_bins(end));

%%  4. INTERPOLACIÓN LOG-LOG EN LOS BINES
%     Solo se evalúan los puntos necesarios, sin rejilla densa.
%
%   Ev: Son los puntos de la muestra (x)
%   Fv: valores correpondientes a Ev (y)
%   Eq: Los puntos que queremos consultar (Los bines)
%
% =========================================================================
interp_ll = @(Ev, Fv, Eq) ...
    10.^interp1(log10(Ev), log10(Fv), log10(Eq), 'pchip'); % Represemtación log-log

h_bins  = interp_ll(vEh,   vh,     E_bins);
Rb_bins = interp_ll(vER,   vRbare, E_bins);
Rc_bins = interp_ll(vER,   vRcore, E_bins);
Rh_bins = interp_ll(vER,   vRhepb, E_bins);
% (phi_bins no entra en la optimización de coeficientes, solo en la dosis final)

%%  5. FORMULACIÓN WLS – RESPUESTA RELATIVA
%
%  Definimos la respuesta relativa en el bin k:
%
%    r_k(x) = [x1·Rb_k + x2·Rc_k + x3·Rh_k] / h_k
%
%  El chi² normalizado (IEC) es:
%
%    chi²(x) = Σ_k  [(r_k - 1) / σ_k]²
%
%  Es equivalente al WLS lineal:
%
%    A_norm · x ≈ b_norm
%
%  donde  A_norm_kj = R_j(k) / [h_k · σ_k]
%         b_norm_k  = 1 / σ_k
% =========================================================================

A_norm = [(Rb_bins ./ h_bins) ./ sigma_rel, ...
          (Rc_bins ./ h_bins) ./ sigma_rel, ...
          (Rh_bins ./ h_bins) ./ sigma_rel];
b_norm = ones(nBins, 1) ./ sigma_rel;

% --- Solución: WLS con positividad de coeficientes (lsqnonneg) ---
%  Si los coeficientes deben ser positivos (por sentido físico),
%  usamos lsqnonneg; de lo contrario, x_wls ya es óptimo.
x_pos = lsqnonneg(A_norm, b_norm);


%%  6. EVALUACIÓN: CHI², RESPUESTA RELATIVA Y VALIDACIÓN IEC

x_sol  = x_pos;
tag    = 'WLS (x≥0)';

r_bins = (Rb_bins*x_sol(1) + Rc_bins*x_sol(2) + Rh_bins*x_sol(3)) ./ h_bins;
chi2_k = ((r_bins - 1) ./ sigma_rel).^2;
chi2   = sum(chi2_k);
nOK    = sum((r_bins >= lim_lo) & (r_bins <= lim_hi));

fprintf('=== %s ===\n', tag);
fprintf('  x = [%.5f,  %.5f,  %.5f]\n', x_sol(1), x_sol(2), x_sol(3));
fprintf('  Chi² = %.4f  (referencia ideal = nBins = %d)\n', chi2, nBins);
fprintf('  Bines IEC OK: %d / %d\n\n', nOK, nBins);

r_pos = (Rb_bins*x_pos(1) + Rc_bins*x_pos(2) + Rh_bins*x_pos(3)) ./ h_bins;
nOK_pos = sum((r_pos >= lim_lo) & (r_pos <= lim_hi));

x_best = x_pos;  r_best = r_pos;  label_best = 'WLS (x≥0)';


% Residuos
res = A_norm*x_pos - b_norm;

% Grados de libertad
nu = nBins - length(x_pos);

% Varianza residual
s2 = sum(res.^2)/nu;

% Matriz de covarianza
Cov = s2 * inv(A_norm' * A_norm);

% Incertidumbres 1 sigma
sigma_x = sqrt(diag(Cov));

fprintf('\nCoeficientes WLS:\n');
for i=1:length(x_pos)
    fprintf('x%d = %.6f ± %.6f\n', ...
            i, x_pos(i), sigma_x(i));
end

% Correlaciones entre coeficientes
Corr = Cov ./ sqrt(diag(Cov)*diag(Cov)');
disp(Corr)


%%  7. TABLA RESUMEN POR BIN

fprintf('%-14s  %-8s  %-10s  %-10s  %-8s  %s\n', ...
        'E (MeV)', 'r_opt', '[lo, hi]', 'chi2_bin', 'sigma', 'Estado');
fprintf('%s\n', repmat('-', 1, 70));
for k = 1:nBins
    r  = r_best(k);
    c2 = ((r-1)/sigma_rel(k))^2;
    ok = (r >= lim_lo(k)) && (r <= lim_hi(k));
    if ok,  estado = 'OK'; else, estado = 'FUERA'; end
    fprintf('%.3e MeV   r=%5.3f  [%.1f, %.1f]  chi2=%6.3f  s=%.2f  %s\n', ...
            E_bins(k), r, lim_lo(k), lim_hi(k), c2, sigma_rel(k), estado);
end
fprintf('%s\n', repmat('-', 1, 70));
fprintf('Chi² TOTAL: %.4f   |   Solución elegida: %s\n', ...
        sum(((r_best-1)./sigma_rel).^2), label_best);

nu = nBins - 3; % Grados de libertad (los 3 coeficientes que buscamos)
chi2_red = sum(((r_best-1)./sigma_rel).^2) / nu;
fprintf('Chi² reducida = %.4f\n\n', chi2_red);


%%  8. REJILLA DENSA PARA GRÁFICAS (solo aquí, no para optimización)

nGrid  = 5000;
E_grid = logspace(log10(Emin_valid), log10(Emax_valid), nGrid)';

h_g  = interp_ll(vEh,   vh,     E_grid);
Rb_g = interp_ll(vER,   vRbare, E_grid);
Rc_g = interp_ll(vER,   vRcore, E_grid);
Rh_g = interp_ll(vER,   vRhepb, E_grid);
phi_g = interp_ll(vEphi, vPhi,  E_grid);

r_grid = (Rb_g*x_best(1) + Rc_g*x_best(2) + Rh_g*x_best(3)) ./ h_g;


%%  9. CÁLCULO DE DOSIS AMBIENTAL INTEGRADA Y TASA DE CONTEO

H_calc = trapz(E_grid, (Rb_g*x_best(1) + Rc_g*x_best(2) + Rh_g*x_best(3)) .* phi_g);
H_ref  = trapz(E_grid, h_g .* phi_g);

fprintf('[DOSIS INTEGRADA]\n');
fprintf('  H*(10) calculado:  %.4e  (pSv/s)\n', H_calc);
fprintf('  H*(10) referencia: %.4e  (pSv/s)\n', H_ref);
fprintf('  Error relativo:  %.4f%% (ideal = 0%%)\n\n', abs((((H_calc - H_ref) / H_ref))*100));

CPS_bare = trapz(E_grid, Rb_g .* phi_g);
CPS_core = trapz(E_grid, Rc_g .* phi_g);
CPS_hepb = trapz(E_grid, Rh_g .* phi_g);

CPS_total = CPS_bare + CPS_core + CPS_hepb;

fprintf('[TASA DE CONTEO]\n');
fprintf('  Bare:   %.4e (cps)\n', CPS_bare);
fprintf('  Core:   %.4e (cps)\n', CPS_core);
fprintf('  HEPB:   %.4e (cps)\n', CPS_hepb);
fprintf('%s\n', repmat('-', 1, 40));
fprintf('  Tc calculado:  %.4e  (cps)\n', CPS_total);


%% 10. GRÁFICA DE VALIDACIÓN NORMATIVA

figure('Name', 'Validación IEC 61005:2014 – Respuesta relativa', ...
       'Color', 'w', 'Position', [80 80 950 520]);

% Zonas de color (fondo) según Tabla 2
E_keV = E_grid * 1e3;  % Convertir a keV para el eje X
patch([1e-5 50 50 1e-5 1e-5], [0.2 0.2 8.0 8.0 0.2], ...
      [1.0 0.88 0.70], 'EdgeColor','none','FaceAlpha',0.4); hold on;
patch([50 1e4 1e4 50 50], [0.5 0.5 2.0 2.0 0.5], ...
      [0.70 0.90 0.70], 'EdgeColor','none','FaceAlpha',0.4);
patch([1e4 2.5e6 2.5e6 1e4 1e4], [0.2 0.2 2.0 2.0 0.2], ...
      [0.70 0.75 1.0], 'EdgeColor','none','FaceAlpha',0.4);

set(gca, 'FontSize', 13, 'FontWeight', 'bold');

% Respuesta continua
semilogx(E_keV, r_grid, 'b-', 'LineWidth', 2.0);

% Bines IEC
semilogx(E_bins*1e3, r_best, 'ko', 'MarkerSize', 9, ...
         'MarkerFaceColor', 'w', 'LineWidth', 1.8);

% Línea de referencia ideal
yline(1.0, 'k--', 'LineWidth', 1.2);
text(1e-4, 1.2, 'r = 1 (ideal)', 'FontSize', 13, 'Color', 'k');

% Límites críticos 50 keV–10 MeV
yline(2.0, 'r--', 'LineWidth', 1.0);
text(1e-2, 2.3, 'Límite sup. zona crítica (×2)', 'FontSize', 13, 'Color', 'r');

yline(0.5, 'r--', 'LineWidth', 1.0);
text(1e-2, 0.35, 'Límite inf. zona crítica (×0.5)', 'FontSize', 13, 'Color', 'r');

% Líneas verticales de frontera de zonas
xline(50,   'k:', 'FontSize', 10);
xline(1e4,  'k:', 'FontSize', 10);

set(gca, 'XScale', 'log', 'YScale', 'log');

% Configurar ejes (tamaños grandes)
xlabel('Energía (keV)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Respuesta relativa r (E)', 'FontSize', 18, 'FontWeight', 'bold');

title(sprintf('ICRU 95 – a = [%.4f, %.4f, %.4f]', ...
      x_best(1), x_best(2), x_best(3)), ...
      'FontSize', 18, 'FontWeight', 'bold');

legend({'Zona térmica–50 keV [0.2–8.0]', ...
        'Zona epitérmica 50 keV–10 MeV [0.5–2.0]', ...
        'Zona >10 MeV [0.2–2.0]', ...
        'Respuesta continua', ...
        'Bines IEC (puntos de validación)'}, ...
       'Location', 'best', 'FontSize', 13);

xlim([E_bins(1)*1e3*0.5, 2.5e6]);
ylim([0.05, 12]);
grid on; grid minor;


%%  11. GRÁFICA AUXILIAR: Funciones de respuesta y h*(10) en los bines

figure('Name', 'Funciones de respuesta en los bines', ...
       'Color', 'w', 'Position', [80 640 950 380]);
loglog(E_grid*1e3, h_g,  'k-',  'LineWidth', 2.5, 'DisplayName', 'h*(10) ICRP-74'); hold on;
loglog(E_grid*1e3, Rb_g, 'b-',  'LineWidth', 1.5, 'DisplayName', 'R_{bare}');
loglog(E_grid*1e3, Rc_g, 'r-',  'LineWidth', 1.5, 'DisplayName', 'R_{core}');
loglog(E_grid*1e3, Rh_g, 'm-',  'LineWidth', 1.5, 'DisplayName', 'R_{hepb}');
% Marcar los bines
loglog(E_bins*1e3, h_bins, 'k+', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Bines IEC');
xlabel('Energía (keV)', 'FontSize', 11);
ylabel('Respuesta / Coeficiente de conversión', 'FontSize', 11);
title('Funciones de respuesta HENSA++ y h*(10) ICRP-95', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 9);
grid on; grid minor;
xline(50,  'k:', '50 keV',  'FontSize', 8, 'LabelVerticalAlignment', 'bottom');
xline(1e4, 'k:', '10 MeV', 'FontSize', 8, 'LabelVerticalAlignment', 'bottom');



%%  12. GRÁFICA h*10 y varias soluciones

A = [Rb_g, Rc_g, Rh_g]; 
b = h_g;
x_normal = A \ b;

% --- Bins de límites para graficar (sin el corte en 1 MeV) ---
E_bins_plot = [
    2.53e-8; 1.00e-7; 1.00e-6; 5.00e-6; 1.00e-5; 5.00e-5;
    1.00e-4; 5.00e-4; 1.00e-3; 1.00e-2; 5.00e-2;
    1.00e-1; 2.50e-1; 5.00e-1; 1.00e+0; 2.00e+0; 5.00e+0;
    1.00e+1; 1.40e+1; 2.00e+1; 5.00e+1; 7.50e+1;
    1.00e+2; 1.50e+2; 2.00e+2; 5.00e+2; 7.50e+2; 1.00e+3;
];

% Filtrar solo al rango válido común (sin el corte en 1 MeV)
mask_plot = (E_bins_plot >= Emin_valid) & (E_bins_plot <= Emax_valid);
E_bins_plot = E_bins_plot(mask_plot);

% Recalcular límites normativos para estos bins
lim_lo_plot = zeros(length(E_bins_plot), 1);
lim_hi_plot = zeros(length(E_bins_plot), 1);
for k = 1:length(E_bins_plot)
    if E_bins_plot(k) < 0.050
        lim_lo_plot(k) = 0.2;  lim_hi_plot(k) = 8.0;
    elseif E_bins_plot(k) <= 10.0
        lim_lo_plot(k) = 0.5;  lim_hi_plot(k) = 2.0;
    else
        lim_lo_plot(k) = 0.2;  lim_hi_plot(k) = 2.0;
    end
end

% Interpolar h*(10) en estos bins
h_bins_plot = interp_ll(vEh, vh, E_bins_plot);


figure('Name', 'Funciones solución y h*(10)', ...
       'Color', 'w', 'Position', [80 640 950 380]);

loglog(E_grid*1e3, h_g,  'k-',  'LineWidth', 2.5, 'DisplayName', 'h* '); hold on;
loglog(E_grid*1e3, A*x_pos, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Solución encontrada');

% Límites graficados desde el inicio del rango válido
loglog(E_bins_plot*1e3, h_bins_plot .* lim_hi_plot, 'g--', ...
'LineWidth', 1.5, 'DisplayName', 'Límites IEC');
loglog(E_bins_plot*1e3, h_bins_plot .* lim_lo_plot, 'g--', ...
'LineWidth', 1.5, 'DisplayName', 'Límite inferior', 'HandleVisibility', 'off');

set(gca, 'FontSize', 13, 'FontWeight', 'bold');

% Bines de optimización (los que sí entraron en el ajuste)
loglog(E_bins*1e3, h_bins, 'k+', 'MarkerSize', 10, 'LineWidth', 2, ...
'DisplayName', 'Bines IEC (optimización)');

xlabel('Energía (keV)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('h* (pSv·cm²)', 'FontSize', 18, 'FontWeight', 'bold');
title('Coeficientes h* y combinación lineal de respuestas', 'FontSize', 18, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 13);

grid on; grid minor;


%%  13. GRÁFICA DIAGNOSTICO

figure;
subplot(2,1,1);
loglog(E_grid, (Rb_g*x_best(1) + Rc_g*x_best(2) + Rh_g*x_best(3)), 'r', E_grid, h_g, 'k');
title('Respuesta Ajustada vs Referencia');
subplot(2,1,2);
semilogx(E_grid, (Rb_g*x_best(1) + Rc_g*x_best(2) + Rh_g*x_best(3)) .* phi_g .* E_grid, 'r'); hold on;
semilogx(E_grid, h_g .* phi_g .* E_grid, 'k');
title('Contribución a la Dosis');