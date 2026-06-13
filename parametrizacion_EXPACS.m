function interactive_phi_fit()
% Herramienta visual para ajustar phi(E) al espectro EXPACS Javalambre
% usando el modelo analítico de FRUIT para las zonas térmicas y epitérmicas 
% y la parametrización de Gordon para Fast y High Energy.
addpath('C:\Users\eduji\Desktop\Master\TFM\Data');

%% 1. CARGA DE DATOS
[vEphi,  vPhi ]   = read_2cols('EXPACSJavalambre_phi.txt');
[vE_Eph, v_EPhi]  = read_2cols('EXPACSJavalambre_Ephi.txt');
[vEh,   vh]       = read_2cols('h_Emax_ICRU95');
[vER, vRbare, vRcore, vRhepb] = read_4cols('Response_Functions_HENSApp_OAJ_Phase0.dat');

%% 2. MALLA DE ENERGÍA
Emin   = max(vEphi(1),   vE_Eph(1));
Emax   = min(vEphi(end), vE_Eph(end));
E_plot = logspace(log10(Emin), log10(Emax), 700);
% Referencia EXPACS en E_plot  (E * phi)
EPhi_ref = exp(interp1(log(vE_Eph), log(v_EPhi), log(E_plot), 'linear','extrap'));

%% 3. RANGOS FÍSICOS
To_range = [1e-10, 1e-6 ];   % Temperatura Maxwell-Boltzmann [MeV]
Ed_range = [1e-9,  1e-5 ];   % Corte inferior epitérmico [MeV]
                             
%% 4. ESPECIFICACIONES DE SLIDERS
% Columnas: { campo, etiqueta, esLog, min, max, valorInicial }
specs = { ...
  'k',  'k  – escala global',              true,  1e-6,  1e2,    4.209e-02 ; ...
  'p1', 'p1 – Int. Térmica (<0.35eV)',     false, 0.01,  0.98,   0.1944    ; ...
  'To', 'To – temp. térmica (MeV)',        true,  1e-10, 1e-6,   2.513e-8  ; ...
  'p2', 'p2 – Int. Epitérmica (->0.1MeV)', false, 0.01,  0.98,   0.2366    ; ...
  'b',  'b  – pendiente epitérmica',       false, -0.45, 0.45,   0.0503    ; ...
  'bp', 'bp – corte epitérmico (MeV)',     false, 0.01,  180,    2.0       ; ...
  'Ed', 'Ed – corte Cd (MeV)',             true,  1e-9,  1e-5,   1.723e-7  ; ...
  'p3', 'p3 – Int. Rápida (->20MeV)',      false, 0.01,  0.98,   0.3162    ; ...
  'bf', 'β2 – beta (Fast Gordon)',         false, 0.05,  1.5,    0.2389    ; ...
  'gf', 'γ2 – gamma (Fast Gordon)',        false, -3.0,  3.0,   -0.5919    ; ...
  'p4', 'p4 – Int. Alta Energía (>20MeV)', false, 0.01,  0.98,   0.2528    ; ...
  'bh', 'β1 – beta (High Gordon)',         false, 0.05,  1.5,    0.6261    ; ...
  'gh', 'γ1 – gamma (High Gordon)',        false, -2.0,  10.0,   5         ; ...
};
nS = size(specs, 1);

%% 5. PARÁMETROS INICIALES
P0 = struct();
for i = 1:nS
    P0.(specs{i,1}) = specs{i,6};
end

%% 6. CONSTRUIR FIGURA
fig = uifigure('Name','Ajuste interactivo phi(E) – Javalambre', ...
               'Position',[30 30 1500 950], 'Color',[0.95 0.95 0.97]);

% Ejes del espectro
ax = uiaxes(fig, 'Position', [20 80 800 840]);
ax.XScale = 'log';  ax.YScale = 'linear';
xlabel(ax,'E (MeV)');  ylabel(ax,'E \cdot \phi(E)');
grid(ax,'on');  hold(ax,'on');

% Panel de sliders
pnl = uipanel(fig, 'Title','Panel de Control de Parámetros', ...
              'Position',[840 90 640 830], ... 
              'BackgroundColor',[0.95 0.95 0.97], ...
              'FontSize',12, 'FontWeight','bold', ...
              'Scrollable', 'on');

% Etiqueta RMS 
lbl_rms = uilabel(fig, 'Position', [20 45 300 30], ...
                  'Text', 'RMSE = —', ...
                  'FontSize', 13, 'FontWeight', 'bold');

% Botón Reset
btn_rst = uibutton(fig, 'push', 'Text', '↺ Reset', ...
                   'Position', [1100 25 80 45], 'FontSize', 12);

% Botón Optimizar (Derecha)
btn_opt = uibutton(fig, 'push', 'Text', '▶  Optimizar (fminsearch)', ...
                   'Position', [1190 25 290 45], 'FontSize', 12, ...
                   'BackgroundColor', [0.9, 1, 0.9], 'FontWeight', 'bold');

%% 7. CREAR SLIDERS
row_h  = 60;    % Más espacio vertical entre sliders
top_y  = 880;   % Punto de inicio superior
sliders   = struct();
vallabels = struct();

for i = 1:nS
    fld   = specs{i,1};
    lbl   = specs{i,2};
    isLog = specs{i,3};
    vmin  = specs{i,4};
    vmax  = specs{i,5};
    vinit = specs{i,6};
    
    y = top_y - i*row_h;
    
    if isLog
        sLim = [log10(vmin), log10(vmax)];
        sVal = log10(vinit);
    else
        sLim = [vmin, vmax];
        sVal = vinit;
    end
    
    uilabel(pnl, 'Position',[15 y+28 350 20], ...
            'Text', lbl, 'FontSize', 10.5, 'FontWeight', 'bold');
            
    vl = uilabel(pnl, 'Position',[450 y+28 160 20], ...
                 'HorizontalAlignment','right', ...
                 'FontSize', 11, 'FontWeight','bold', 'FontColor', [0 0.4 0.7]);
    
    if isLog
        vl.Text = sprintf('%.3e', vinit);
    else
        vl.Text = sprintf('%.4f', vinit);
    end
    vallabels.(fld) = vl;
    
    sl = uislider(pnl, 'Position',[25 y+15 580 3], ...
                  'Limits', sLim, 'Value', sVal);
    sliders.(fld) = sl;
end

%% 8. ESTADO COMPARTIDO EN UserData
st.vER       = vER;
st.vRbare    = vRbare;
st.vRcore    = vRcore;
st.vRhepb    = vRhepb;
st.P         = P0;
st.E_plot    = E_plot;
st.EPhi_ref  = EPhi_ref;
st.specs     = specs;
st.sliders   = sliders;
st.vallabels = vallabels;
st.To_range  = To_range;
st.Ed_range  = Ed_range;
st.P0        = P0;
st.vEh       = vEh;      
st.vh        = vh;       
fig.UserData = st;

%% 9. DIBUJO INICIAL
refresh_plot(ax, lbl_rms, fig.UserData);

%% 10. CONECTAR CALLBACKS
for i = 1:nS
    fld   = specs{i,1};
    isLog = specs{i,3};
    sl    = sliders.(fld);
    vl    = vallabels.(fld);
    sl.ValueChangingFcn = @(~,evt) do_slider(evt.Value, fig, ax, lbl_rms, fld, isLog, vl);
    sl.ValueChangedFcn  = @(src,~) do_slider(src.Value, fig, ax, lbl_rms, fld, isLog, vl);
end
btn_opt.ButtonPushedFcn = @(~,~) do_optimize(fig, ax, lbl_rms);
btn_rst.ButtonPushedFcn = @(~,~) do_reset(fig, ax, lbl_rms);


end  % ── FIN interactive_phi_fit ────────────────────────────────────────────


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  CALLBACKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function do_slider(sval, fig, ax, lbl_rms, fld, isLog, vl)
    st = fig.UserData;
    if isLog
        val = 10^sval;
        vl.Text = sprintf('%.3e', val);
    else
        val = sval;
        vl.Text = sprintf('%.4f', val);
    end
    st.P.(fld) = val;
    fig.UserData = st;
    refresh_plot(ax, lbl_rms, st);
end

function do_optimize(fig, ax, lbl_rms)
    st = fig.UserData;
    a0 = pack_params(st.P, st.To_range, st.Ed_range);
    
    fprintf('Iniciando optimización con p1 = %.4f y k = %.4e\n', st.P.p1, st.P.k);
    cost = @(a) rms_cost(a, st.E_plot, st.EPhi_ref, ...
                         st.To_range, st.Ed_range);
    opts = optimset('Display','iter', ...
                    'MaxIter',5000, 'MaxFunEvals',25000, ...
                    'TolFun',1e-30, 'TolX',1e-30);
    fprintf('\n[Optimizar] Arrancando desde RMS = %.10f ...\n', sqrt(cost(a0)));
    [a_fit, fval] = fminsearch(cost, a0, opts);
    P_fit = unpack_params(a_fit, st.To_range, st.Ed_range);
    st.P  = P_fit;
    fig.UserData = st;
    sync_sliders(fig);
    refresh_plot(ax, lbl_rms, fig.UserData);
    fprintf('[Optimizar] Terminado.  RMS = %.10f\n\n', sqrt(fval));
    print_params(P_fit);
end

function do_reset(fig, ax, lbl_rms)
    st   = fig.UserData;
    st.P = st.P0;
    fig.UserData = st;
    sync_sliders(fig);
    refresh_plot(ax, lbl_rms, st);
end

function sync_sliders(fig)
    st        = fig.UserData;
    P         = st.P;
    specs     = st.specs;
    sliders   = st.sliders;
    vallabels = st.vallabels;
    for i = 1:size(specs,1)
        fld   = specs{i,1};
        isLog = specs{i,3};
        val   = P.(fld);
        sl    = sliders.(fld);
        vl    = vallabels.(fld);
        if isLog
            sv = log10(val);
        else
            sv = val;
        end
        sv = max(sl.Limits(1), min(sl.Limits(2), sv));
        sl.Value = sv;
        if isLog
            vl.Text = sprintf('%.3e', val);
        else
            vl.Text = sprintf('%.4f', val);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  VISUALIZACIÓN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refresh_plot(ax, lbl_rms, st)
    P   = st.P;
    E   = st.E_plot;
    ref = st.EPhi_ref;
    [Phi_tot, phi_t, phi_e, phi_f, phi_h] = compute_phi(E, P);
    EPhi_mod = Phi_tot .* E;

    umbral_minimo = 1e-20; 
    ok = (E < 1000) & (EPhi_mod > umbral_minimo) & (ref > 0) & isfinite(EPhi_mod);
    rms_val = sqrt(mean( (EPhi_mod(ok) - ref(ok)).^2 ));

    p_val = max([P.p1, P.p2, P.p3, P.p4], 1e-10);
    p_val = p_val / sum(p_val);
    
    cla(ax);
    hold(ax,'on');
    
    semilogx(ax, E, ref,      'r-',  'LineWidth', 2.2, 'DisplayName', 'EXPACS Javalambre');
    semilogx(ax, E, EPhi_mod, 'k-',  'LineWidth', 2.2, 'DisplayName', 'Ajuste total');
    semilogx(ax, E, phi_t.*E, 'b--', 'LineWidth', 1.4, 'DisplayName', sprintf('p1=%.3f Térmica', p_val(1)));
    semilogx(ax, E, phi_e.*E, 'g--', 'LineWidth', 1.4, 'DisplayName', sprintf('p2=%.3f Epitérmica', p_val(2)));
    semilogx(ax, E, phi_f.*E, 'm--', 'LineWidth', 1.4, 'DisplayName', sprintf('p3=%.3f Rápida', p_val(3)));
    semilogx(ax, E, phi_h.*E, 'c--', 'LineWidth', 1.4, 'DisplayName', sprintf('p4=%.3f Alta Energía', p_val(4)));
    
    set(ax, 'FontSize', 18, 'FontWeight', 'bold');
    legend(ax, 'Location', 'best', 'FontSize', 13);
    grid(ax, 'on');
    grid(ax, 'minor');

    xlabel(ax, 'E (MeV)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel(ax, 'E · ϕ̇  (E) (cm^{-2} · s^{-1})', ...
       'FontSize', 18, ...
       'FontWeight', 'bold');

    ax.Title.String = 'Ajuste del Espectro Neutrónico EXPACS - Javalambre';
    ax.Title.FontSize = 18;
    ax.Title.FontWeight = 'bold';
    lbl_rms.Text = sprintf('RMSE = %.6f', rms_val);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  FÍSICA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Phi_tot, phi_t, phi_e, phi_f, phi_h] = compute_phi(E, P)
    E = E(:)';
    % 1. Definir funciones de forma base puras
    u_t = (E ./ P.To^2) .* exp(-E ./ P.To);
    u_e = (1 - exp(-(E./P.Ed).^2)) .* max(E,1e-40).^(P.b-1) .* exp(-E./P.bp);
    u_f = exp( -P.bf .* (log(E)).^2 + P.gf .* log(E) );
    u_h = exp( -P.bh .* (log(E)).^2 + P.gh .* log(E) );
    
    % Elimina NaNs e Infinitos sustituyéndolos por cero
    u_t(~isfinite(u_t)) = 0; u_e(~isfinite(u_e)) = 0;
    u_f(~isfinite(u_f)) = 0; u_h(~isfinite(u_h)) = 0;
    
    % Normalización inicial de las formas
    phi_t = u_t ./ max(trapz(E, u_t), 1e-35);
    phi_e = u_e ./ max(trapz(E, u_e), 1e-35);
    phi_f = u_f ./ max(trapz(E, u_f), 1e-35);
    phi_h = u_h ./ max(trapz(E, u_h), 1e-35);
    
    % 2. Definir las máscaras para los 4 rangos físicos
    m1 = (E <= 3.5e-7);
    m2 = (E > 3.5e-7) & (E <= 0.1);
    m3 = (E > 0.1)    & (E <= 20);
    m4 = (E > 20);
    
    masks  = {m1, m2, m3, m4};
    shapes = {phi_t, phi_e, phi_f, phi_h};
    
    % 3. Matriz de transferencia M
    M = zeros(4,4);
    for i = 1:4
        idx = masks{i};
        if sum(idx) > 1
            for j = 1:4
                M(i,j) = trapz(E(idx), shapes{j}(idx));
            end
        end
    end
    
    % 4. Vector de integrales objetivo
    p_target = max([P.p1, P.p2, P.p3, P.p4]', 1e-10);
    p_target = p_target / sum(p_target); 
    
    % 5. Resolver el sistema M * w = p_target 
    opts = optimset('Display','off');
    w = lsqnonneg(M, p_target, opts);
    
    % 6. Construir espectro normalizado
    Phi_norm = w(1)*phi_t + w(2)*phi_e + w(3)*phi_f + w(4)*phi_h;
    Phi_norm = Phi_norm ./ max(trapz(E, Phi_norm), 1e-35);
    
    % 7. Escalar por el factor global k
    Phi_tot = P.k * Phi_norm;
    
    phi_t = P.k * w(1) * phi_t;
    phi_e = P.k * w(2) * phi_e;
    phi_f = P.k * w(3) * phi_f;
    phi_h = P.k * w(4) * phi_h;
    
    E_corte = 10000; 
    Phi_tot(E > E_corte) = 0;
    phi_t(E > E_corte) = 0;
    phi_e(E > E_corte) = 0;
    phi_f(E > E_corte) = 0;
    phi_h(E > E_corte) = 0;
end

function c = rms_cost(a, E, ref, To_r, Ed_r)
    try
        P = unpack_params(a, To_r, Ed_r);
        [Phi_tot,~,~,~,~] = compute_phi(E, P);
        EPhi_mod = Phi_tot .* E;
        ok = (EPhi_mod > 0) & (ref > 0) & isfinite(EPhi_mod);
        if sum(ok) < 10
            c = 1e10; return;
        end
        c = sqrt(mean( (EPhi_mod(ok) - ref(ok)).^2 ));
    catch
        c = 1e10;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  EMPAQUETADO / DESEMPAQUETADO DE PARÁMETROS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = pack_params(P, To_r, Ed_r)
    ll   = @(v,lo,hi) logit_log_fn(v, lo, hi);
    lln  = @(v,lo,hi) log( max((v-lo)/(hi-v), 1e-15) );
    
    w = max([P.p1, P.p2, P.p3, P.p4], 1e-9);
    w = w / sum(w);
    
    a = [ log(w(1)), log(w(2)), log(w(3)), log(w(4)), ... % 1-4
          lln(P.b,  -0.45, 0.45), ...                     % 5
          ll(P.bp,   0.01, 2.0 ), ...                     % 6
          log(P.k),               ...                     % 7
          ll(P.To,  To_r(1), To_r(2)), ...                % 8
          ll(P.Ed,  Ed_r(1), Ed_r(2)), ...                % 9
          lln(P.bf,  0.05, 1.5), ...                      % 10 (Fast beta)
          lln(P.gf, -3.0,  3.0), ...                      % 11 (Fast gamma)
          lln(P.bh,  0.05, 1.5), ...                      % 12 (High beta)
          lln(P.gh, -2.0,  5.0) ];                        % 13 (High gamma)
end

function P = unpack_params(a, To_r, Ed_r)
    sl   = @(x,lo,hi) siglog_fn(x, lo, hi);
    sln  = @(x,lo,hi) lo + (hi-lo) ./ (1 + exp(-x));
    
    w = exp(a(1:4)); w = w / sum(w);
    P.p1 = w(1); P.p2 = w(2); P.p3 = w(3); P.p4 = w(4);
    P.b  = sln(a(5), -0.45, 0.45);
    P.bp = sl( a(6),  0.01, 2.0 );
    P.k  = exp(a(7));
    P.To = sl( a(8),  To_r(1),  To_r(2));
    P.Ed = sl( a(9),  Ed_r(1),  Ed_r(2));
    P.bf = sln(a(10), 0.05,  1.5);
    P.gf = sln(a(11), -3.0,  3.0);
    P.bh = sln(a(12), 0.05,  1.5);
    P.gh = sln(a(13), -2.0,  5.0);
end

function y = logit_log_fn(v, lo, hi)
    v  = max(lo*1.000001, min(hi*0.999999, v));
    y  = log( (log(v) - log(lo)) / (log(hi) - log(v)) );
end

function y = siglog_fn(x, lo, hi)
    y = exp( log(lo) + (1/(1+exp(-x))) * (log(hi) - log(lo)) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  UTILIDADES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_params(P)
    w = max([P.p1, P.p2, P.p3, P.p4], 1e-10);  
    w = w / sum(w);
    fprintf('\n==================================================\n');
    fprintf('          RESULTADOS DEL AJUSTE (PHI)\n');
    fprintf('==================================================\n');
    fprintf('  k   (Escala global)      : %.4e\n', P.k);
    fprintf('--------------------------------------------------\n');
    fprintf('  TÉRMICA:\n');
    fprintf('    p1 (Int. Térmica)      : %.4f (%.1f%%)\n', w(1), w(1)*100);
    fprintf('    To (Temperatura)       : %.4e MeV\n', P.To);
    fprintf('--------------------------------------------------\n');
    fprintf('  EPITÉRMICA:\n');
    fprintf('    p2 (Int. Epitérmica)   : %.4f (%.1f%%)\n', w(2), w(2)*100);
    fprintf('    b  (Pendiente)         : %.4f\n', P.b);
    fprintf('    bp (Corte exponencial) : %.4f MeV\n', P.bp);
    fprintf('    Ed (Corte Cd)          : %.4e MeV\n', P.Ed);
    fprintf('--------------------------------------------------\n');
    fprintf('  FAST (Evaporación / Gordon j=2):\n');
    fprintf('    p3 (Int. Rápida)       : %.4f (%.1f%%)\n', w(3), w(3)*100);
    fprintf('    β2 (Beta Fast)         : %.4f\n', P.bf);
    fprintf('    γ2 (Gamma Fast)        : %.4f\n', P.gf);
    fprintf('--------------------------------------------------\n');
    fprintf('  HIGH ENERGY (Cascada / Gordon j=1):\n');
    fprintf('    p4 (Int. Alta Energía) : %.4f (%.1f%%)\n', w(4), w(4)*100);
    fprintf('    β1 (Beta High)         : %.4f\n', P.bh);
    fprintf('    γ1 (Gamma High)        : %.4f\n', P.gh);
    fprintf('==================================================\n\n');
end