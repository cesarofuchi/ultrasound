%%

clc; close all; clear all;

%%
% bifasico ================================================================
samples = 5000;
channels = 1;

% filepath = 'E:\2_Doutorado\Atividades\Rocking_Flow\raw_data\';

% ref_tubo_cheio
% filepath = 'E:\2_Doutorado\Atividades\Rocking_Flow\raw_data\tubo_cheio\';
% filename{1,1}  = '2020-12-17 19_03_53_,372tuboCheio.xml';
% filename{1,2}  = 'tubo_cheio';

% ref_tubo_vazio
% filepath = 'E:\2_Doutorado\Atividades\Rocking_Flow\raw_data\0_5_rad_s\';
% filename{1,1}  = '2020-12-17 19_01_11_,281testeSync2_0_5rad.xml';
% filename{1,2}  = 'oscilation _vel = 0.5';

% oscilation _vel = 0,5 rad/s
% filepath = 'E:\2_Doutorado\Atividades\Rocking_Flow\raw_data\0_5_rad_s\';
% filename{1,1}  = '2020-12-17 19_01_11_,281testeSync2_0_5rad.xml';
% filename{1,2}  = 'oscilation _vel = 0.5';

% oscilation _vel = 1,0 rad/s
% filepath = 'E:\2_Doutorado\Atividades\Rocking_Flow\raw_data\1_0_rad_s\';
% filename{1,1}  = '2020-12-17 18_54_55_,577testeSync2_1_0rad.xml';
% filename{1,2}  = 'oscilation _vel = 1.0';

% oscilation _vel = 1,5 rad/s
filepath = 'E:\2_Doutorado\Atividades\Rocking_Flow\raw_data\1_5_rad_s\';
filename{1,1}  = '2020-12-17 18_51_12_,250testeSync3_1_5rad.xml';
filename{1,2}  = 'oscilation _vel = 1.0';

% =========================================================================

% =========================================================================
% fill doppler properties

temperatura = 25;

% velocidade do som na água por Lubers and Graaaf
% c = 1404.3 + 4.7 * temperatura - 0.04 * temperatura^2;
c = 1332; % m/s -> para o óleo

dec = 1; % fator de decimação
ns = [16 32 64 128 256]; % janela espacial (amostras) de medição de velocidade
nc = [16 32 64 128 256]; % janela temporal (pulsos) de medição de velocidade
dif = 0; 
ovt = 1; ovs = 1; % step (passo) janelas temporal e espacial

% =========================================================================
%%

% Create and initialize usObj and usdObj for each experimental point

for i = 1 : 1 : size(filename,1) % percorre todos os pontos experimentais
    
    file = [filepath filename{i}];
    
    usObj(i) = sensor.Ultrasonic.loadData(file,channels,samples);
    usObj(i).fc = 4e6;
    usObj(i).channels = 1;

    usObj(i); % show Ultrasonic Object

    % crop the data of this experiment

    % data = usObj(i).data(1800:4550,1:30000);
    data = usObj(i).data(1900:4550,1:30000);
    % usObj(i).workData = data;
        
    
    % Matched Filter ======================================================
    % Is the optimal linear filter for maximizing the SNR in the presence of
    % additive stochastic noise.

    impulse_response = sin(2*pi*usObj(i).fc*(0:1/usObj(i).fs:(usObj(i).cycles)/usObj(i).fc));

    % Hamming windowing
    impulse_response = impulse_response.*hamming(max(size(impulse_response)))';

    b = flipud(impulse_response);
    dataf1 = filter(b,1,data);
    
%     figure;
%     subplot 211
%     hold all    
%     plot(dataf1(:,1))
%     plot(data(:,1))
%     xlabel('samples')
%     ylabel('amplitude(V)')
%     legend('Matched Filter signal','Raw signal')   
%     title('Ultrasound timeseries')
%     grid on
% 
%     subplot 212
%     plot(b)
%     title(['Impulse response with ' num2str(usObj(i).cycles) 'cycles for Matched Filter'])
    
    % create doppler object ===============================================
    usObj(i).workData = dataf1;
    usdObj(i) = sensor.UltrasonicDoppler(usObj(i),ns,nc,ovs,ovt);
    
    % Show demodulation results ===========================================

%     figure
%     subplot 311
%     hold all
%     plot(real(usdObj(i).iq(:,1)))
%     plot(imag(usdObj(i).iq(:,1)))
%     plot(abs(usdObj(i).iq(:,1)))
%     legend('real','imag','abs')
%     grid on
%     title('Hilbert')
% 
%     subplot 312
%     hold all
%     plot(real(usdObj(i).iq(:,1)))
%     plot(imag(usdObj(i).iq(:,1)))
%     legend('real','imag')
%     grid on
%     title('Complete Demodulation')
% 
%     subplot 313
%     hold all
%     plot(real(usdObj(i).iq(1:dec:end,1)))
%     plot(imag(usdObj(i).iq(1:dec:end,1)))
%     legend('real','imag')
%     grid on
%     title(['Complete Demodulation decimation:' num2str(dec)])
    
end % usObj and usdObj created and initialized
    
%% Bubble Detection ===================================================

%     clear dbubble dtemp
% 
%     dtemp = (abs(dataf1));
%     threshold = 5e3;
% 
%     for i = 1 : 1 : size(dtemp,2) % todas as colunas (ondas ou pulsos)
% 
%             idx_temp = find(dtemp(:,i) > threshold);
% 
%             if size(idx_temp,1) > 0
%                 idx = idx_temp(1,1);
%                 dbubble(:,i) = dataf1(:,i);
%                 dbubble(idx:end,i) = 0;
%             else
%                 dbubble(:,i) = dataf1(:,i);
%             end
% 
%     end
% =====================================================================
    


%% ACM ====================================================================
disp('============================');
disp('ACM Running'); tic;

for i = 1 : 1 : size(usdObj,2) % varre todos os objetos (pontos experimentais)
    for j = 1 : 1 : size(nc,2) % varre todos os nc
        for k = 1 : 1 : size(ns,2) % varre todos os ns
            
            subsampling = usdObj(i).ns(k);
            vr_acm(i,j,k) = usdObj(i).AutoCorrelation(c,subsampling,j,k);
            
        end
    end
end
p_time = toc; disp([num2str(p_time) ' seconds']);
save(fullfile('D:\Área de Trabalho\OnePhase\velocities', 'vr_acm.mat'), 'vr_acm');

%% CCM ====================================================================
disp('============================');
disp('CCM Running'); tic;

% dbubble = zeros(size(data)); % variável usada para exp. bifásico

for i = 1 : 1 : size(usdObj,2) % varre todos os objetos (pontos experimentais)
    for j = 1 : 1 : size(nc,2) % varre todos os nc
        for k = 1 : 1 : size(ns,2) % varre todos os ns
            
            d_temp = usdObj(i).usObj.workData;
            vr_ccm(i,j,k) = usdObj(i).TimeShiftEstimator(double(d_temp),c,j,k);
            
        end
    end
end
p_time = toc; disp([num2str(p_time) ' seconds']);
save(fullfile('D:\Área de Trabalho\OnePhase\velocities', 'vr_ccm.mat'), 'vr_ccm');

%% EAM ====================================================================
disp('============================');
disp('EAM Running'); tic;

Np_temp{1} = [-2 -1 0 1 2]; Np = Np_temp{1};
% Np_temp{2} = [-3 -2 -1 0 1 2 3];
% Np_temp{3} = [-4 -3 -2 -1 0 1 2 3 4];

% thld = ([4 5 4 7 10 10 5 8 9 12 12 12 5 12 14 14 14 15 4.5 15 10 20 20 25])*1e3;
thld = 20e3;

% esse 'for' é utilizado para testar diferentes ordens de Np (range de busca)
% for l = 1 : 1 : size(Np_temp,2)
% txt = ['NpOrd = ', num2str(l)];
% disp(txt);

% Np = Np_temp{l};

clear vr_eam;

for i = 1 : 1 : size(usdObj,2) % varre todos os objetos (pontos experimentais)
    txt = ['usdObj = ', num2str(i)];
    disp(txt)    
    % detecção da bolha ===================================================
    clear dbubble dtemp

    dtemp = abs(usdObj(i).usObj.workData);
    %threshold = thld(i);
    threshold = thld;

    for z = 1 : 1 : size(dtemp,2) % todas as colunas (ondas ou pulsos)

            idx_temp = find(dtemp(:,z) > threshold);

            if size(idx_temp,1) > 0
                idx = idx_temp(1,1);
                dbubble(:,z) = usdObj(i).usObj.workData(:,z);
                dbubble(idx:end,z) = 0;
            else
                dbubble(:,z) = usdObj(i).usObj.workData(:,z);
            end
    end
%     dbbl_temp(i,:,:) = dbubble;
%     b_temp(:,:) = dbbl_temp(i,:,:);
%     figure; imagesc(flipud(abs(hilbert(b_temp(i,:,:))))); colormap jet; colorbar;
    % =====================================================================
    
    % "for's" aninhados usados para teste de varredura em (ns) e (nc).
    % "resolução do pixel"
    % for j = 1 : 1 : size(nc,2) % varre todos os nc
        % for k = 1 : 1 : size(ns,2) % varre todos os ns
            j = 2; k = 2;
            vr_eam = usdObj(i).ExtAutoCorrelation(dbubble,c,dif,Np,j,k);
            %vr_eam(i,j,k) = usdObj(i).ExtAutoCorrelation(dbubble,c,dif,Np,j,k);            
        % end    
    % end   
end

p_time = toc; disp([num2str(p_time) ' seconds']);

%     if l == 1
%         save(fullfile('D:\Área de Trabalho\TwoPhase\velocities', 'vr_eam_NpOrd2.mat'), 'vr_eam');
%     elseif l == 2
%         save(fullfile('D:\Área de Trabalho\TwoPhase\velocities', 'vr_eam_NpOrd3.mat'), 'vr_eam');
%     else
%         save(fullfile('D:\Área de Trabalho\TwoPhase\velocities', 'vr_eam_NpOrd4.mat'), 'vr_eam');
%     end

% end

%% Staggered ==============================================================

% converte os dados para como se tivesse sido amostrados no método
% staggered m/n
disp('============================');
disp('STT Running'); tic;

m = 2; n = m + 1;

for i = 1 : 1 : size(usdObj,2) % varre todos os objetos (pontos experimentais)

    noise = usObj(i).data(2000:5000,1:20000);

    iqn = hilbert(noise);
    tn  = (0:1/usObj(i).fs:size(iqn,1)*(1/usObj(i).fs)-(1/usObj(i).fs))';
    tn  = repmat(tn,1,size(iqn,2));
    
    noised = iqn.*exp(-1i*2*pi*(usObj(i).fc)*tn);
    
    datad = usdObj(i).iq.*exp(-1i*2*pi*(usObj(i).fc)*usdObj(i).t);

    for z = 1 : size(datad,2)
        if z == 1
            indiceST(z) = 1;
        else
            if (rem(z,2) == 0) % desloca m
                indiceST(z) = indiceST(z-1) + m;
            else
                indiceST(z) = indiceST(z-1) + n;            
            end
            if indiceST(z) > size(datad,2)
                break;
            end
        end
    end

    indiceST = indiceST(1:end-1);
    datadST = datad(:,indiceST);
    clear datad;
    
    for j = 1 : 1 : size(nc,2) % varre todos os nc
        for k = 1 : 1 : size(ns,2) % varre todos os ns
            
            noisesample = noised(1:usdObj(i).ns(k),1:usdObj(i).nc(j));
            vnoise = sum(noisesample,1);
            auto0noise = (sum(abs(vnoise(1:usdObj(i).nc(j))).^2)+sum(abs(vnoise(1:usdObj(i).nc(j))).^2)) ...
                         / (2*(length(vnoise(1:usdObj(i).nc(j)))));
            
            T2 = n/usObj(i).fprf; T1 = m/usObj(i).fprf;
                     
            vr_stt(i,j,k) = usdObj(i).StaggeredTriggerTorres(datadST,c,T1,T2,(auto0noise),m,n,j,k);
            
        end
    end    
end
p_time = toc; disp([num2str(p_time) ' seconds']);
save(fullfile('D:\Área de Trabalho\OnePhase\velocities', 'vr_stt.mat'), 'vr_stt');

%% Pós-Processamento

% vr = load('E:\2_Doutorado\Atividades\Doppler\Exps\Experimento 20_02_2020\Experimento 20_02_2020 Resultados\Figs_and_velocities\TwoPhase\velocities\vr_eam_NpOrd2.mat');
% vr = vr.vr_eam;
% vr = vr.vr_acm;
% vr = vr.vr_ccm;
% vr = vr.vr_stt;

vr = vr_eam;

angle = 10;

d_real = 26; % pipe diameter (mm)

% in samples
pipe_ini = 1753; % 0 mm
pipe_end = 4590; % 25.4 mm
size_pipe = pipe_end - pipe_ini + 1;

% how many samples represent 1 mm
factor = size_pipe/d_real;

% portion of the wave to be analyzed
data_ini = 1800;
data_end = 4550;

y_ini = (data_ini - pipe_ini)/factor;
y_end = (data_end - data_ini + 1)/factor;

method = 'EAM';

%%
% fig_path = (['D:\Área de Trabalho\TwoPhase\',method]);
% fig_path = ([fig_path,'\double_medfilt\JG_100_JL_150']);
% fig_path = ([fig_path,'\single_medfilt\medfilt2D(3x3)\JG_100_JL_150']);

% i = 24;

% i = 1;  j = 5; k = 3; x_time = 10;  excel_sheet = 'JG_025_JL_025';
% i = 2;  j = 3; k = 3; x_time = 2.5; excel_sheet = 'JG_025_JL_050';
% i = 3;  j = 2; k = 4; x_time = 1.5; excel_sheet = 'JG_025_JL_075';
% i = 4;  j = 2; k = 4; x_time = 1;   excel_sheet = 'JG_025_JL_100';
% i = 5;  j = 2; k = 5; x_time = 1;   excel_sheet = 'JG_025_JL_125';
% i = 7;  j = 4; k = 2; x_time = 10;  excel_sheet = 'JG_050_JL_025';
% i = 8;  j = 3; k = 3; x_time = 2.5; excel_sheet = 'JG_050_JL_050';
% i = 9;  j = 3; k = 4; x_time = 1.5; excel_sheet = 'JG_050_JL_075';
% i = 11; j = 2; k = 5; x_time = 1;   excel_sheet = 'JG_050_JL_125';
% i = 13; j = 4; k = 3; x_time = 10;  excel_sheet = 'JG_075_JL_025';
% i = 14; j = 3; k = 3; x_time = 2.5; excel_sheet = 'JG_075_JL_050';
% i = 15; j = 3; k = 4; x_time = 2.5; excel_sheet = 'JG_075_JL_075';
% i = 16; j = 3; k = 4; x_time = 1.5; excel_sheet = 'JG_075_JL_100';
% i = 19; j = 4; k = 3; x_time = 10;  excel_sheet = 'JG_100_JL_025';
% i = 20; j = 4; k = 4; x_time = 4;   excel_sheet = 'JG_100_JL_050';
% i = 21; j = 3; k = 4; x_time = 2;   excel_sheet = 'JG_100_JL_075';

% vr_temp(:,:) = vr(i,:,:); method = 'EAM';
vr_temp = vr_eam; method = 'EAM';

JG_JL = usObj(i).filename;
JG_JL_full = extractAfter(JG_JL,"JG_");
JG_JL_full = ['JG_' JG_JL_full];

% thld = ([4 5 4 7 10 10 5 8 9 12 12 12 5 12 14 14 14 15 4.5 15 10 20 20 25])*1e3;

thld = 20e3;

% loop's utilizados para testes de varredura ==============================
% for j = 1 : 1 : size(nc,2) 
%    for k = 1 : 1 : size(ns,2) 

        % corrige a velocidade (na direção do fluxo)
        flowr = vr_temp.vel/sind(angle); 
        flow = (-flowr);

        % aplica filtro de mediana (suavização)
        %vv = medfilt1(flow,3,[],1);
        %vv = medfilt2(vv,[3 3]); filt = 'medfilt1D(3)_medfilt2D(3x3)';
        vv = medfilt1(flow,3,[],1);
        vv = medfilt2(vv,[3 3]); filt = 'medfilt2D(3x3)';        
        
        vv = vv.*vr_temp.msk;
        
        % congiguration of velocity maps axis =============================
        clear x_vv y_vv
        
        x_vv = linspace(0,size(usObj(i).workData,2)*(1/usObj(i).fprf),size(vv,2));
        y_vv = linspace(y_ini,y_end,size(vv,1));
        % =================================================================
        
        dtemp = abs(usdObj(i).usObj.workData);
        % threshold = thld(i);
        threshold = thld;

        % bubble detection (on raw data) ==================================
        for z = 1 : 1 : size(dtemp,2) % todas as colunas (ondas ou pulsos)

                idx_temp = find(dtemp(1:end,z) > threshold);

                if size(idx_temp,1) > 0
                    idx = idx_temp(1,1);
                    dbubble(:,z) = usdObj(i).usObj.workData(:,z);
                    dbubble(idx:end,z) = 0;
                else
                    dbubble(:,z) = usdObj(i).usObj.workData(:,z);
                end
        end

%%        
        vtemp = flipud(vv);
        for z = 1 : 1 : size(wdw,1) % varre todas as janelas de um velocity map

            % extracts windows from the velocity map to a temporary variable
            wdw_temp = vtemp(wdw(z,2):wdw(z,4),wdw(z,1):wdw(z,3));

            % evaluates row(spatial) and column(temporal) means of a window
            mean_col_wdw = mean(wdw_temp,1);
            mean_row_wdw = mean(wdw_temp,2);
            
            % max velocity of each window (from profile)
            max_mean_row(z,:) = max(mean_row_wdw);
            
            % temporal and spatial durations
            wdw_duration(z,:) = x_vv(wdw(z,3)) - x_vv(wdw(z,1));
            wdw_width(z,:) = y_vv(wdw(z,4)) - y_vv(wdw(z,2));
            
            % FIGURES =====================================================
            % spatial axis configuration
%             y_ini_row = (wdw(z,2) * ns(k))/factor + (y_ini/2);
%             y_end_row = (wdw(z,4) * ns(k))/factor;
            y_ini_row = (y_vv(wdw(z,2)));
            y_end_row = (y_vv(wdw(z,4)));
            y_axis_row = linspace(y_end_row,y_ini_row,size(mean_row_wdw,1));
            
            % time axis configuration
%             tx_ini = (wdw(z,1) * nc(j) * 1/usdObj(i).usObj.fprf);
%             tx_end = (wdw(z,3) * nc(j) * 1/usdObj(i).usObj.fprf);
            tx_ini = (x_vv(wdw(z,1)));
            tx_end = (x_vv(wdw(z,3)));
            x_axis_col = linspace(tx_ini,tx_end,size(mean_col_wdw,2));

            figure; set(gcf, 'Position',  [500, 45, 225, 640]);
            
            subplot 311 % Row Mean ========================================
            plot(flipud(mean_row_wdw),y_axis_row); grid on;
            set(gca, 'XDir','reverse');
            set(gca, 'YDir','reverse');
            title('Row Mean');
            xlabel('velocity (m/s)');
            ylabel('distance (mm)');
%             xlim([0.5 0.9])
%             ylim([0 21])
            
            subplot 312 % Column Mean =====================================
            mean_col_fit = fit(x_axis_col',flipud(mean_col_wdw)','poly1');
            coefs_mean_col(z,:) = (mean_col_fit.p1);
            
            fig(1) = plot(x_axis_col,flipud(mean_col_wdw)); grid on; hold on;
            fig(2) = plot(mean_col_fit);
            legend([fig(2)],num2str(mean_col_fit.p1));
            
            % set(gca, 'XDir','reverse');
            % set(gca, 'YDir','reverse');
            title('Column Mean');
            xlabel('time (s)');
            ylabel('velocity (m/s)');
%             ylim([0.5 0.9])
            
            subplot 313 % Histogram =======================================
            histogram(wdw_temp); grid on;
            title('Velocity Distribution');
            xlabel('velocity (m/s)');
            ylabel('frequency');
%             xlim([0 1.3])
%             ylim([0 200])
            % xticks([0 0.5 1.0 1.3])
            
            % END FIGURES =================================================
            
        end % varredura de todas as janelas do mapa de velocidade
        
        % SAVE VAR's TO EXCEL FILE ========================================
            
%             excel_file = 'selected_velocities.xlsx';
%             
%             % build table with the vars to be saved
%             T = table(max_mean_row,coefs_mean_col,wdw_duration,wdw_width);
% 
%             writetable(T,excel_file,'Sheet',excel_sheet,'Range','A1');
%             
%             clear max_mean_row coefs_mean_col wdw_duration wdw_width
            
        % =================================================================
         
        % FIGURE ==========================================================
        
%         fig_name = ([method '_' JG_JL_full '_nc_' num2str(usdObj(i).nc(j))...
%                     '_ns_' num2str(usdObj(i).ns(k))]);
        
        fig = figure; set(gcf, 'Position',  [250, 80, 700, 400])
        
        subplot 311 % raw data after bubble detection
        % rescale the time axis (raw data - bubble detected) and plot =====
        x_ax = linspace(0,size(usObj(i).workData,2)*(1/usObj(i).fprf),size(usObj(i).workData,2));
        % xl = [x_ax(1) x_ax(end)];

        y_ax = linspace(y_ini,y_end,size(usObj(i).workData,1));
        % yl = ceil([y_ax(1) y_ax(end)]);
                
        %imagesc(xl,yl,flipud(abs(hilbert(dbubble))));
        imagesc(x_ax,y_ax,flipud(abs(hilbert(dbubble)))); ax1 = gca;
        cb = colorbar; colormap jet;
        cb.Title.String = "amplitude (a.u)";
        title('US Data - Bubble Detection');
        xlabel('time (s)');
        ylabel('distance (mm)');
        xlim([0 x_time])
        % =================================================================
        
        subplot 312 % velocity map
        % rescale the time axis (velocity map) and plot ===================
        x_vv = linspace(0,size(usObj(i).workData,2)*(1/usObj(i).fprf),size(vv,2));
        % x_vv = [x_vv(1) x_vv(end)]; 

        y_vv = linspace(y_ini,y_end,size(vv,1));
        % y_vv = ceil([y_vv(1) y_vv(end)]);
        
        imagesc(x_vv,y_vv,flipud(vv)); ax2 = gca;
        cb = colorbar; colormap jet;
        cb.Title.String = "velocity (m/s)";
        title([method, ' - ', JG_JL_full, ' - nc = ', num2str(nc(j)), ' - ', 'ns = ', num2str(ns(k))],'Interpreter', 'none');
        xlabel('time (s)');
        ylabel('distance (mm)');
        xlim([0 x_time]); 
        
        subplot 313 % flow profile
        plot(x_axis_col,flipud(mean_col_wdw)); grid on; ax3 = gca;
        title(['Flow Profile (red window - #' num2str(z) ')']);
        xlabel('time (s)');
        ylabel('velocity (m/s)');
        xlim([0 x_time]);
        % ylim([0.4 1.2])
        % =================================================================
        
        for z = 1 : 1 : (size(wdw,1)) % varre todas as janelas de um velocity map
            x1 = x_vv(wdw(z,1));
            y1 = y_vv(wdw(z,2));
            x2 = x_vv(wdw(z,3));
            y2 = y_vv(wdw(z,4));
            
            if z ~= size(wdw,1)
                %rectangle(ax1,'Position',[x1 y1 (x2-x1) (y2-y1)],'Curvature',0.2,'LineWidth',3);
                rectangle(ax2,'Position',[x1 y1 (x2-x1) (y2-y1)],'Curvature',0.2,'LineWidth',3);
            else
                %rectangle(ax1,'Position',[x1 y1 (x2-x1) (y2-y1)],'Curvature',0.2,'EdgeColor','r','LineWidth',3);
                %rectangle(ax2,'Position',[x1+0.01 y1-1.45 (x2-x1)-0.025 (y2-y1)+1.8],'Curvature',0.2,'EdgeColor','r','LineWidth',3);
            end
        end
        
%%               
%         saveas(gca, fullfile(fig_path, fig_name), 'tif');
            
%         close(fig);
        
        % =================================================================
%    end
% end