%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% (1) Read pulse library and apply principal component analysis (PCA). Save
% principal component (PC) matrix and mean for the pulse library and PC
% weight matrix describing the pulses.
%
% (2) Extract speech parameters (optional if already extracted)
%
% (3) Read analysis parameters and construct PC weight matrices describing
% the voice source of each speech file.
%
% 19.1.2015
% (c) Tuomo Raitio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you use GlottHMM within an article, paper, report or any other work
% that you wish to publish, please cite it as:
% 
% T. Raitio, A. Suni, J. Yamagishi, H. Pulakka, J. Nurminen, M. Vainio,
% and P. Alku, "HMM-based speech synthesis utilizing glottal inverse
% filtering," IEEE Trans. on Audio, Speech and Lang. Proc., vol. 19, no. 1,
% pp. 153-165, Jan. 2011.
% 
% T. Raitio, A. Suni, H. Pulakka, M. Vainio, and P. Alku, "Utilizing
% glottal source pulse library for generating improved excitation signal
% for HMM-based speech synthesis," in Proc. ICASSP, 2011, pp. 4564-4567.
% 
% T. Raitio, A. Suni, L. Juvela, M. Vainio, and P. Alku, "Deep neural
% network based trainable voice source model for synthesis of speech with
% varying vocal effort," in Proc. of Interspeech, 2014, pp. 1969-1973.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pulselib_pca()

    clc
    clear
    close all

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Procedures: (1)/(2)/(3)
    perform_pca = 1;
    perform_analysis = 0; % Set 0 if Analysis is already performed
    convert_parameters_to_pcweights = 1;
    
    % Save data
    save_data = 1;
    
    % Shifting of the pulse peaks
    shift_pulse_peaks = 1;

    % (1) Settings for PCA
    pca_order = 12;
    pulselen = 800;
    pulselibname = '../pulse_libraries/pulselib1/pulselib1';
    select_pulse_closest_to_mean = 0;
    data_format = 1;
    
    % (2) Settings for parameter extraction
    number_of_files = 1;
    wavdir = '../wav/';
    glotthmm_analysis = '../bin/Analysis';
    config_file = '../cfg/config_default';
    
    % (3) Settings for parameter conversion
    filedir = '../wav/';
        
    % Plot options
    plot_principal_components = 1;
    plot_reconstructed_pulses = 1;
    
    % Plot option for each file to be converted
    plot_reconstructed_pulses_file = 0;
    plot_pc_matrix_file = 0;
    smooth_pc_matrix_file = 0;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % (1) PCA for pulse library
    if perform_pca == 1
        
        disp('Perform PCA for pulse library')

        % Load data and interpolate pulses to constant length
        info = load_data([pulselibname '.infofile'],data_format);
        data = load_data([pulselibname '.pulses'],data_format);
        lengths = load([pulselibname '.pulselengths']);
        pulsemaxlen = round(info(11)/1000*info(14));
        N = length(data)/pulsemaxlen;
        pulses = zeros(N,pulselen);

        % Interpolate
        for i = 1:N
            p = data((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+lengths(i));
            pulses(i,:) = interp1(1:length(p),p,linspace(1,length(p),pulselen),'cspline');
        end

        % Shift pulse excitation peak to the center
        if shift_pulse_peaks == 1
            for i = 1:N
                center = find(diff(pulses(i,:)) == max(diff(pulses(i,:))))-1;
                shift = round(pulselen/2) - center;
                lim = 30*pulselen/400;
                if center > pulselen/2-lim && center < pulselen/2+lim
                    pulses(i,:) = circshift(pulses(i,:),[1 shift]);
                end
            end
        end

        % Mean
        m = mean(pulses,1);

        % Select pulse closest to mean
        if select_pulse_closest_to_mean == 1
            proto_i = 0;
            min_err = 100000;
            for i = 1:N
                err = norm(m-pulses(i,:));
                if err < min_err
                   min_err = err;
                   proto_i = i;
                end
            end
            m = pulses(proto_i,:);
        end

        % Remove mean / closest natural pulse
        x = pulses-repmat(m,N,1);
        %x = pulses-repmat(m_natural,N,1);
        % Note that this is not exactly the same as performing PCA without
        % removing the mean (the first eigenvector being similar to mean)!
        
        % PCA
        [~,~,pc] = svd(x,'econ');

        % PC weights
        w = x*pc;
        
        % Reduce dimensions
        w = w(:,1:pca_order);
        pc = pc(:,1:pca_order);

        % Reconstruct pulses from reduced set of PCs and PC weights
        pulses_rc = repmat(m,N,1) + w*pc';
        
        % Extract noise part from pulses
        error = pulses - pulses_rc;
        
        % Save PCs, PC weights, and mean
        if save_data == 1
            filename_w = [pulselibname '.pca_w'];
            filename_pc = [pulselibname '.pca_pc'];
            filename_mean = [pulselibname '.pca_mean'];
            disp(['   Writing file ' filename_w]);
            disp(['   Writing file ' filename_pc]);
            disp(['   Writing file ' filename_mean]);
            if data_format == 1
                fid1 = fopen(filename_w,'wt');
                fid2 = fopen(filename_pc,'wt');
                fid3 = fopen(filename_mean,'wt');
                fprintf(fid1,'%1.9f\n',w');
                fprintf(fid2,'%1.9f\n',pc');
                fprintf(fid3,'%1.9f\n',m);
            elseif data_format == 2
                fid1 = fopen(filename_w,'w');
                fid2 = fopen(filename_pc,'w');
                fid3 = fopen(filename_mean,'w');
                fwrite(fid1,w','double');
                fwrite(fid2,pc','double');
                fwrite(fid3,m_natural,'double');
            end
            fclose(fid1);
            fclose(fid2);
            fclose(fid3);
        end
    end
    
    
    
    % (2) Extract parameters
    if perform_analysis == 1
        disp(' ')
        disp('Extracting parameters from files')
        files = dir([wavdir '*.wav']);
        files = files(1:min(length(files),number_of_files));
        for i = 1:length(files)
            disp(['   Extracting file ' files(i).name ' (' int2str(i) '/' int2str(number_of_files) ')'])
            command_extract = [glotthmm_analysis ' ' wavdir files(i).name ' ' config_file];
            [~,result] = system(command_extract);
            if ~isempty(regexp(result,'[Ee][Rr][Rr][Oo][Rr]', 'once')) || ~isempty(regexp(result,'[Ff][Aa][Uu][Ll][Tt]', 'once'))
                disp(result);
                return
            end
        end
    end
    
    
    
    % (3) Convert parameters
    if convert_parameters_to_pcweights == 1
        
        disp(' ')
        disp('Converting extracted pulses to corresponding principal component weights')

        % Start conversion for each speech file
        files = dir([filedir '*.PULSELIB.pulses']);
        for f = 1:length(files)
            
            % Read pulse files
            filename = files(f).name(1:end-16);
            info = load_data([filedir '/' filename '.infofile'],data_format);
            pulsedata = load_data([filedir '/' filename '.PULSELIB.pulses'],data_format);
            lengths = load_data([filedir '/' filename '.PULSELIB.pulselengths'],data_format);
            pulsepos = load_data([filedir '/' filename '.PULSELIB.pulsepos'],data_format);
            f0 = load_data([filedir '/' filename '.F0'],data_format);
            
            % Read PCA files
            pc = load_data([pulselibname '.pca_pc'],data_format);
            pc = reshape(pc,pca_order,pulselen);
            pc = pc';
            m = load_data([pulselibname '.pca_mean'],data_format);
            
            % Create pulse matrix
            pulsemaxlen = info(11)/1000*info(14);
            Nf = length(pulsedata)/pulsemaxlen;
            pulses_file = zeros(Nf,pulselen);
            for i = 1:Nf
                p = pulsedata((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+lengths(i));
                pulses_file(i,:) = interp1(1:length(p),p,linspace(1,length(p),pulselen),'cspline');
            
                % Shift pulse excitation peak to the center
                if shift_pulse_peaks == 1
                    center = find(diff(pulses_file(i,:)) == max(diff(pulses_file(i,:))))-1;
                    shift = round(pulselen/2) - center;
                    lim = 30*pulselen/400;
                    if center > pulselen/2-lim && center < pulselen/2+lim
                        pulses_file(i,:) = circshift(pulses_file(i,:),[1 shift]);
                    end
                end
            end

            % Convert pulses to principal component weights using the PCs
            % from the original pulse library
            x = pulses_file-repmat(m',Nf,1);
            w = x*pc;
            
            % Reduce dimensions
            w = w(:,1:pca_order);
            pc = pc(:,1:pca_order);

            % Create PC weight matrices for the speech files
            pcw = zeros(length(f0),pca_order);
            pulsepos = pulsepos + 3; % Indexing starts from 1 + 2 preframes
            for i = 1:length(pcw)
                idx = find(pulsepos == i);
                if ~isempty(idx)
                    
                    % Save the pulse pcw that is closest to the mean
                    err = zeros(size(idx));
                    for j = 1:length(idx)
                        err(j) = sum((pulses_file(idx(j),:)-m').^2);
                    end
                    pcw(i,:) = w(idx(err == min(err)),:);
                end
            end

            % Fill empty frames
            for i = 1:length(pcw)
                if f0(i) ~= 0 && sum(pcw(i,:)) == 0
                    j = 1;
                    stop_flag = 0;
                    while(stop_flag == 0)
                        ind = min(i + j,length(pcw));
                        if f0(ind) ~= 0 && sum(pcw(ind,:)) ~= 0
                        	pcw(i,:) = pcw(ind,:);
                            stop_flag = 1;
                        end
                        ind = max(i - j,1);
                        if f0(ind) ~= 0 && sum(pcw(ind,:)) ~= 0
                            pcw(i,:) = pcw(ind,:);
                             stop_flag = 1;
                        end
                        j = j + 1;
                    end
                end
            end

            % Smooth PC weight matrix
            if smooth_pc_matrix_file == 1
                for i = 1:pca_order
                    pcw(:,i) = smooth(pcw(:,i));
                end
            end

            % Save PC data
            if save_data == 1
                filename_pc = [filedir filename '.pca_w'];
                disp(['   Writing file ' filename_pc ' (' int2str(f) '/' int2str(number_of_files) ')'])
                if data_format == 1
                    fid = fopen(filename_pc,'wt');
                    fprintf(fid,'%1.9f\n',pcw');
                elseif data_format == 2
                    fid = fopen(filename_pc,'w');
                    fwrite(fid,pcw','double');
                end
                fclose(fid);
                disp(' ')
            end
            
            % Reconstruct pulses from the reduced set of PCs and PC weights
            if plot_reconstructed_pulses_file == 1
                disp('Plotting reconstruct pulses from the reduced set of PCs and PC weights...');
                disp('  Press any key to continue')
                disp('  Press Ctrl-c to stop plotting');
                figure(4)
                pulses_red = repmat(m',Nf,1) + w*pc';
                for i = 1:N
                    clf
                    plot(-pulses_file(i,:))
                    hold on
                    plot(-pulses_red(i,:),'r')
                    legend('Original',['PCA (' int2str(pca_order) ')'])
                    pause
                end
            end

            % Plot PC weight matrix
            if plot_pc_matrix_file == 1
                disp('Plotting PC weight matrix');
                disp('  Press any key to continue')
                disp('  Press Ctrl-c to stop plotting');
                figure(5)
                clf
                plot(pcw)
                xlabel('Frame')
                ylabel('Weight')
                h = title(['Weights of the ' int2str(pca_order) ' principal components for file ' filename '    (paused)']);
                set(h,'interpreter','none')
                xlim([1 length(pcw)])
                pause
            end
        end
    end

    % Plot mean and principal components
    if plot_principal_components == 1
        disp(['Plotting mean and principal components' char(10)]);
        
        % Plot mean
        figure(1);
        plot(m)
        title('Mean')
        
        % Plot principal components
        figure(2)
        for i = 1:pca_order
            a = ceil(sqrt(pca_order));
            subplot(a,a,i)
            plot(pc(:,i))
            xlim([1 pulselen]);
            title(['Component ' int2str(i)])
        end
    end

    % Plot reconstructed pulses
    if plot_reconstructed_pulses == 1
        disp('Plotting reconstructed pulses...');
        disp('  Press any key to continue')
        disp('  Press Ctrl-c to stop plotting');
        figure(3)
        for i = 1:N
            clf
            plot(pulses(i,:))
            hold on
            plot(pulses_rc(i,:),'r')
            legend('Original',['PCA (' int2str(pca_order) ')'],'Location','SouthEast')
            %plot(error(i,:),'m')
            %legend('Original',['PCA (' int2str(pca_order) ')'],'Error','Location','SouthEast')
            pause
        end
    end
end



% Function for reading data either in ascii (1) or binary (2) format
function data = load_data(filename,data_format)
    if data_format == 1
        data = load(filename);
    elseif data_format == 2
        file = fopen(filename,'r');
        data = fread(file,'double');
        fclose(file);
    end
end
