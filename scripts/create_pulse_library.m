%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Create Pulse Library
%
% Read glottal flow pulses and speech parameters extracted by GlottHMM
% vocoder and construct a pulse library.
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

function create_pulse_library()

    clc
    clear
    close all
    version = 1.7;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % CONSTRUCT LIBRARY
    extract_params = 1;
    construct_library = 1;
    save_to_file = 1;
    
    % PULSE LIBRARY NAME AND SPEECH FILES
    pulse_library_name = 'pulselib1';
    wavdir = '../wav/';
    number_of_files = 1;
    
    % Remove outlier pulses (optional)
    remove_outliers = 1;
    outlier_limit = 0.5;
    interpolation_samples = 800;

    % Shift excitation peaks to center (refinement)
    shift_excitation_peaks = 1;

    % PATHS AND PARAMS
    homedir = [pwd '/../'];
    glotthmm_analysis = '../bin/Analysis';
    config_def = '../cfg/config_default';
    config_usr = '../cfg/config_default';
    param_names = {'lsf','gain','hnr','lsfsource','pulselengths',...
                   'pulses','rspulses','pulsepos','infofile'};
    
    % Select additional parameters (h1h2, naq, harmonics, waveform)
    %param_names = [param_names,'h1h2','naq','harmonics','waveform'];
    use_h1h2 = 0;
    use_naq = 0;
    use_harmonics = 0;
    use_waveform = 0;

    % UTILITIES
    plot_pulses = 1;
    plot_params_only = 0;
    index_pulses = 0;
    cluster_pulses = 0;
    number_of_clusters = 100;
    remove_extracted_files = 0;
    user = 'User Name';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
       
    % Print
    disp('<><><><><><><><><><><><><><><><><><><><><><><><>');
    disp(['<>      Pulse library constructor v. ' num2str(version) '      <>']);
    disp('<><><><><><><><><><><><><><><><><><><><><><><><>');
    disp(' ');

    % Defaults
    ascii = 1;
    binary = 2;


    
    %%%%%%%%%%%%%%%%%%%%%% Extract parameters %%%%%%%%%%%%%%%%%%%%%%%
    if extract_params == 1
        
        % Load file names
        files = dir([wavdir '*.wav']);
        if isempty(files)
            disp(['No audio files found in directory ' wavdir]);
            return;
        end

        % Check that no synthesized files are included
        discard = ones(size(files));
        for i = 1:length(files)
            if strcmp(files(i).name(end-6:end),'syn.wav') == 1
                discard(i) = 0;
            end 
        end
        if sum(discard) < length(files)
            disp('Already synthesized files found in the data (*.syn.wav)')
            disp(['  Discarding synthesized data' char(10)]);
            files = files(discard == 1);
        end
       
        disp('Extracting parameters:')

        % Analysis
        number_of_files = min(number_of_files,length(files));
        for i = 1:number_of_files
            
            % Create folder for parameter files
            if ~exist([homedir 'parameters'],'dir')
                system(['mkdir ' homedir 'parameters']);
            end
            
            % Show progress
            disp(['  Extracting file ' files(i).name ' (' int2str(i) '/' int2str(number_of_files) ')']);

            % Extract
            command_extract = [glotthmm_analysis ' ' wavdir files(i).name ' ' config_def ' ' config_usr];
            [~,result] = system(command_extract);
            if ~isempty(regexp(result,'[Ee][Rr][Rr][Oo][Rr]', 'once')) || ~isempty(regexp(result,'[Ff][Aa][Uu][Ll][Tt]', 'once'))
                disp(result);
                return
            end
            %system(command_extract);pause;

            % Copy parameter files and delete others
            for j = 1:length(param_names)
                if strcmp(param_names{j},'infofile') == 1  % Infofile does not have _PULSELIB extension
                    filename1 = [files(i).name(1:end-4) '.' param_names{j}];
                    filename2 = [files(i).name(1:end-4) '.PULSELIB.' param_names{j}];
                    system(['cp ' wavdir filename1 ' ' homedir 'parameters/' filename2]);
                else
                    filename = [files(i).name(1:end-4) '.PULSELIB.' param_names{j}];
                    system(['cp ' wavdir filename ' ' homedir 'parameters/' filename]);
                end
            end
            if remove_extracted_files == 1
                system(['rm ' wavdir files(i).name(1:end-4) '.*[!*.wav]']);
            end
        end
    end

    


    
    

    %%%%%%%%%%%%%%%%%%%%%%%%% Index pulses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if index_pulses == 1

        disp(' ')
        disp('Indexing pulses...')

        % Create folder for pulse position files
        if ~exist([homedir 'pulse_libraries/' pulse_library_name],'dir')
            system(['mkdir ' homedir 'pulse_libraries/' pulse_library_name]);
        end
        
        % Read files
        files = dir([wavdir '*.wav']);
        if isempty(files)
            disp(['No audio files found in directory ' wavdir]);
            return;
        end
        
        % Check that no synthesized files are included
        discard = ones(size(files));
        for i = 1:length(files)
            if strcmp(files(i).name(end-6:end),'syn.wav') == 1
                discard(i) = 0;
            end 
        end
        if sum(discard) < length(files)
            disp('Already synthesized files found in the data (*.syn.wav)')
            disp(['  Discarding synthesized data' char(10)]);
            files = files(discard == 1);
        end

        % Index each pulse individually
        INDEX = 0;
        for i = 1:min(number_of_files,length(files))

            % Read number of pulses
            params = load([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.infofile']);
            npulses = params(10);
            data_format = params(15);

            % Read pulse position info
            pulsepos = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.pulsepos'],data_format);

            % Group pulse indices by frame and create individual indices
            infoind = 1;
            ind = 1;
            pulsepos_info = [];
            while(ind < npulses+1)
                str = int2str(INDEX);
                INDEX = INDEX + 1;
                k = 1;
                while(ind+k < npulses+1 && pulsepos(ind) == pulsepos(ind+k))
                    str = [str ',' int2str(INDEX)];
                    INDEX = INDEX + 1;
                    k = k + 1;
                end
                pulsepos_info{infoind} = [int2str(pulsepos(ind)) '     ' str];
                ind = ind + k;
                infoind = infoind + 1;
            end
            
            % Save to file
            filename = [homedir 'pulse_libraries/' pulse_library_name '/pulsepos'];
            if data_format == ascii
                fid = fopen(filename,'wt');
                for k = 1:length(pulsepos_info)
                    fprintf(fid,'%s\n',pulsepos_info{k});
                end
            elseif data_format == binary
                if exist(filename,'file')
                    system(['rm ' filename]);
                end
                fid = fopen(filename,'a');
                for k = 1:length(pulsepos_info)
                    fwrite(fid,pulsepos_info{1},'double');
                end
            end
            fclose(fid);
        end
        disp('  Finished indexing.');
    end

    
    
    
    

    %%%%%%%%%%%%%%%%%%%%%% Construct pulse library %%%%%%%%%%%%%%%%%%%%%%%
    if construct_library == 1

        disp(' ');
        disp('Constructing pulse library:');
        
        % Load file names
        files = dir([wavdir '*.wav']);
        if isempty(files)
            disp(['No audio files found in directory ' wavdir]);
            return;
        end
        
        % Check that no synthesized files are included
        discard = ones(size(files));
        for i = 1:length(files)
            if strcmp(files(i).name(end-6:end),'syn.wav') == 1
                discard(i) = 0;
            end 
        end
        if sum(discard) < length(files)
            disp([char(10) 'Already synthesized files found in the data (*.syn.wav)']  )
            disp(['  Discarding synthesized data' char(10)]);
            files = files(discard == 1);
        end

        % Read basic params
        params = load([homedir 'parameters/' files(1).name(1:end-4) '.PULSELIB.infofile']);
        harm_degree = params(9);
        hnr_degree = params(8);
        tilt_degree = params(5);
        lsf_degree = params(4);
        fs = params(14);
        pulsemaxlen = floor(params(11)/1000*fs);
        rspulsemaxlen = floor(params(12)/1000*fs);
        waveform_samples = params(13);
        data_format = params(15);

        % Allocate memory for variables (faster speed)
        l = 0;
        for i = 1:min(number_of_files,length(files))
            params = load([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.infofile']);
            l = l + params(10);
        end
        pulses = zeros(l*pulsemaxlen,1);
        rspulses = zeros(l*rspulsemaxlen,1);
        pulselengths = zeros(l,1);
        gain = zeros(l,1);
        lsf = zeros(l*lsf_degree,1);
        hnr = zeros(l*hnr_degree,1);
        lsfsource = zeros(l*tilt_degree,1);
        if use_h1h2;h1h2 = zeros(l,1);end;
        if use_naq;naq = zeros(l,1);end;
        if use_waveform;waveform = zeros(l*waveform_samples,1);end;
        if use_harmonics;harmonics = zeros(l*harm_degree,1);end;
        
        % Read and save data
        pulseindex = 0;
        for i = 1:min(number_of_files,length(files))

            % Progress
            disp(['  Processing parameters of ' files(i).name ' (' int2str(i) '/' int2str(number_of_files) ')']);

            % Read number of pulses
            params = load([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.infofile']);
            npulses = params(10);

            % Read pulse and parameter data
            pulses_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.pulses'],data_format);
            rspulses_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.rspulses'],data_format);
            pulselengths_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.pulselengths'],data_format);
            lsf_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.lsf'],data_format);
            gain_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.gain'],data_format);
            hnr_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.hnr'],data_format);
            lsfsource_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.lsfsource'],data_format);
            if use_h1h2;h1h2_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.h1h2'],data_format);end;
            if use_naq;naq_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.naq'],data_format);end;
            if use_harmonics;harmonics_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.harmonics'],data_format);end;
            if use_waveform;waveform_temp = load_data([homedir 'parameters/' files(i).name(1:end-4) '.PULSELIB.waveform'],data_format);end;
            
            % Save pulse and parameter data
            pulses(pulseindex*pulsemaxlen+1:(pulseindex+npulses)*pulsemaxlen) = pulses_temp;
            rspulses(pulseindex*rspulsemaxlen+1:(pulseindex+npulses)*rspulsemaxlen) = rspulses_temp;
            pulselengths(pulseindex+1:pulseindex+npulses) = pulselengths_temp;
            gain(pulseindex+1:pulseindex+npulses) = gain_temp;
            lsf(pulseindex*lsf_degree+1:(pulseindex+npulses)*lsf_degree) = lsf_temp;
            hnr(pulseindex*hnr_degree+1:(pulseindex+npulses)*hnr_degree) = hnr_temp;
            lsfsource(pulseindex*tilt_degree+1:(pulseindex+npulses)*tilt_degree) = lsfsource_temp;
            if use_h1h2;h1h2(pulseindex+1:pulseindex+npulses) = h1h2_temp;end;
            if use_naq;naq(pulseindex+1:pulseindex+npulses) = naq_temp;end;
            if use_harmonics;harmonics(pulseindex*harm_degree+1:(pulseindex+npulses)*harm_degree) = harmonics_temp;end;
            if use_waveform;waveform(pulseindex*waveform_samples+1:(pulseindex+npulses)*waveform_samples) = waveform_temp;end;

            % Increment pulse index
            pulseindex = pulseindex + npulses;
        end

        % Shift excitation peak of the pulses to the center
        if shift_excitation_peaks == 1
            for i = 1:length(pulselengths)
                plen = pulselengths(i);
                pulse_temp = pulses((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+plen);
                center = find(diff(pulse_temp) == max(diff(pulse_temp)))-1;
                shift = round(plen/2) - center;
                lim = 30;
                if center > plen/2-lim && center < plen/2+lim
                    pulse_temp = circshift(pulse_temp,[shift 1]);
                end
                pulses((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+plen) = pulse_temp;
            end
        end
        
        % Shift excitation peak of the resampled pulses to the center
        if shift_excitation_peaks == 1
            for i = 1:length(pulselengths)
                pulse_temp = rspulses((i-1)*rspulsemaxlen+1:i*rspulsemaxlen);
                center = find(diff(pulse_temp) == max(diff(pulse_temp)))-1;
                shift = round(rspulsemaxlen/2) - center;
                lim = 30;
                if center > rspulsemaxlen/2-lim && center < rspulsemaxlen/2+lim
                    pulse_temp = circshift(pulse_temp,[shift 1]);
                end
                rspulses((i-1)*rspulsemaxlen+1:i*rspulsemaxlen) = pulse_temp;
            end
        end
        
        % Remove pulses that are far from average (possibly bad pulses)
        if remove_outliers == 1

            % Interpolate to constant length
            pulses_cl = zeros(length(pulselengths),interpolation_samples);
            for i = 1:length(pulselengths)
                p = pulses((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+pulselengths(i));
                pulses_cl(i,:) = interp1(1:length(p),p,linspace(1,length(p),interpolation_samples),'cspline');
            end

            % Remove outliers (pulses that are far from mean
            remove = zeros(size(pulselengths));
            m = mean(pulses_cl,1);            
            rmse = zeros(size(pulselengths));
            for i = 1:length(pulselengths)
                rmse(i) = sum((pulses_cl(i,:)-m).^2);
            end
            remove(rmse > outlier_limit) = 1;
            
            % Select n-dimensional parameters
            pulses = remoutliers(pulses,remove);
            rspulses = remoutliers(rspulses,remove);
            lsf = remoutliers(lsf,remove);
            hnr = remoutliers(hnr,remove);
            lsfsource = remoutliers(lsfsource,remove);
            if use_harmonics;harmonics = remoutliers(harmonics,remove);end;
            if use_waveform;waveform = remoutliers(waveform,remove);end;
            
            % Select 1-dimensional parameters
            pulselengths = pulselengths(remove == 0);
            gain = gain(remove == 0);
            if use_h1h2;h1h2 = h1h2(remove == 0);end;
            if use_naq;naq = naq(remove == 0);end;

            % Save the total number of pulses
            params(10) = length(pulselengths);
        end
    end





    %%%%%%%%%%%%%%%%%%%%%% Save pulse library %%%%%%%%%%%%%%%%%%%%%%%
    if save_to_file == 1

        % Progress
        disp(' ');
        disp('Saving pulse library...');
        
        % Set params to infofile
        infofile = params;

        % Save
        if ~exist([homedir 'pulse_libraries'],'dir')
            system(['mkdir ' homedir 'pulse_libraries']);
        end
        if ~exist([homedir 'pulse_libraries/' pulse_library_name],'dir')
            system(['mkdir ' homedir 'pulse_libraries/' pulse_library_name]);
        end
        for i = 1:length(param_names)
            if strcmp(param_names{i},'pulsepos') == 1  % Exclude pulse_pos
                continue;
            end
            filename = [homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.' param_names{i}];
            temp1 = 'temp';
            temp2 = param_names{i};
            evalc(sprintf('%s=%s',temp1,temp2));
            
            % Write to file
            if data_format == ascii
                fid = fopen(filename,'wt');
                fprintf(fid,'%1.9f\n',temp);
            elseif data_format == binary
                fid = fopen(filename,'w');
                fwrite(fid,temp,'double');
            end
            fclose(fid);
        end
        
        % Write data to infofile (always as ASCII)
        filename = [homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.infofile'];
        fid = fopen(filename,'wt');
        fprintf(fid,'%1.9f\n',infofile);
        fclose(fid);
        
        % Write additional data to readme (always as ASCII)
        filename = [homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.readme'];
        fid = fopen(filename,'wt');
        fprintf(fid,'%s\n',['Pulse library: ' pulse_library_name]);
        for i = 1:length(pulse_library_name)+15;fprintf(fid,'%s','-');end;
        fprintf(fid,'\n%s\n',['Created: ' datestr(now)]);
        fprintf(fid,'%s\n',['Created by: ' user]);
        fprintf(fid,'%s\n',['Pulse library constructor version: ' num2str(version)]);
        fprintf(fid,'%s\n',['Number of pulses: ' int2str(length(pulselengths))]);
        fprintf(fid,'%s\n','Files:');
        for i = 1:min(number_of_files,length(files))
            fprintf(fid,' %s\n',files(i).name);
        end
        
        % Copy configs
        system(['cp ' config_def ' ' homedir 'pulse_libraries/' pulse_library_name '/config']);
        system(['cp ' config_usr ' ' homedir 'pulse_libraries/' pulse_library_name '/config_usr']);
        disp(['  Created pulse library "' pulse_library_name '" (' int2str(number_of_files) ' wav-file(s), ' int2str(length(pulselengths)) ' pulses)']);
    end

    
     %%%%%%%%%%%%%%%%%%%%%% Cluster library %%%%%%%%%%%%%%%%%%%%%%%%
    if cluster_pulses == 1

        % Check pulse library
        if ~exist([homedir 'pulse_libraries/' pulse_library_name],'dir')
            disp(' ');
            disp(['Pulse library "' pulse_library_name '" does not exist']);
            return
        else
            disp(' ');
            disp('Clustering pulse library...')
        end

        % Load infofile
        params = load([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.infofile']);
        data_format = params(15);
        hnr_degree = params(8);
        tilt_degree = params(5);
        lsf_degree = params(4);
        fs = params(14);
        pulsemaxlen = round(params(11)/1000*fs);
        rspulsemaxlen = round(params(12)/1000*fs);
        harm_degree = params(9);
        waveform_samples = params(13);
        
        % Load data
        pulses = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.pulses'],data_format);
        rspulses = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.rspulses'],data_format);
        pulselengths = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.pulselengths'],data_format);
        gain = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.gain'],data_format);
        hnr = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.hnr'],data_format);
        lsfsource = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.lsfsource'],data_format);
        lsf = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.lsf'],data_format);
        if use_naq;naq = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.naq'],data_format);end
        if use_h1h2;h1h2 = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.h1h2'],data_format);end
        if use_waveform;waveform = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.waveform'],data_format);end
        if use_harmonics;harmonics = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.harmonics'],data_format);end
        
        % Create NxP data matrix
        N = length(pulselengths);
        p = [1 1 hnr_degree tilt_degree lsf_degree];
        if use_harmonics;p = [p harm_degree];end
        if use_naq;p = [p 1];end
        if use_h1h2;p = [p 1];end
        P = sum(p);
        data = zeros(N,P);
        for i = 1:N
            data(i,1) = pulselengths(i);
            data(i,sum(p(1:1))+1:sum(p(1:2))) = gain((i-1)*p(2)+1:i*p(2));
            data(i,sum(p(1:2))+1:sum(p(1:3))) = hnr((i-1)*p(3)+1:i*p(3));
            data(i,sum(p(1:3))+1:sum(p(1:4))) = lsfsource((i-1)*p(4)+1:i*p(4));
            data(i,sum(p(1:4))+1:sum(p(1:5))) = lsf((i-1)*p(5)+1:i*p(5));
            if use_harmonics;data(i,sum(p(1:5))+1:sum(p(1:6))) = harmonics((i-1)*p(6)+1:i*p(6));end
            if use_naq;data(i,sum(p(1:6))+1:sum(p(1:7))) = naq((i-1)*p(7)+1:i*p(7));end
            if use_h1h2;data(i,sum(p(1:7))+1:sum(p(1:8))) = h1h2((i-1)*p(8)+1:i*p(8));end
        end
    
        % K-means clustering
        k = number_of_clusters;
        [idx,~,~,d] = kmeans(data, k);

        % Select one pulse closest to the centroid from each cluster
        min_dist = 2*max(max(d))*ones(1,k);
        min_ind = zeros(1,k);
        for i = 1:N
            dist = d(i,idx(i));
            if dist < min_dist(idx(i))
                min_dist(idx(i)) = dist;
                min_ind(idx(i)) = i;
            end
        end

        % Save selected pulses to pulse library
        pulses_k = zeros(k*pulsemaxlen,1);
        rspulses_k = zeros(k*rspulsemaxlen,1);
        pulselengths_k = zeros(k,1);
        gain_k = zeros(k,1);
        lsf_k = zeros(k*lsf_degree,1);
        hnr_k = zeros(k*hnr_degree,1);
        lsfsource_k = zeros(k*tilt_degree,1);
        if use_naq;naq_k = zeros(k,1);end
        if use_h1h2;h1h2_k = zeros(k,1);end
        if use_harmonics;harmonics_k = zeros(k*harm_degree,1);end
        if use_waveform;waveform_k = zeros(k*waveform_samples,1);end
        for i = 1:k
            C = min_ind(i);
            pulses_k((i-1)*pulsemaxlen+1:i*pulsemaxlen) = pulses((C-1)*pulsemaxlen+1:C*pulsemaxlen);
            rspulses_k((i-1)*rspulsemaxlen+1:i*rspulsemaxlen) = rspulses((C-1)*rspulsemaxlen+1:C*rspulsemaxlen);
            pulselengths_k((i-1)+1:i) = pulselengths((C-1)+1:C);
            gain_k((i-1)+1:i) = gain((C-1)+1:C);
            lsf_k((i-1)*lsf_degree+1:i*lsf_degree) = lsf((C-1)*lsf_degree+1:C*lsf_degree);
            hnr_k((i-1)*hnr_degree+1:i*hnr_degree) = hnr((C-1)*hnr_degree+1:C*hnr_degree);
            lsfsource_k((i-1)*tilt_degree+1:i*tilt_degree) = lsfsource((C-1)*tilt_degree+1:C*tilt_degree);
            if use_naq;naq_k((i-1)+1:i) = naq((C-1)+1:C);end
            if use_h1h2;h1h2_k((i-1)+1:i) = h1h2((C-1)+1:C);end
            if use_harmonics;harmonics_k((i-1)*harm_degree+1:i*harm_degree) = harmonics((C-1)*harm_degree+1:C*harm_degree);end
            if use_waveform;waveform_k((i-1)*waveform_samples+1:i*waveform_samples) = waveform((C-1)*waveform_samples+1:C*waveform_samples);end
            
        end
        params(10) = k;

        % Save files
        pulse_library_name = [pulse_library_name '_' int2str(k) 'cl'];
        if ~exist([homedir 'pulse_libraries/' pulse_library_name],'dir')
            system(['mkdir ' homedir 'pulse_libraries/' pulse_library_name]);
        end
        for i = 1:length(param_names)
            if strcmp(param_names{i},'pulsepos') == 1 || strcmp(param_names{i},'infofile') == 1 % Exclude pulse_pos and infofile
                continue;
            end
            filename = [homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.' param_names{i}];
            temp1 = 'temp';
            temp2 = [param_names{i} '_k'];
            evalc(sprintf('%s=%s',temp1,temp2));
            
            % Write to file
            if data_format == ascii
                fid = fopen(filename,'wt');
                fprintf(fid,'%1.9f\n',temp);
            elseif data_format == binary
                fid = fopen(filename,'w');
                fwrite(fid,temp,'double');
            end
            fclose(fid);
        end
        
        % Write infofile always as ASCII
        filename = [homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.infofile'];
        fid = fopen(filename,'wt');
        fprintf(fid,'%1.9f\n',params);
        fclose(fid);
        
        disp(['  Created clustered pulse library "' pulse_library_name '" (' int2str(k) ' pulses)']);
    end


    
    


    %%%%%%%%%%%%%%%%%%%%%% Plot pulse library %%%%%%%%%%%%%%%%%%%%%%%
    if plot_pulses == 1

        % Check pulse library
        if ~exist([homedir 'pulse_libraries/' pulse_library_name],'dir')
            disp(' ');
            disp(['Pulse library "' pulse_library_name '" does not exist']);
            return
        else
            disp(' ');
            disp('Plotting pulse library...');
            disp('  Press any key to go to the next pulse');
            disp('  Press Ctrl-c to stop plotting');
        end

        % Load infofile
        params = load([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.infofile']);
        data_format = params(15);
        fs = params(14);
        
        % Load data
        pulses = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.pulses'],data_format);
        rspulses = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.rspulses'],data_format);
        pulselengths = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.pulselengths'],data_format);
        gain = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.gain'],data_format);
        hnr = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.hnr'],data_format);
        lsfsource = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.lsfsource'],data_format);
        lsf = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.lsf'],data_format);
        if use_naq;naq = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.naq'],data_format);end;
        if use_h1h2;h1h2 = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.h1h2'],data_format);end;
        if use_harmonics;harmonics = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.harmonics'],data_format);end;
        if use_waveform;waveform = load_data([homedir 'pulse_libraries/' pulse_library_name '/' pulse_library_name '.waveform'],data_format);end;
        pulsemaxlen = floor(params(11)/1000*fs);
        
        % Plot
        if plot_params_only == 0
            for i = 1:length(pulselengths)

                % Pulses
                plen = pulselengths(i);
                pulse = pulses((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen + plen);

                % Parameters
                cur_phnr = hnr((i-1)*params(8)+1:i*params(8));
                cur_ptilt = lsfsource((i-1)*params(5)+1:i*params(5));
                cur_plsf = lsf((i-1)*params(4)+1:i*params(4));
                if use_waveform;cur_waveform = waveform((i-1)*params(13)+1:i*params(13));end;
                if use_harmonics;cur_pharm = harmonics((i-1)*params(9)+1:i*params(9));end;

                % Modify VT parameters
                nfft = 1024;
                a_vt = lsf2poly(cur_plsf);
                [h_vt,w_vt] = freqz(1,a_vt,nfft,fs);
                H_vt = db(h_vt);
                lambda_vt = params(6);
                if lambda_vt ~= 0
                    omega = linspace(0,pi,nfft);
                    w_vt = (fs/2)*(omega + 2*atan((-lambda_vt*sin(omega))./(1 + lambda_vt*cos(omega))))/pi;
                end
                
                % Modify GL parameters
                a_gl = lsf2poly(cur_ptilt);
                [h_gl,w_gl] = freqz(1,a_gl,nfft,fs);
                H_gl = db(h_gl);
                lambda_gl = params(7);
                if lambda_gl ~= 0
                    omega = linspace(0,pi,nfft);
                    w_gl = (fs/2)*(omega + 2*atan((-lambda_gl*sin(omega))./(1 + lambda_gl*cos(omega))))/pi;
                end
                
                % Plot
                figure(1)
                clf

                subplot(2,2,1)
                plot(pulse)
                axis([1 length(pulse) min(pulse) max(pulse)])
                title('Waveform')
                ylabel('Amplitude')
                xlabel('Time (samples)')

                subplot(2,2,2)
                plot(w_vt,H_vt,'b')
                axis([0 fs/2 -50 60]);
                title('Vocal tract spectrum')
                ylabel('Magnitude (dB)')
                xlabel('Frequency (Hz)')
                
                subplot(2,2,3)
                plot(cur_phnr)
                axis([1 length(cur_phnr) min(hnr) max(hnr)])
                set(gca,'XTick',[1 2 3 4 5]);
                set(gca,'XTickLabel',{'1','2','3','4','5'});
                title('HNR')
                ylabel('HNR (dB)')
                xlabel('ERB bin')

                subplot(2,2,4)
                plot(w_gl,H_gl)
                axis([0 fs/2 -50 60]);
                title('Voice source spectrum')
                ylabel('Magnitude (dB)')
                xlabel('Frequency (Hz)')
                
                if use_harmonics
                    figure(2)
                    plot(cur_pharm)
                    axis([1 length(cur_pharm) min(harmonics) max(harmonics)])
                    title('Harmonics')
                    ylabel('Difference from H0 (dB)')
                    xlabel('Harmonic index')
                end
                
                if use_waveform
                    figure(3)
                    plot(cur_waveform)
                    axis([1 params(13) min(waveform) max(waveform)])
                    title(['Resampled waveform (' int2str(params(13)) ')'])
                    ylabel('Amplitude')
                    xlabel('Time (samples)')
                end

                pause
            end
        else
            % Plot parameters at once
            harm_degree = params(9);
            hnr_degree = params(8);
            tilt_degree = params(5);
            lsf_degree = params(4);
            cur_phnr = zeros(length(pulselengths),hnr_degree);
            cur_ptilt = zeros(length(pulselengths),tilt_degree);
            cur_plsf = zeros(length(pulselengths),lsf_degree);
            if use_harmonics;cur_pharm = zeros(length(pulselengths),harm_degree);end;
            for i = 1:length(pulselengths)
                cur_phnr(i,:) = hnr((i-1)*params(8)+1:i*params(8));
                cur_ptilt(i,:) = lsfsource((i-1)*params(5)+1:i*params(5));
                cur_plsf(i,:) = lsf((i-1)*params(4)+1:i*params(4));
                if use_harmonics;cur_pharm(i,:) = harmonics((i-1)*params(9)+1:i*params(9));end
            end

            figure(1)
            plot(pulselengths);
            title('Pulse lengths')

            figure(2)
            plot(pulses);
            hold on
            plot(rspulses+5,'r');
            legend('Pulses','Pulses (resampled)')

            figure(3)
            plot(gain);
            title('Gain')

            figure(4)
            plot(cur_phnr)
            title('HNR')

            figure(5)
            plot(cur_ptilt)
            title('Tilt')

            figure(6)
            plot(cur_plsf)
            title('LSF')
            
            if use_harmonics == 1
                figure(7);
                plot(cur_pharm)
                title('Harmonics')
            end

            if use_naq == 1
                figure(8)
                plot(naq)
                title('NAQ')
            end
            
            if use_h1h2 == 1
                figure(9)
                plot(h1h2)
                title('H1H2')
            end
        end
    end
end
 




   

% Function for removing outlier pulses and parameters
function par = remoutliers(par,rem)
    par = reshape(par,length(par)/length(rem),length(rem));
    par = par(:,rem == 0);
    par = reshape(par,size(par,1)*size(par,2),1);
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


