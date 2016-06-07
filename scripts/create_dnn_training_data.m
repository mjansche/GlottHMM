function create_dnn_training_data()

    clc
    clear
    close all

    % Pulse library name
    pulselib_dir = '../pulse_libraries/';
    pulselib_name = 'pulselib1';

    % Output directory
    out_dir = '../DNNTrain/';

    % Interpolated pulse length (samples)
    pulselen = 400;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Check pulse library
    if ~exist([pulselib_dir pulselib_name],'dir')
        disp(['Pulse library "' [pulselib_dir pulselib_name] '" does not exist']);
        return
    else
        disp('Creating data for DNN training...')
    end

    % Load infofile
    params = load([pulselib_dir pulselib_name '/' pulselib_name '.infofile']);
    data_format = params(15);
    hnr_degree = params(8);
    tilt_degree = params(5);
    lsf_degree = params(4);
    fs = params(14);
    pulsemaxlen = round(params(11)/1000*fs);
    rspulsemaxlen = round(params(12)/1000*fs);
    waveform_samples = params(13);

    % Load data
    pulses_raw = load_data([pulselib_dir pulselib_name '/' pulselib_name '.pulses'],data_format);
    pulselengths = load_data([pulselib_dir pulselib_name '/' pulselib_name '.pulselengths'],data_format);
    gain = load_data([pulselib_dir pulselib_name '/' pulselib_name '.gain'],data_format);
    hnr = load_data([pulselib_dir pulselib_name '/' pulselib_name '.hnr'],data_format);
    lsfsource = load_data([pulselib_dir pulselib_name '/' pulselib_name '.lsfsource'],data_format);
    lsf = load_data([pulselib_dir pulselib_name '/' pulselib_name '.lsf'],data_format);

    % Create Nxp training parameter matrix
    N = length(pulselengths);
    p = [1 1 hnr_degree tilt_degree lsf_degree];
    P = sum(p);
    data = zeros(N,P);
    for i = 1:N
        data(i,1) = 2*fs/pulselengths(i);
        data(i,sum(p(1:1))+1:sum(p(1:2))) = gain((i-1)*p(2)+1:i*p(2));
        data(i,sum(p(1:2))+1:sum(p(1:3))) = hnr((i-1)*p(3)+1:i*p(3));
        data(i,sum(p(1:3))+1:sum(p(1:4))) = lsfsource((i-1)*p(4)+1:i*p(4));
        data(i,sum(p(1:4))+1:sum(p(1:5))) = lsf((i-1)*p(5)+1:i*p(5));
    end

    % Interpolate pulses to constant length
    pulses = zeros(N,pulselen);
    for i = 1:N
        p = pulses_raw((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+pulselengths(i));
        pulses(i,:) = interp1(1:length(p),p,linspace(1,length(p),pulselen),'cspline');
    end

    % Shift pulse excitation peak to the center
    for i = 1:N
        center = find(diff(pulses(i,:)) == max(diff(pulses(i,:))))-1;
        if length(center) > 1
            continue;
        end
        shift = round(pulselen/2) - center;
        lim = 30*pulselen/400;
        if center > pulselen/2-lim && center < pulselen/2+lim
            pulses(i,:) = circshift(pulses(i,:),[1 shift]);
        end
    end

    % Normalize energy
    for i = 1:N
        pulses(i,:) = pulses(i,:)/sqrt(sum(pulses(i,:).^2));
    end

    % Rename variables
    data_in = data;
    data_out = pulses;

    % Save
    save([out_dir pulselib_name '.mat'],'data_in','data_out','-v7.3');
    disp('Done')
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


