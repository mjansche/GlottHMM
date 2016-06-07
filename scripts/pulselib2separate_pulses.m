clc
clear
close all

% Save pulses of a pulse library into individual files and also save the
% pulse closest to the mean (a prototype pulse)

% Settings
mean_pulse_length = 400; % samples

% Define pulse library
plibdir = '../pulse_libraries/';
plibname = 'pulselib1';

% Load data
disp('Saving pulses of a pulse library into individual files...')
pulses = load([plibdir plibname '/' plibname '.pulses']);
pulselengths = load([plibdir plibname '/' plibname '.pulselengths']);
params = load([plibdir plibname '/' plibname '.infofile']);
fs = params(14);
pulsemaxlen = floor(params(11)/1000*fs);

% Reshape pulse data to matrix
P = reshape(pulses,pulsemaxlen,length(pulses)/pulsemaxlen)';
Npulses = size(P,1);

% Make directory for individual pulses
if ~exist([plibdir plibname '/separate_pulses'],'dir')
    system(['mkdir ' plibdir plibname '/separate_pulses']);
end

% Write pulses to files
for i = 1:Npulses
    p = P(i,1:pulselengths(i));
    fid = fopen([plibdir plibname '/separate_pulses/pulse' int2str(i)],'wt');
    fprintf(fid,'%1.6f\n',p);
    fclose(fid);
end

% Interpolate pulses in order to evaluate the mean
pulses_interp = zeros(Npulses,mean_pulse_length);
for i = 1:Npulses
    p = pulses((i-1)*pulsemaxlen+1:(i-1)*pulsemaxlen+pulselengths(i));
    pulses_interp(i,:) = interp1(1:length(p),p,linspace(1,length(p),mean_pulse_length),'cspline');
end

% Find pulse closest to the mean
pulse_mean = mean(pulses_interp,1);            
rmse = zeros(size(pulselengths));
for i = 1:length(pulselengths)
    rmse(i) = sum((pulses_interp(i,:)-pulse_mean).^2);
end
minind = find(rmse == min(rmse),1);
pulse_proto = pulses((minind-1)*pulsemaxlen+1:(minind-1)*pulsemaxlen+pulselengths(minind));

% Save mean pulse
fid = fopen([plibdir plibname '/separate_pulses/pulse_mean'],'wt');
fprintf(fid,'%1.9f\n',pulse_mean);
fclose(fid);

% Save prototype pulse
fid = fopen([plibdir plibname '/separate_pulses/pulse_proto'],'wt');
fprintf(fid,'%1.9f\n',pulse_proto);
fclose(fid);
disp('Done')
