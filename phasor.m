% -- initialize program
close all; format compact; clc; clear;

% -- set constants
thresh = 0.3;               % threshold [0,1]
harmonic = 2;               % harmonic number (1,2,3,4,...)
singleorbatch = 'batch';   % 'single' or 'batch' processing
lifetimefitting = 'true';   % 'true' for linear fitting of phasor points for lifetime
% 'false' for calculation of lifetime with mean of phasor points

shift_bins = 0;             % number of bins shifted to right from max for minimizing effect of IRF

freq0 = 8E+7;                   % laser frequency, here 80 MHz
delta_t = 1.5625E-10;           % width of one time channel
time_bins = 1 / (freq0 * delta_t) / harmonic;
time_bins = round(time_bins);
% number of time bins for one period (not 64!). 80 for harmonic no 1

% ------ initial calculations ------ %

freq = harmonic * freq0;        % virtual frequency for higher harmonics
w = 2 * pi * freq ;             % angular frequency

% ------ single calculation or batch processing ------ %

switch lower(singleorbatch)
    case {'single'}
        [FileName, PathName] = uigetfile({'*.txt'; '*.asc'}, 'Select FLIM-ASCII-file');
        datadir = [0 0 0] ;     % dummy for for-loop
    case {'batch' }
        PathName = uigetdir('E:\', 'Select FLIM-ASCII-folder');
        datadir = dir(PathName);
        % folder must not contain other files or folders beside FLIM data!
        parts = strsplit(PathName, '\');
        DirPart = parts{end};
end

% ------  ------ %

% ------ begin for-loop for batch data processing ------ %

for jj = 3:numel(datadir) % start with 3, because 1 and 2 are "." and ".." in folder
    
    switch lower(singleorbatch)
        case {'single'}
            decayfile = FileName;
        case {'batch' }
            decayfile = datadir(jj).name;
    end
    
    tic
    
    % ------ load data ------ %
    
    fidd = [PathName '\' decayfile];
    file_dec = fopen(fidd,'rt');
    initials = fscanf(file_dec, ['Number of pixels in line: %d \n Number of lines: %d \n Output sequence: \n Pixel(x=%*d,y=%*d,t=1)...Pixel(x=%*d,y=%*d,t=%d)']);
    fclose(file_dec);
    rawlines = initials(1,1);
    rawcolum = initials(2,1);
    tc_data = initials(3,1);
    
    s = repmat('%d ', 1, tc_data);
    file_dec = fopen(fidd,'rt');
    D = textscan(file_dec, s, 'HeaderLines', 10);
    fclose(file_dec);
    E = [D{:}];
    
    % Find Max and remove data before max
    maxE = sum(E,1);
    [maxdata,I] = max(maxE);
    decdata = E(:,I+shift_bins:end);
    timechannels_data = size(decdata,2);
    decdatatooshort = 0;
    
    if timechannels_data < time_bins
        decdata(1,time_bins) = 0;
        decdatatooshort = 1;
    elseif  timechannels_data > time_bins
        decdata = decdata(:,1:time_bins);
    end
    
    %decdata = int32(medfilt1(double(decdata)));        %smoothing
    
    % remove offset from data
    data_off = mean(E(:,round(I/3):round(2*I/3)),2);
    data_off = repmat(data_off,1,size(decdata,2));
    decdata = decdata - int32(data_off) ;
    decdata(decdata<0) = 0 ;
    
    % Threshold (not peak vs peak. Used whole intensity vs max whole intensity)
    maxmax = max(sum(decdata,2));
    rows_to_remove = any(sum(decdata,2) < (thresh * maxmax), 2);
    decdata(rows_to_remove,:) = [];
    numberofpixels = size(E,1) - sum(rows_to_remove,1);
    
    % calculate G/S-sin/cos-matrix
    tb_vec = linspace(1,time_bins,time_bins)'; % time bin indexing vector
    
    Gn_ma = cos(w * delta_t * (tb_vec - 0.5));
    Sn_ma = sin(w * delta_t * (tb_vec - 0.5));
    
    % calculate data phasor
    Gn = double(decdata) * double(Gn_ma) ;  %take decdata from preprocessing
    Sn = double(decdata) * double(Sn_ma) ;
    area = sum(decdata, 2) ;
    
    % normalization
    G_f = Gn ./ area;
    S_f = Sn ./ area;
    
    % data
    Z = [G_f,S_f];
    
    % ------ Binning for 2D histogram ------ %
    
    % bin centers
    steps = 0.005;
    G_f_bins = 0:steps:1;
    S_f_bins = 0:steps:1;
    G_f_NumBins = numel(G_f_bins);
    S_f_NumBins = numel(S_f_bins);
    
    % map X/Y values to bin indices
    G_f_i = round( interp1(G_f_bins, 1:G_f_NumBins, G_f, 'linear', 'extrap') );
    S_f_i = round( interp1(S_f_bins, 1:S_f_NumBins, S_f, 'linear', 'extrap') );
    
    % limit indices to the range [1,numBins]
    G_f_i = max( min(G_f_i,G_f_NumBins), 1);
    S_f_i = max( min(S_f_i,S_f_NumBins), 1);
    
    % count number of elements in each bin
    H = accumarray([S_f_i(:) G_f_i(:)], 1, [S_f_NumBins G_f_NumBins]);
    
    % plot 2D histogram (w/ contour)
    figure('Name','Phasor Plot');
    contourf(G_f_bins, S_f_bins, H);
    imagesc(G_f_bins, S_f_bins, H)
    j = jet;
    j(1,:) = [ 1 1 1 ];     %comment that out for blue background in plot
    colormap(j);
    cb=colorbar;
    ylabel(cb,'phasor counts')
    set(gca,'Ydir','Normal')
    hold on
    %plot(G_f, S_f, 'b.', 'MarkerSize',1);  % scatter plot of all phasor points
    
    % ------ ellipse ------ %
    
    % Calculate the eigenvectors and eigenvalues
    covariance = cov(Z);
    [eigenvec, eigenval ] = eig(covariance);
    
    % Get the index of the largest eigenvector
    [largest_eigenvec_ind_c, r] = find(eigenval == max(max(eigenval)));
    largest_eigenvec = eigenvec(:, largest_eigenvec_ind_c);
    
    % Get the largest eigenvalue
    largest_eigenval = max(max(eigenval));
    
    % Get the smallest eigenvector and eigenvalue
    if(largest_eigenvec_ind_c == 1)
        smallest_eigenval = max(eigenval(:,2));
        smallest_eigenvec = eigenvec(:,2);
    else
        smallest_eigenval = max(eigenval(:,1));
        smallest_eigenvec = eigenvec(1,:);
    end
    
    % Calculate the angle between x-axis and largest eigenvector
    angle = atan2(largest_eigenvec(2), largest_eigenvec(1));
    
    % This angle is between -pi and pi.
    % shift between 0 and 2pi:
    if(angle < 0)
        angle = angle + 2*pi;
    end
    
    % Get the coordinates of the data mean
    avg = mean(Z);
    
    % Get the 95% confidence interval error ellipse
    chisquare_val = 2.4477;
    theta_grid = linspace(0,2*pi);
    phi = angle;
    X0 = avg(1);
    Y0 = avg(2);
    a = chisquare_val * sqrt(largest_eigenval);
    b = chisquare_val * sqrt(smallest_eigenval);
    
    % the ellipse in x and y coordinates
    ellipse_x_r = a*cos( theta_grid );
    ellipse_y_r = b*sin( theta_grid );
    
    % rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];
    
    % rotate the ellipse to angle phi
    r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;
    
    % Draw the error ellipse
    plot(r_ellipse(:,1) + X0, r_ellipse(:,2) + Y0,'-g','LineWidth',1);
    
    % ------ uni. circle, fit, etc. ------ %
    
    % Plotuniversal circle
    x = 0:0.005:1;
    circle = sqrt(0.25 - (x - 0.5) .^ 2);
    plot(x,circle,'r','LineWidth',1);
    axis equal tight
    axis([0 1 0 0.6]);
    xlabel('G');
    ylabel('S');
    
    % ------ calculate lifetimes ------ %
    
    switch lower(lifetimefitting)
        case {'true'}
            x = 0:0.001:1;
            ell_fit = @(x) (largest_eigenvec(2) / largest_eigenvec(1)) * (x - X0) + Y0;
            circle = @(x) sqrt(0.25 - (x - 0.5) .^ 2);
            
            lincir = @(x) sqrt(0.25 - (x - 0.5) .^ 2) - (largest_eigenvec(2) / largest_eigenvec(1)) * (x - X0) - Y0;
            [maxdummy,X_max_lincir] = max(lincir(x));
            X_inter1 = fzero(lincir,[0 x(X_max_lincir)]);
            X_inter2 = fzero(lincir,[x(X_max_lincir) 1]);
            
            Y_inter1 = ell_fit(X_inter1);
            Y_inter2 = ell_fit(X_inter2);
            
            tau_1 = 1 / w * (Y_inter2 / X_inter2) * 1E12;       % lifetimes in ps
            tau_2 = 1 / w * (Y_inter1 / X_inter1) * 1E12;       % lifetimes in ps
            
            plot(X_inter1,circle(X_inter1),'ro','LineWidth',1);
            plot(X_inter2,circle(X_inter2),'ro','LineWidth',1);
            
            plot(x,ell_fit(x),'c','LineWidth',1)
            
            % ------ KOS-Trafo for Histogram ------- %
            
            % Translation
            translationmatrix = [X_inter1,Y_inter1];
            Z_tr = Z - repmat(translationmatrix,size(Z,1),1);
            
            % Rotation
            alpha = atan((Y_inter1 - Y_inter2)/(X_inter1 - X_inter2));
            rotmat = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
            Z_rot = Z_tr * rotmat;
            
            % Scaling
            scale = sqrt((X_inter1 - X_inter2) ^ 2 + (Y_inter1 - Y_inter2) ^ 2);
            Z_scale = Z_rot ./ scale;
            
            % Histogram
            figure('Name','Histogram')
            h = histogram(Z_scale(:,1),'Normalization','pdf');
            h.BinLimits = [0:1];
            
            % Values Histogram
            alpha_1 = mean(Z_scale(:,1),1);
            alpha_1_std = std(Z_scale(:,1),0,1);
            
        case {'false' }
            tau_1 = 1 / w * (Y0 / X0) * 1E12;
            tau_2 = 0;
            alpha_1 = 1;
            alpha_1_std = 0;
    end
    
    % ------ save temporarily variables to array ------ %
    
    kk = jj-2;
    decayfile_exp{kk} = decayfile;
    tau_1_exp(kk) = tau_1;
    tau_2_exp(kk) = tau_2;
    alpha_1_exp(kk) = alpha_1;
    alpha_1_std_exp(kk) = alpha_1_std;
    a_exp(kk) = a;
    b_exp(kk) = b;
    atob(kk) = a/b;
    twoexpconf(kk) = 1-exp((1-atob(kk))/1.5) ; %1.5 is  a random factor, can bechanged
    thresh_exp(kk) = thresh;
    harmonickk(kk) = harmonic;
    rawlines_exp(kk) = rawlines;
    rawcolum_exp(kk) = rawcolum;
    tc_data_exp(kk) = timechannels_data;
    time_bins_needed(kk) = time_bins;
    if decdatatooshort == 1;
        toolesstc(kk) = {'Too less time channels. Try again with higher harmonic'};
    else toolesstc(kk) = {'fine'};
    end
    exec_time1_exp(kk) = toc;
    pixelnumber(kk) = numberofpixels;
    
    counter = [num2str(jj-2),' of ',num2str(numel(datadir)-2),' done'];
    disp(counter)
    
    switch lower(singleorbatch)
        case {'batch'}
            close all
    end
    
end

% ------ collect and export results ------ %

results_phasor = table(decayfile_exp', tau_1_exp', tau_2_exp', ...
    alpha_1_exp', alpha_1_std_exp', a_exp', b_exp', atob', twoexpconf', ...
    thresh_exp', harmonickk', exec_time1_exp', rawlines_exp', ...
    rawcolum_exp', pixelnumber', tc_data_exp', time_bins_needed', toolesstc',...
    'VariableNames', {'filename' 'tau_1' 'tau_2' 'alpha_1' 'alpha_1_std' ...
    'a_ellipse' 'b_ellipse' 'a_b_ratio' 'twoexp_confidence' 'threshold' ...
    'harmonic' 'exec_time1' 'pixel_lines' 'pixel_columns' 'pixelnumber' ...
    'time_channels_data' 'time_bins_needed' 'no_of_tc'})

% ------ write output ------ %

switch lower(singleorbatch)
    case {'single'}
        parts = strsplit(FileName, '.');
        FileName = parts{1};
        Filenameout_sug = [PathName '\' FileName '_phasorresult.txt'];
    case {'batch' }
        lengthofpath = numel(PathName) - numel(DirPart);
        Filenameout_sug = [PathName(1:lengthofpath) '\' DirPart '_phasorresult.txt'];
end

[FileNameout, PathNameout] = uiputfile(Filenameout_sug, 'Save phasor results as...' );
outputfile = [PathNameout '\' FileNameout];
writetable(results_phasor,outputfile,'Delimiter','\t');

