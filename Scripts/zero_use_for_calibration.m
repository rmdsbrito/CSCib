function [offset, gain]= zero_use_for_calibration
%
% Strain Gauge - compute offset correction to apply before the gain K_GAGE
%
% Let x=ADC, x0=strain_gauge_offset, k0=strain_gauge_gain_after_offset
% g= K*k0*(x-x0)
% x0= g/K/k0-g

% April 2016, Sept 2023 (new Simulink files), JG

flg1= evalin('base',['exist(''K_GAGE'',''var'')']);
flg2= evalin('base',['exist(''zero_g'',''var'')']);
if flg1 && flg2
    %error('missing vars, please run "setup_ceelab.m" and "zero_*.mdl"')
    % ^assume init work was well
    [offset, gain]= calc_data_flexbar_old;
else
    [offset, gain]= calc_data_flexbar;
end

if nargout<1
    assignin('base', 'strain_gauge_offset', offset);
    assignin('base', 'strain_gauge_gain_after_offset', gain);
    fprintf(1, 'Strain gauge, offset=%f, gain=%f\n', offset, gain);
end


function [offset, gain]= calc_data_flexbar_old
g= evalin('base', 'median(zero_g(:,2))');
K= evalin('base', 'K_GAGE');
offset= g/K/1.0-g;
gain= 1;  % do nothing


function [offset, gain]= calc_data_flexbar
% ret= anim_arm('getSignals'); % get ScopeData or data in a file
% offset= median(ret.g);
ScopeData= view_io('getSignals');
[~,~,~,g]= get_data_from_scope(ScopeData);
offset= median(g);
gain= 1;  % do nothing


function [t,r,p,g]= get_data_from_scope(ScopeData)

% TO DO: accept students creating struct with time and just 3 columns of
% data, namely [ref pot straing], or even just 2 columns [pot strain]

% ScopeData is an array of values
if ~isstruct(ScopeData)
    error_if_lesser_4columns( ScopeData );
    t= ScopeData(:,1);
    r= ScopeData(:,2);
    p= ScopeData(:,3);
    g= ScopeData(:,4);
    return
end

% ScopeData is a struct
error_if_lesser_4columns( ScopeData.signals.values );
t= ScopeData.signals.values(:,1);
r= ScopeData.signals.values(:,2);
p= ScopeData.signals.values(:,3);
g= ScopeData.signals.values(:,4);

return


function error_if_lesser_4columns( ScopeData )
if size(ScopeData, 2) < 4
    % in the zero_*.mdl time is included in the signals list
    error('expecting time added into the signals, making an Nx4 array [time ref pot strain]')
end
