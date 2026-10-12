function zero_pot
%
% Open a Simulink RTW model that drives the Quanser SRV02 potentiometer to
% zero value. Assumes the potentiometer is in input 1.

% Sept2023, J. Gaspar

% ( useMIO16e4 )
% if nargin<1
%     useMIO16e4= 0;
% end
% 
% if ~useMIO16e4
%     % usage: zero_pot
%     open('zero_pci6221.mdl')
% else
%     % usage: zero_pot(1)
%     open('zero_mio16e4.mdl')
% end

NIBoardModel= niboard('get');
mname= ['zero_' NIBoardModel];
fname= [mname '.mdl'];

str= {'Confirm rotate flexbar to set potentiometer to zero value?'; ...
    'Important: S1 cable must be connected so the potentiometer can be read.'};
str= questdlg(str, 'move confirm', 'Yes','No','No');
if strcmp(str, 'No')
    return;
end

open(fname);
set_param(mname, 'SimulationCommand', 'start')
