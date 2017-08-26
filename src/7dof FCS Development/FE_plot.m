function  FE_plot
  
modelName = 'F16ASYM_Controlled';
% Do some simple error checking on the input
if ~localValidateInputs(modelName)
    estr = sprintf('The model %s.mdl cannot be found.',modelName);
    errordlg(estr,'Model not found error','modal');
    return
end
% Do some simple error checking on varargout
error(nargoutchk(0,1,nargout));

% Create the UI if one does not already exist.
% Bring the UI to the front if one does already exist.
hfi = findall(0,'Name',sprintf('Plot of Flight envelope position at %s.mdl',modelName));

if isempty(hfi)
    % Create a UI
    hfi = localCreateUI(modelName);
    figure(hfi);
else
    % Bring it to the front
    figure(hfi);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create the user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hf = localCreateUI(modelName)

    % Create the figure, setting appropriate properties
    hf = figure('IntegerHandle','off',...
        'Units','normalized',...
        'Resize','on',...
        'NumberTitle','off',...
        'HandleVisibility','callback',...
        'Name',sprintf('Plot of Flight envelope position at %s.mdl',modelName),...
        'CloseRequestFcn',@localCloseRequestFcn,...
        'Visible','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Main Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MDL = 'F16ASYM_Controlled';

    Main = uipanel('Parent',hf,'Title','Flight Envelope plot','FontSize',15,...
        'BackgroundColor','white','Units','normalized',...
        'Position',[0.02 0.02 0.95 0.95]);

    % Plots and axis
    global hplot htext hlist

    AX=axes('Parent',Main,'Units','normalized ','Position',[0.1 0.1 0.85 0.85]);
    load('FlightEnvelope.mat')

    plot(AX,F_envelope.M_1g ,F_envelope.H_1g ,'r')
    xlabel(AX,'Mach')
    ylabel(AX,'Alt (m)')
    % PLot th epoint
    hold(AX)
    hplot=scatter(AX,0,0,'o','filled');
    hChildren = get(hplot, 'Children');
    set(hChildren, 'Markersize', 10)
    htext = text(0,0,['    YOU   ';'\downarrow'],'HorizontalAlignment','Center','VerticalAlignment','Bottom','FontSize',18,'Parent', AX);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Add Listener
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Block_search = [MDL,'/Environment params Estimations/ALT_and_M_GAIN_LIST'];

hlist = add_exec_event_listener(Block_search, 'PostOutputs', @localEventListener);  




function modelExists = localValidateInputs(modelName)

num = exist(modelName,'file');
if num == 4
    modelExists = true;
else
    modelExists = false;
end

    
  
function localCloseRequestFcn(hObject,eventdata,ad) %#ok
    global hplot hlist

    try
        delete(hlist)
        delete(hObject)
    catch
        close all Force
    end