function  Control_GUI
  
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
hfi = findall(0,'Name',sprintf('UI for faults and damages at %s.mdl',modelName));

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
    hf = figure('Toolbar','none',...
        'MenuBar','none',...
        'IntegerHandle','off',...
        'Units','normalized',...
        'Resize','on',...
        'NumberTitle','off',...
        'HandleVisibility','callback',...
        'Name',sprintf('UI for faults and damages at %s.mdl',modelName),...
        'CloseRequestFcn',@localCloseRequestFcn,...
        'Visible','off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Main Panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MDL = 'F16ASYM_Controlled';

    global oversize_param
    oversize_param=0;
    Main = uipanel('Parent',hf,'Title','Main Panel','FontSize',20,...
        'BackgroundColor','white','Units','normalized',...
        'Position',[0.02 -oversize_param 0.95 1+oversize_param]);
    
    uicontrol('Style','Slider','Parent',hf,...              % Slider
      'Units','normalized','Position',[0.97 0 0.03 1],...
      'Value',1,'Callback',{@slider_callback1,Main});

        uicontrol('Parent',Main,...
            'Style','text','FontSize',16,...
            'Units','normalized',...
            'Position',[0 1-0.06 0.95 0.05],...
            'String','This GUI allows the user to inputs faults, changes in control surfaces and airframe, as well as, control the learning proces of ANNs',...
            'Backgroundcolor',[1 1 1]);


      
    %% Control surfaces 
    panel_control = uipanel('Parent',Main,'Title','Control surfaces','FontSize',12, 'Units','normalized','Position',[.02 .25 .45 .7]);
        % Blockades
        panel_control_blockades = uipanel('Parent',panel_control,'Title','Blockades','FontSize',10, 'Units','normalized','Position',[.05 .5 .9 .49]);  
        string={'Elev_L','Elev_R','Ail_L','Ail_R','Rudd','LEF_L','LEF_R'};

                panel_control_blockades_at = uipanel('Parent',panel_control_blockades,'Title','Blockade at:','FontSize',10, 'Units','normalized','Position',[.01 .4 .98 .55]);  
                    % Individual
                    for i=1:7
                            uicontrol('Parent',panel_control_blockades_at,'Style','text','FontSize',10,'BackgroundColor','white','String',string{i},'Units','normalized','Position',[0.02+0.95/7*(i-1),0.8,0.95/7,0.15]);
                            Val_h{i} = uicontrol('Parent',panel_control_blockades_at,'Style','edit','FontSize',10,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/7*(i-1),0.5,0.95/7,0.3]); 
                            uicontrol('Parent',panel_control_blockades_at,'Style','togglebutton','Tag',['Block_at',num2str(i)],'String','Block','Units','normalized','Position',[0.02+0.95/7*(i-1),0.1,0.95/7,0.3],...
                                'Callback',{@Blockat_CurrentValue_Callback,Val_h{i},i}) ;
                    end
                panel_control_blockades_now = uipanel('Parent',panel_control_blockades,'Title','Blockade NOW','FontSize',10,'Units','normalized','Position',[.01 .02 .98 .36]);
                    % Individual
                    for i=1:7
                            uicontrol('Parent',panel_control_blockades_now,'Style','text','FontSize',10,'BackgroundColor','white','String',string{i},'Units','normalized','Position',[0.02+0.95/7*(i-1),0.8,0.95/7,0.2]);
                            uicontrol('Parent',panel_control_blockades_now,'Style','togglebutton','Tag',['Block_now',num2str(i)],'String','Block now','Units','normalized','Position',[0.02+0.95/7*(i-1),0.1,0.95/7,0.6],...
                                'Callback',{@Blocknow_CurrentValue_Callback,i}) ;
                    end
        % Floating or loss            
        panel_control_Float = uipanel('Parent',panel_control,'Title','Floating or loss','FontSize',10,'Units','normalized','Position',[.05 .30 .9 .17]);  
            for i=1:7
                uicontrol('Parent',panel_control_Float,'Style','text','FontSize',10,'BackgroundColor','white','String',string{i},'Units','normalized','Position',[0.02+0.95/7*(i-1),0.8,0.95/7,0.2]);
                uicontrol('Parent',panel_control_Float,'Style','togglebutton','Tag',['Float',num2str(i)],'String','Lose or float','Units','normalized','Position',[0.02+0.95/7*(i-1),0.1,0.95/7,0.6],...
                    'Callback',{@Float_loss_Callback,i}) ;
            end
        % Lose effectivenes or surface    
        panel_control_Area = uipanel('Parent',panel_control,'Title','Loss effectiveness or Surface: Input % lost','FontSize',10,'Units','normalized','Position',[.05 .02 .9 .27]);  
            for i=1:7
                uicontrol('Parent',panel_control_Area,'Style','text','FontSize',10,'BackgroundColor','white','String',string{i},'Units','normalized','Position',[0.02+0.95/7*(i-1),0.8,0.95/7-0.02,0.15]);
                Val_h{i} = uicontrol('Parent',panel_control_Area,'Style','edit','FontSize',10,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/7*(i-1),0.5,0.95/7-0.03,0.3]);
                uicontrol('Parent',panel_control_Area,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/7*(i)-0.03,0.4,0.03,0.3]);
                uicontrol('Parent',panel_control_Area,'Style','togglebutton','Tag',['Block_at',num2str(i)],'String','Affect','Units','normalized','Position',[0.02+0.95/7*(i-1),0.1,0.95/7,0.3],...
                    'Callback',{@Loss_area_Callback,Val_h{i},i}) ;
            end
    
    %% Structural changes
   panel_control2 = uipanel('Parent',Main,'Title','Airframe and global changes','FontSize',12,'Units','normalized','Position',[.5 .25 .48 .7]);

        panel_control_massprop = uipanel('Parent',panel_control2,'Title','Mass properties changes','FontSize',10,...
            'Units','normalized','Position',[.05 .57 .9 .4]); 
                Val_h={};
                % mass and xcg
                i=1;
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta mass(%) of 9295.4kg','Units','normalized','Position',[0.02+0.95/2*(i-1),0.8,0.95/2-0.03,0.15]);
                Val_h{1}=uicontrol('Parent',panel_control_massprop,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/2*(i-1),0.5,0.95/2-0.03,0.3]);
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/2*(i)-0.03,0.55,0.03,0.15]);
                i=2;
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'BackgroundColor','white','String','x_cg position in% of CMA(0.3%)','Units','normalized','Position',[0.02+0.95/2*(i-1),0.8,0.95/2-0.03,0.15]);
                Val_h{2}=uicontrol('Parent',panel_control_massprop,'Style','edit','FontSize',12,'BackgroundColor','white','String','0.3','Units','normalized','Position',[0.02+0.95/2*(i-1),0.5,0.95/2-0.03,0.3]);
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/2*(i)-0.03,0.55,0.03,0.15]);
                % Moments of inertia
                i=1;
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta J_y(%) of 75673.6kg.m^2','Units','normalized','Position',[0.02+0.95/4*(i-1),0.3,0.95/4-0.03,0.2]);
                Val_h{3}=uicontrol('Parent',panel_control_massprop,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/4*(i-1),0.0,0.95/4-0.03,0.3]);
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/4*(i)-0.03,0.05,0.03,0.15]);
                i=2;
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta J_xz(%) of 1331.4kg.m^2','Units','normalized','Position',[0.02+0.95/4*(i-1),0.3,0.95/4-0.03,0.2]);
                Val_h{4}=uicontrol('Parent',panel_control_massprop,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/4*(i-1),0.0,0.95/4-0.03,0.3]);
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',12,'String','%','Units','normalized','Position',[0.02+0.95/4*(i)-0.03,0.05,0.03,0.15]);
                i=3;
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta J_z(%) of 85552.1kg.m^2','Units','normalized','Position',[0.02+0.95/4*(i-1),0.3,0.95/4-0.03,0.2]);
                Val_h{5}=uicontrol('Parent',panel_control_massprop,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/4*(i-1),0.0,0.95/4-0.03,0.3]);
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/4*(i)-0.03,0.05,0.03,0.15]);
                i=4;
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta J_x(%) of 12874.8kg.m^2','Units','normalized','Position',[0.02+0.95/4*(i-1),0.3,0.95/4-0.03,0.2]);
                Val_h{6}=uicontrol('Parent',panel_control_massprop,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/4*(i-1),0.0,0.95/4-0.03,0.3]);
                uicontrol('Parent',panel_control_massprop,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/4*(i)-0.03,0.05,0.03,0.15]);
                
        panel_control_AeroAndArea = uipanel('Parent',panel_control2,'Title','Aerodynamics and surfaces changes','FontSize',10,...
            'Units','normalized','Position',[.05 .17 .9 .4]);  
                % Surface changes
                i=1;
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'BackgroundColor','white','String','Left wing surface loss(%) of  11.14m^2','Units','normalized','Position',[0.02+0.95/3*(i-1),0.75,0.95/3-0.03,0.2]);
                Val_h{7}=uicontrol('Parent',panel_control_AeroAndArea,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/3*(i-1),0.5,0.95/3-0.03,0.3]);
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/3*(i)-0.03,0.55,0.03,0.15]);
                i=2;
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'BackgroundColor','white','String','Right wing surface loss(%) of  11.14m^2','Units','normalized','Position',[0.02+0.95/3*(i-1),0.75,0.95/3-0.03,0.2]);
                Val_h{8}=uicontrol('Parent',panel_control_AeroAndArea,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/3*(i-1),0.5,0.95/3-0.03,0.3]);
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',12,'String','%','Units','normalized','Position',[0.02+0.95/3*(i)-0.03,0.55,0.03,0.15]);
                i=3;
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'BackgroundColor','white','String','Fin surface loss(%) of  6.56m^2','Units','normalized','Position',[0.02+0.95/3*(i-1),0.75,0.95/3-0.03,0.2]);
                Val_h{9}=uicontrol('Parent',panel_control_AeroAndArea,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/3*(i-1),0.5,0.95/3-0.03,0.3]);
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'String','%','Units','normalized','Position',[0.02+0.95/3*(i)-0.03,0.55,0.03,0.15]);
                % DElta coeffs
                i=1;
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta C_Lift (additionally to ~0.2)','Units','normalized','Position',[0.02+0.95/3*(i-1),0.3,0.95/3-0.03,0.2]);
                Val_h{10}=uicontrol('Parent',panel_control_AeroAndArea,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/3*(i-1),0.0,0.95/3-0.03,0.3]);
                i=2;
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta C_Drag (additionally to ~0.036)','Units','normalized','Position',[0.02+0.95/3*(i-1),0.3,0.95/3-0.03,0.2]);
                Val_h{11}=uicontrol('Parent',panel_control_AeroAndArea,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/3*(i-1),0.0,0.95/3-0.03,0.3]);
                i=3;
                uicontrol('Parent',panel_control_AeroAndArea,'Style','text','FontSize',10,'BackgroundColor','white','String','Delta C_pitch_moment (additionally to ~ - 0.054)','Units','normalized','Position',[0.02+0.95/3*(i-1),0.3,0.95/3-0.03,0.2]);
                Val_h{12}=uicontrol('Parent',panel_control_AeroAndArea,'Style','edit','FontSize',12,'BackgroundColor','white','String','0','Units','normalized','Position',[0.02+0.95/3*(i-1),0.0,0.95/3-0.03,0.3]);

        
        
         uicontrol('Parent',panel_control2,'String','Update Changes','FontSize',16,'Units','normalized',...
            'Position',[0.65 0.02 0.3 0.13],'Callback',{@Surf_and_Mass_changes,Val_h});
         %% ANN control
   panel_controlANN = uipanel('Parent',Main,'Title','ANNs control','FontSize',12,'Units','normalized','Position',[.02 .02 .96 .22]);
            % Roll
           panel_control_Roll = uipanel('Parent',panel_controlANN,'Title','Roll chanel','FontSize',12,...
                'Units','normalized','Position',[.01 .02 .98/4-0.01 0.95]);
                % Learning rate
                uicontrol('Parent',panel_control_Roll,'Style','text','FontSize',11,'BackgroundColor','white','String','Learning Rate (~20)','Units','normalized',...
                        'Position',[0.02 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'Roll_learning rate','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Roll,'Tag','Roll_learn','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.28 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Roll,'FontSize',10,'String','Speed-up','Units','normalized',...
                        'Position',[0.02 0.43 0.25 0.35],'Callback',{@Change_params,{'Roll_learn','up'}});
                uicontrol('Parent',panel_control_Roll,'FontSize',10,'String','Slow-down','Units','normalized',...
                        'Position',[0.02 0.05 0.25 0.35],'Callback',{@Change_params,{'Roll_learn','down'}});  
                % Regularization param
                uicontrol('Parent',panel_control_Roll,'Style','text','FontSize',11,'BackgroundColor','white','String','Reg. Param (~1)','Units','normalized',...
                        'Position',[0.52 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'Roll_reg param','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Roll,'Tag','Roll_reg','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.78 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Roll,'FontSize',10,'String','Underfit','Units','normalized',...
                        'Position',[0.52 0.43 0.25 0.35],'Callback',{@Change_params,{'Roll_reg','up'}});
                uicontrol('Parent',panel_control_Roll,'FontSize',10,'String','Overfit','Units','normalized',...
                        'Position',[0.52 0.05 0.25 0.35],'Callback',{@Change_params,{'Roll_reg','down'}});  
   
            % Long chanel
           panel_control_Long  = uipanel('Parent',panel_controlANN,'Title','Long chanel','FontSize',12,...
                'Units','normalized','Position',[0.98/4+0.005 .02 .98/4-0.005 0.95]); 
                % Learning rate
                uicontrol('Parent',panel_control_Long,'Style','text','FontSize',11,'BackgroundColor','white','String','Learning Rate (~20)','Units','normalized',...
                        'Position',[0.02 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'Long_learning rate','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Long,'Tag','Long_learn','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.28 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Long,'FontSize',10,'String','Speed-up','Units','normalized',...
                        'Position',[0.02 0.43 0.25 0.35],'Callback',{@Change_params,{'Long_learn','up'}});
                uicontrol('Parent',panel_control_Long,'FontSize',10,'String','Slow-down','Units','normalized',...
                        'Position',[0.02 0.05 0.25 0.35],'Callback',{@Change_params,{'Long_learn','down'}});  
                % Regularization param
                uicontrol('Parent',panel_control_Long,'Style','text','FontSize',11,'BackgroundColor','white','String','Reg. Param (~1)','Units','normalized',...
                        'Position',[0.52 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'Long_reg param','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Long,'Tag','Long_reg','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.78 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Long,'FontSize',10,'String','Underfit','Units','normalized',...
                        'Position',[0.52 0.43 0.25 0.35],'Callback',{@Change_params,{'Long_reg','up'}});
                uicontrol('Parent',panel_control_Long,'FontSize',10,'String','Overfit','Units','normalized',...
                        'Position',[0.52 0.05 0.25 0.35],'Callback',{@Change_params,{'Long_reg','down'}});  
            % Yaw chanel
           panel_control_Yaw  = uipanel('Parent',panel_controlANN,'Title','Yaw chanel chanel','FontSize',12,...
                'Units','normalized','Position',[2*.98/4+0.005 .02 .98/4 0.95]); 
                % Learning rate
                uicontrol('Parent',panel_control_Yaw,'Style','text','FontSize',11,'BackgroundColor','white','String','Learning Rate (~20)','Units','normalized',...
                        'Position',[0.02 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'Yaw_learning rate','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Yaw,'Tag','Yaw_learn','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.28 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Yaw,'FontSize',10,'String','Speed-up','Units','normalized',...
                        'Position',[0.02 0.43 0.25 0.35],'Callback',{@Change_params,{'Yaw_learn','up'}});
                uicontrol('Parent',panel_control_Yaw,'FontSize',10,'String','Slow-down','Units','normalized',...
                        'Position',[0.02 0.05 0.25 0.35],'Callback',{@Change_params,{'Yaw_learn','down'}});  
                % Regularization param
                uicontrol('Parent',panel_control_Yaw,'Style','text','FontSize',11,'BackgroundColor','white','String','Reg. Param (~1)','Units','normalized',...
                        'Position',[0.52 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'Yaw_reg param','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Yaw,'Tag','Yaw_reg','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.78 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Yaw,'FontSize',10,'String','Underfit','Units','normalized',...
                        'Position',[0.52 0.43 0.25 0.35],'Callback',{@Change_params,{'Yaw_reg','up'}});
                uicontrol('Parent',panel_control_Yaw,'FontSize',10,'String','Overfit','Units','normalized',...
                        'Position',[0.52 0.05 0.25 0.35],'Callback',{@Change_params,{'Yaw_reg','down'}});  
                    
            % Complete ANN
           panel_control_Com  = uipanel('Parent',panel_controlANN,'Title','Complete ANN','FontSize',12,...
                'Units','normalized','Position',[3*.98/4+0.01 .02 .98/4 0.95]); 
                % Learning rate
                uicontrol('Parent',panel_control_Com,'Style','text','FontSize',11,'BackgroundColor','white','String','Learning Rate (~20)','Units','normalized',...
                        'Position',[0.02 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'C learning rate','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Com,'Tag','C_learn','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.28 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Com,'FontSize',10,'String','Speed-up','Units','normalized',...
                        'Position',[0.02 0.43 0.25 0.35],'Callback',{@Change_params,{'C_learn','up'}});
                uicontrol('Parent',panel_control_Com,'FontSize',10,'String','Slow-down','Units','normalized',...
                        'Position',[0.02 0.05 0.25 0.35],'Callback',{@Change_params,{'C_learn','down'}});  
                % Regularization param
                uicontrol('Parent',panel_control_Com,'Style','text','FontSize',11,'BackgroundColor','white','String','Reg. Param (~1)','Units','normalized',...
                        'Position',[0.52 .82 .45 0.22]);
                            Block_search = find_system(MDL, 'Name', 'C reg param','Blocktype','Constant');
                            Actual_val = get_param(Block_search,'Value');
                uicontrol('Parent',panel_control_Com,'Tag','C_reg','Style','text','FontSize',14,'BackgroundColor','white','String',Actual_val,'Units','normalized',...
                        'Position',[0.78 .35 .15 0.25]);
                    
                uicontrol('Parent',panel_control_Com,'FontSize',10,'String','Underfit','Units','normalized',...
                        'Position',[0.52 0.43 0.25 0.35],'Callback',{@Change_params,{'C_reg','up'}});
                uicontrol('Parent',panel_control_Com,'FontSize',10,'String','Overfit','Units','normalized',...
                        'Position',[0.52 0.05 0.25 0.35],'Callback',{@Change_params,{'C_reg','down'}});  
                                        
                    
function Change_params(hObject, eventdata,handles)
MDL = 'F16ASYM_Controlled';

Str_to_change= findall(get(hObject,'parent'),'Tag',handles{1});
switch handles{1}
% Roll    
    case {'Roll_learn'}
        Block_search = find_system(MDL, 'Name', 'Roll_learning rate','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);
    case{'Roll_reg'}
        Block_search = find_system(MDL, 'Name', 'Roll_reg param','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);
 % Long       
    case {'Long_learn'}
        Block_search = find_system(MDL, 'Name', 'Long_learning rate','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);
    case{'Long_reg'}
        Block_search = find_system(MDL, 'Name', 'Long_reg param','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);
        
  % Yaw       
    case {'Yaw_learn'}
        Block_search = find_system(MDL, 'Name', 'Yaw_learning rate','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);
    case{'Yaw_reg'}
        Block_search = find_system(MDL, 'Name', 'Yaw_reg param','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);       

  % Completed       
    case {'C_learn'}
        Block_search = find_system(MDL, 'Name', 'C learning rate','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);
    case{'C_reg'}
        Block_search = find_system(MDL, 'Name', 'C reg param','Blocktype','Constant');
        Actual_val = get_param(Block_search,'Value');
        Actual_val = eval(char(Actual_val));
        switch handles{2}
            case {'up'}
                Actual_val = Actual_val+0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
                
            case {'down'}
                Actual_val = Actual_val-0.1;
                set_param(char(Block_search),'Value',num2str(Actual_val));
        end
        Actual_val = get_param(Block_search,'Value');
        set(Str_to_change,'String',Actual_val);             
end



      
function Surf_and_Mass_changes(hObject, eventdata,handles)

% Data from mdl
MDL = 'F16ASYM_Controlled';
str= {'delta_mass','xcg','Delta_J','Delta_J','Delta_J','Delta_J','delta_S_L','delta_S_R','delta_S_fin','delta_coef','delta_coef','delta_coef'};
% 'delta_mass','xcg'
for i=1:2
    Block_search = find_system([MDL,'/Faults injection'], 'Name', str{i} );
    Actual_val = get_param(Block_search,'Value');
    Actual_val = eval(char(Actual_val));

    switch str{i}
        case {'delta_mass'}
            NewStrVal = get(handles{i}, 'String');
            Actual_val =  str2double(NewStrVal)/100;
        case {'xcg'}
            NewStrVal = get(handles{i}, 'String');
            Actual_val =  str2double(NewStrVal);
    end
    set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);
end

% Delta_J
Block_search = find_system([MDL,'/Faults injection'], 'Name', str{5} );
Actual_val = get_param(Block_search,'Value');
Actual_val = eval(char(Actual_val));
for i=3:6

    NewStrVal = get(handles{i}, 'String');
    Actual_val(i-2) =  str2double(NewStrVal)/100;
    
end
set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);
% 'delta_S_L','delta_S_R','delta_S_fin'
for i=7:9
    Block_search = find_system([MDL,'/Faults injection'], 'Name', str{i} );
    Actual_val = get_param(Block_search,'Value');
    Actual_val = eval(char(Actual_val));

    NewStrVal = get(handles{i}, 'String');
    Actual_val =  str2double(NewStrVal)/100;

    set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);
end
% delta_coef
Actual_val=[];
Block_search = find_system([MDL,'/Faults injection'], 'Name', str{10} );
Actual_val = get_param(Block_search,'Value');
Actual_val = eval(char(Actual_val));
for i=10:12

    NewStrVal = get(handles{i}, 'String');
    Actual_val(i-9) =  str2double(NewStrVal);

end

set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);

function Loss_area_Callback(hObject, eventdata,handles, num)
% Data from mdl
MDL = 'F16ASYM_Controlled';
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'effectiveness');
Actual_val = get_param(Block_search,'Value');
Actual_val = eval(char(Actual_val));

% Disappear the other blocking option
hf = get(hObject,'parent');
hf = get(hf,'parent');hf = get(hf,'parent');hf = get(hf,'parent');
Elem = findall(hf,'Tag',['Float',num2str(num)]);


on=get(hObject, 'Value');
if on
    NewStrVal = get(handles, 'String'); 
    Actual_val(1,num) =  1-str2double(NewStrVal)/100;
    set(hObject,'String','Affected');
    set(Elem,'Visible','off')
else
    Actual_val(1,num) =  1;
    set(hObject,'String','Affect');
    set(Elem,'Visible','on')

end
set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);

function Float_loss_Callback(hObject, eventdata, num)
% Data from mdl
MDL = 'F16ASYM_Controlled';
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'effectiveness');
Actual_val = get_param(Block_search,'Value');
Actual_val = eval(char(Actual_val));

% Disappear the other blocking option
hf = get(hObject,'parent');hf = get(hf,'parent');

Elem = findall(hf,'Tag',['Block_now',num2str(num)]);
Elem2 = findall(hf,'Tag',['Block_at',num2str(num)]);


on=get(hObject, 'Value');
if on
    Actual_val(1,num) =  0;
    set(Elem,'Visible','off')
    set(Elem2,'Visible','off')
else
    Actual_val(1,num) =  1;
    set(Elem,'Visible','on')
    set(Elem2,'Visible','on')
end
set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);

function Blockat_CurrentValue_Callback(hObject, eventdata, handles,num)
% Data from mdl
MDL = 'F16ASYM_Controlled';
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'Blockades at');
Actual_val = get_param(Block_search,'Value');
Actual_val = eval(char(Actual_val));

% Disappear the other blocking option
hf = get(hObject,'parent');
hf = get(hf,'parent');hf = get(hf,'parent');hf = get(hf,'parent');
Elem = findall(hf,'Tag',['Block_now',num2str(num)]);
Elem2 = findall(hf,'Tag',['Float',num2str(num)]);


on=get(hObject, 'Value');
if on
    NewStrVal = get(handles, 'String'); 
    Actual_val(1,num) =  str2double(NewStrVal);
    set(hObject,'String','Blocked');
    set(Elem,'Visible','off')
    set(Elem2,'Visible','off')
else
    Actual_val(1,num) =  NaN;
    set(hObject,'String','Block');
    set(Elem,'Visible','on')
    set(Elem2,'Visible','on')

end
set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);

function Blocknow_CurrentValue_Callback(hObject, eventdata, num)
% Data from mdl

MDL = 'F16ASYM_Controlled';
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'Blockades now');
Actual_val = get_param(Block_search,'Value');
Actual_val = eval(char(Actual_val));


% Disappear the other blocking option
hf = get(hObject,'parent');
hf = get(hf,'parent');hf = get(hf,'parent');hf = get(hf,'parent');
Elem = findall(hf,'Tag',['Block_at',num2str(num)]);
Elem2 = findall(hf,'Tag',['Float',num2str(num)]);

on=get(hObject, 'Value');
if on
    Actual_val(1,num) = 1;
    set(hObject,'String','Blocked');
    set(Elem,'Visible','off')
    set(Elem2,'Visible','off')

else
    Actual_val(1,num) =  0;
    set(hObject,'String','Block now');
    set(Elem,'Visible','on')
    set(Elem2,'Visible','on')

end

set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);

function modelExists = localValidateInputs(modelName)

num = exist(modelName,'file');
if num == 4
    modelExists = true;
else
    modelExists = false;
end

function slider_callback1(hObject, eventdata, handles)
    global oversize_param

val = get(hObject,'Value');
set(handles,'Position',[0.02 -oversize_param*val 0.95 1+oversize_param])

function localCloseRequestFcn(hObject,eventdata,ad) %#ok
MDL = 'F16ASYM_Controlled';

% Reseting effectiveness
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'effectiveness');
set_param(char(Block_search),'Value',['[',num2str(ones(1,7)),']']);
% Blockades at
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'Blockades at');
set_param(char(Block_search),'Value',['[',num2str(ones(1,7)*NaN),']']);
% Blockades now
Block_search = find_system([MDL,'/Faults injection'], 'Name', 'Blockades now');
set_param(char(Block_search),'Value',['[',num2str(zeros(1,7)),']']);


% Airframe parameters
str= {'delta_mass','xcg','Delta_J','Delta_J','Delta_J','Delta_J','delta_S_L','delta_S_R','delta_S_fin','delta_coef','delta_coef','delta_coef'};
% 'delta_mass','xcg'
Actual_val=[];
for i=1:2
    Block_search = find_system([MDL,'/Faults injection'], 'Name', str{i} );
    switch str{i}
        case {'delta_mass'}
            Actual_val =  0;
        case {'xcg'}
            Actual_val =  0.3;
    end
    set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);
end
% 'delta_S_L','delta_S_R','delta_S_fin'
Actual_val=[];
for i=7:9
    Block_search = find_system([MDL,'/Faults injection'], 'Name', str{i} );

    Actual_val =  0;

    set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);
end
% Delta_J
Actual_val=[];
Block_search = find_system([MDL,'/Faults injection'], 'Name', str{5} );
for i=3:6
    Actual_val(i-2) =  0;
    
end
set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);

% delta_coef
Actual_val=[];
Block_search = find_system([MDL,'/Faults injection'], 'Name', str{10} );
for i=10:12
    Actual_val(i-9) = 0;
end
set_param(char(Block_search),'Value',['[',num2str(Actual_val),']']);




delete(hObject)
% close all Force