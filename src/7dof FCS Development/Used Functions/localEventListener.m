function localEventListener(block, eventdata) 
global hplot htext hlist

% % fprintf('\n Event has occured!')

% Gets the time and output value
Data = block.OutputPort(1).Data;



% Update point coordinates
set(hplot,'XData',Data(2),'YData',Data(1));
set(htext,'Position',[Data(2),Data(1)]);

drawnow;