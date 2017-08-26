% Roll_response_3L=Roll_response;
% Roll_response_undam
% Roll_response_no_add=Roll_response;

% Roll_response_2L=Roll_response;

figure
plot(Responses_out.Time,[Roll_response_no_add,Roll_response_3L(:,2),Roll_response_2L(:,2)])
legend('\fontsize{12} Roll rate comm (º/s)','\fontsize{12} Roll rate actual (º/s)')
ylabel('\fontsize{12} p (º/s)')
title('\fontsize{14} Roll rates')
