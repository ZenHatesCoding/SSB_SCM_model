function output = SAW(input)

% if (input <= 1) && (input >= -1)
%     output = input;
% elseif (input > 1) && (input < 3)
%     output = input - 2;
% elseif (input < -1) && (input > -3)
%     output = input + 2;
% end



% function output = SAW(input)
% 
% if (input <= 0.5) && (input >= -0.5)
%     output = input;
% elseif (input > 0.5) && (input < 1.5)
%     output = input - 1;
% elseif (input < -0.5) && (input > -1.5)
%     output = input + 1;
% end


output = 2*(0.5*input+0.5) - 2*floor(0.5*input+0.5) - 1;