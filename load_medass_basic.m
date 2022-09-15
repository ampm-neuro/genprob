function [medass_cell] = load_medass_basic(filepath)
%function [medass_cell, char_vars] = load_medass(filepath)
% Loads med associates datafiles (specified for the GEN programs)
%

% -----------------------------------
% ------------- INPUTS --------------
% -----------------------------------

% A = Variable vector (see below)
% B = Trial start times
% C = Reward availability
% D = Trial end times
% E = House light off times array
% F = Nose poke light off times array
% G = Nose_poke onset times array
% H = Nose_poke offset times array
% I = Tone onset times array
% J = Tone frequency array
% K = Delay duration vector
% L = Head_entry onset times array
% M = Head_entry offset times array
% N = Reward delivery times array
% O = Lickometer times array
% P = Tone offset array
% Q = Opto on array (probe opto sessions only) or valid minimic recording
% T = ChiSqr4 Delay durations (100)
% U = Speaker frequencies
% V = Reward probability distribution
% W = Randomizer (100)
% X = Randomizer (50)
% Y = Speaker ID
% Z = Speaker holding variable

% A(2) = Tone duration in seconds
% A(3) = Reward size (delivery duration)
% A(4) = Tone rise/fall
% A(5) = Tone amplitude
% A(6) = Delay duration (fixed component)
% A(7) = Timer
% A(8) = Trial counter
% A(9) = Nose_poke counter (Toggle-filtered)
% A(10) = Head_entry counter (Toggle-filtered)
% A(12) = Reward counter
% A(13) = Lickometer counter
% A(15) = Correct trial counter
% A(16) = Valid trial counter
% A(17) = Nose-poke level counter
% A(18) = Head-entry level counter
% A(19) = Randomization holder
% A(20) = Current delay duration (variable component)
% A(21) = Tone frequency
% A(22) = Reward Availability (momentarily probability)
% A(23) = Trial initiation counter
% A(24) = Variable reward timer
% A(25) = Response timer duration
% A(26) = Time out duration


% -----------------------------------
% ------------ OUTPUTS --------------
% -----------------------------------

% A medass_cell{1} = Assorted variable values (see A above)
% B medass_cell{2} = Trial start nose-poke times
% C medass_cell{3} = Reward availability (1=available; 0=probe)
% D medass_cell{4} = Trial end times
% E medass_cell{5} = House light off times array / start of random delay
% F medass_cell{6} = Nose poke light off times array / trial-start head entry
% G medass_cell{7} = Nose_poke onset times array
% H medass_cell{8} = Nose_poke offset times array
% I medass_cell{9} = Tone onset times array
% J medass_cell{10} = Tone frequency array
% K medass_cell{11} = Delay duration vector
% L medass_cell{12} = Head_entry onset times array
% M medass_cell{13} = Head_entry offset times array
% N medass_cell{14} = Reward delivery times array
% O medass_cell{15} = Lickometer times array
% P medass_cell{16} = Tone offset times array
% Q medass_cell{17} = Minimic video valid index / opto on
% R medass_cell{18} = [UNUSED]
% S medass_cell{19} = [UNUSED]
% T medass_cell{20} = Random delay duration vector
% U medass_cell{21} = Tone frequency vector
% V medass_cell{22} = Reward probability vector
% W medass_cell{23} = Randomizer (100)
% X medass_cell{24} = Randomizer (50)
% Y medass_cell{25} = Speaker ID
% Z medass_cell{26} = Speaker holding variable



%% load medass cell

    %read output (table)
    fid = fopen(filepath);
    medass_output = textscan(fid, '%s%s%s%s%s%s');

    %isolate session ID information
    end_id = find(strcmp(medass_output{1,1}, 'MSN:'));
    medass_id = [];
    for mi = 1:size(medass_output,2)
        medass_id = [medass_id medass_output{1,mi}(1:end_id)];
        medass_output{1,mi} = medass_output{1,mi}(end_id+1:end);
    end

    
    %parse by identifying med associates variables
    mav_idx = [];
    count = 0;
    medass_char_vars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    char_vars = [];

    %iterate through first column of output
    for imav = 1:size(medass_output{1,1},1)
        count = count+1;

        %identify med associates variables
        if (contains(medass_char_vars, medass_output{1,1}{imav,1}(1))...
                && strcmp(medass_output{1,1}{imav,1}(2), ':') &&...
                length(medass_output{1,1}{imav,1})==2) || ...
                imav==size(medass_output{1,1},1)
            
            %record variable char
            char_vars = [char_vars medass_output{1,1}{imav,1}(1)];

            char = medass_output{1,1}{imav,1}(1);
            
            %keep track of row and variable id
            if imav<size(medass_output{1,1},1)
                mav_idx = [mav_idx; count];
            else
                %weirdness for last variable
                mav_idx = [mav_idx; count+1]; 
            end

            if length(mav_idx)>1
                %gather data related to that variable
                hold_var_data = [];
                for hvdi = 2:size(medass_output,2)
                    hold_var_data = [hold_var_data ...
                        medass_output{1,hvdi}(mav_idx(end-1):mav_idx(end)-1)];
                end

                hold_var_data = hold_var_data';
                hold_var_data = str2double(hold_var_data(:));
                hold_var_data = hold_var_data(~isnan(hold_var_data));

                %load output cell sequentially with each med
                %associates character variable (A is 1, etc.)
                medass_cell{length(mav_idx)-1} = hold_var_data;
            end

        end
    end

    %alphabetize character variables
    char_vars = char_vars(1:end-1);
    [char_vars, sort_idx] = sort(char_vars);
    medass_cell = medass_cell(sort_idx);
    
   
    
