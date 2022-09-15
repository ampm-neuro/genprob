function [problem_tones] = rich_bounds_prob(problem_set_str, problem_numbers, varargin)
% plot rich tone locations for all 6 problems of an input problem set
% use zero as problem_numebrs input to avoid plot

if nargin == 1
    problem_numbers = 1:6; % all problem numbers
    freq_or_num = 0; % freq
elseif nargin == 2
    freq_or_num = 0; % freq
elseif nargin == 3
    freq_or_num = varargin{1}; % 0 freq, 1 num
end

freq_or_num

% load tones (rich and poor)
load('unqfrq41', 'unqfrq41');
problem_tones = load('rich_tones_mtx_2_var_5.mat', ['rwd1_nrwd_seqence_' problem_set_str]);
problem_tones = cell2mat(struct2cell(problem_tones));
%problem_tones = problem_tones(problem_numbers,:);

% use zero to avoid plot
if problem_numbers ~= 0
    
    % colors
    rich_colors = [153 0 0 ; 255 0 0] ./ 255;
    poor_colors = [0 0 153 ; 0 0 255] ./ 255;


    % plot
    hold on;
    for iprob = problem_numbers

        if iprob > 6
            continue
        end

        if freq_or_num == 0
            % rich
            plot(problem_tones(iprob, 1).*[1 1], ylim, '-', 'color', rich_colors(1,:));

            % poor
            plot(problem_tones(iprob, 2).*[1 1], ylim, '-', 'color', poor_colors(1,:));

        elseif freq_or_num == 1

            % rich
            plot(find(unqfrq41==problem_tones(iprob, 1)).*[1 1], ylim, '-', 'color', rich_colors(2,:));

            % poor
            plot(find(unqfrq41==problem_tones(iprob, 2)).*[1 1], ylim, '-', 'color', poor_colors(2,:));
        end




    end

end

if freq_or_num == 1
    for itone = 1:numel(problem_tones)
        problem_tones(itone) = find(unqfrq41==problem_tones(itone));
    end
end