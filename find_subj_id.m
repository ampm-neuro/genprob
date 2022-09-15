function subj_id_str = find_subj_id(string_in)
% outputs a string of 6 consecutive digits found in string input

% digit locations
dig_locs = regexp(string_in,'\d');

% find sequence of string size
str_size = 6;

N = str_size-1; % Required number of consecutive numbers following a first one
x = diff(dig_locs)==1;
f = find([false,x]~=[x,false]);
g = find(f(2:2:end)-f(1:2:end-1)==N,1,'first');
first_loc = dig_locs(f(2*g-1)); % First t followed by >=N consecutive numbers

% gather id
subj_id_str = string_in(first_loc : first_loc + N + 2);