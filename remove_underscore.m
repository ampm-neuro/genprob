function string_in = remove_underscore(string_in)
% replaces underscores with spaces in a string
% e.g. "dance_in_the_rain" becomes 'dance in the rain'

string_in(strfind(string_in, '_')) = ' ';