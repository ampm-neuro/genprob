function output = percent_change(original_number, new_number, percentage_change)
% input 2 values and one empty bracket, function returns third value

if isempty(original_number)
    output = new_number/(1+(percentage_change/100));
elseif isempty(new_number)
    output = original_number + original_number*(percentage_change/100);
elseif isempty(percentage_change)
    percentage_change = (new_number/original_number);
    percentage_change = (1 - percentage_change) * 100;
    output = percentage_change;
end