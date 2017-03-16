function varargout = parseKeyValuePairs(userKeyNamePairs,valid_keys,default_values)
% PARSEKEYVALUEPAIRS accepts a calling function's 'varargin' cell and interprets its key-value pairs. 
%
%	INPUT
%       userKeyNamePairs: The 'varargin' cell from another function
%       valid_keyNames: A cell array of valid key names. Passed in from calling function.
%       default_values: A cell array of the default values for each key name.
%
%	OUTPUT
%       varargout: Values of the key-name pairs, in the order specified by 'valid_keyNames'.
%
%   TODO
%       Combine the 'valid_keyNames' and 'default_values' into a single variable
%
% John W. Miller
% 2015-07-31

% Check for the proper dimensions
if mod(length(userKeyNamePairs),2)
    error('Variable names and their values must come in pairs.')
end
valid_keys = lower(valid_keys);

% Extract the key names and values
keyNames_FromUser = lower(userKeyNamePairs(1:2:end));
values_FromUser   =       userKeyNamePairs(2:2:end);

% Check that the user-supplied key names are valid (Excessive, I know...)
name_check = ismember(keyNames_FromUser,valid_keys);
if ~all(name_check)    
    errmsg = sprintf('Invalid key names: %s', printcell(keyNames_FromUser(~name_check),', '));
    error('%s\nValid key names:\n%s',errmsg,printcell(valid_keys))
end

% Replace the default values with user-supplied values
processed_vars  = containers.Map(valid_keys,default_values,'UniformValues',false);
n_vars_FromUser = length(keyNames_FromUser);
for n_var_FromUser = 1:n_vars_FromUser
    processed_vars(keyNames_FromUser{n_var_FromUser}) = values_FromUser{n_var_FromUser};
end

% Varargout
varargout = values(processed_vars,valid_keys); % Must specify the order of keys to pull out

end % End of main