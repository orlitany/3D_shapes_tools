function [ST] = set_options(varargin)

if mod(length(varargin),2) ~= 0,
   error('Malformed fields, exiting.');
end

ST = [];
for p = 1:2:length(varargin),
    if ~isstr(varargin{p})
        error('Malformed fields. Field names must be strings.');
    else%if ~isfield(ST,varargin{p})
        ST = setfield(ST,varargin{p},varargin{p+1});
    end
end

return