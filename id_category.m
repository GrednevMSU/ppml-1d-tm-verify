function c = id_category(id)
if     ~isempty(strfind(id,'scenario')),     c='scenario';
elseif ~isempty(strfind(id,'branch')),       c='branch';
elseif ~isempty(strfind(id,'sigma')),        c='conducting_interface';
elseif ~isempty(strfind(id,'out_of_domain')),c='out_of_domain';
elseif ~isempty(strfind(id,'rand')),         c='randomized';
else                                         c='other';
end
end
