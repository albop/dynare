function [PI,VI,EI] = read_model_indices()
    global M_
    global oo_

    s = struct;
    for i = 1:length(M_.param_names(:,1))
	name = M_.param_names(i,:);
	value = i;
	s = setfield(s, name, value);
    end
    PI = s;
    
    s = struct;

    for i = 1:length(M_.endo_names(:,1))
	name = M_.endo_names(i,:);
	value = i;
	s = setfield(s, name, value);
    end
    VI = s;


    s = struct;
    for i = 1:length(M_.exo_names(:,1))
        name = M_.exo_names(i,:);
        value = i;
        s = setfield(s, name, value);
    end    
    EI = s;
    
end
