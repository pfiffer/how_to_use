function DBM=FDKK(nw,TFL)


add_function_paths();

       DBM=FDK(nw,TFL) ;      

end

function add_function_paths()
try
    functionname='FDK.m';
    functiondir=which(functionname);
    functiondir=functiondir(1:end-length(functionname));
%     addpath([functiondir '/FDK'])
addpath([functiondir ])
%     addpath([functiondir '/shortestpath'])
catch me
    disp(me.message);
end
end