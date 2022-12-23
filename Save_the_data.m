%% store all the vector in the table
clc
clear
a = [];
save('tension.mat','a');
if isfile('tension.mat')
    delete('tension.mat');
    a = [];
    save('tension.mat','a');
else
     % File does not exist
     a = [];
     save('tension.mat','a');
end
 

for i = 1:10
   y  = myfunc(i) ;
end

load('tension.mat','a');



function y = myfunc(i)
     load('tension.mat','a');
     a = [a, i]; 
     save('tension.mat','a');
     y = i;
end