%%%% definition du probleme %%%%

A= [6 3 ; 3 4]
f= [-4; 2]
H= [3 2; -1 2]
B= [6 ;4]

%exemple du cours
%A= [2 -2 ; -2 4]
%f= [-2; -6]
%H= [1 1; -1 2]
%B= [2 ;2]

tstart = tic;
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA] = quadprog(A,f,H,B)
telapse = toc(tstart)
