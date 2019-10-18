%Resolution de probleme quadratiques par l'algorithme de Lemke
%
% Felix CHOI, 2019
%
% le probleme est suppose de forme : | min (1/2)<Ax,x> + <b,x>
%                                    |      <H,x> <= B
%L'algorithme sortira 2 vecteurs z=[u; x] et w=[y; v]

%%%% definition du probleme %%%%

A= [6 3 ; 3 4]
b= [-4; 2]
H= [3 2; -1 2]
B= [6 ;4]



%exemple du cours
%A= [2 -2 ; -2 4]
%b= [-2; -6]
%H= [1 1; -1 2]
%B= [2 ;2]


%%%% Reecriture sous forme LCP
[mA,nA] = size(A); %taille de A
[mH,nH] = size(H); %taille de H

%%% M
zero = zeros(nA);
M = [zero -H; H' A]
[mM,nM] = size(M);

%%% q
q= [B; b]
[mq,nq] = size(q);

%variable a trouver
w = zeros(mq,1);
z = zeros(mq,1);


%---------------------------%
%%%% Algorithme de Lemke %%%%
%---------------------------%
tstart = tic;

z0 = -1 * ones(mq);
Z = -M;

W = eye(mq);

%matrice indiquant les entrees et sorties en base
temp = zeros(mq,1);
for i=1:mq
    temp(i) = i;
end
base = [ones(mq,1)  temp] %colonne 1 pour le vecteur correspondant 
                            %(1 pour w et 0 pour z)
                            %colonne 2 pour l'indice
                                  

tableau = [ W Z z0 q]

%--------------------------%
%%%%%%% iteration 1: %%%%%%%
k = find( q == min(q(:)))  %ligne du pivot
k0 = k;                     %ligne initial du pivot -> utilise pour critere
                            %d'arret

                            % pivot -> z0(k)
base(k,:) = [0 0];
% -- operation pivot -- %

%reduction de la ligne du pivot
a = 1/z0(k);       % a =-1
W(k,:) = a*W(k,:);
Z(k,:) = a*Z(k,:);
z0(k)  = a*z0(k);
q(k)   = a*q(k);

for i=1:mq
    if i ~= k               % ~= -> !=
        a = -z0(i)/z0(k) ; %on doit avoir "z0(i) + a*z0(k) = 0"        
        W(i,:) = W(i,:) + a*W(k,:);
        Z(i,:) = Z(i,:) + a*Z(k,:);
        z0(i)  = z0(i)  + a*z0(k);
        q(i)   = q(i)   + a*q(k);
    end
end

tableau = [ W Z z0 q];


%--------------------------%
% iterations jusqu'a condition d'arret -> z0 sortant, donc si on tombe sur
% k == k0

k_mem = k;       %on garde en memoire la ligne du pivot precedent
k=0    ;         %on entre dans la boucle k0 ne vaut jamais 0

nb_iter = 1;
while k ~= k0,
    vec_temp = q./Z(:,k_mem);
    vec_temp(vec_temp  <= 0) = inf  ;  % on ecarte les valeurs <= 0
    k = find( vec_temp == min(vec_temp(:)));%on cherche le pivot=Z(k,k_mem)
    Z(k,k_mem);
    % -- operation pivot -- %

    %reduction de la ligne du pivot
    a = 1/Z(k,k_mem);       %on doit avoir "a*Z(k,k_mem) = 1"
    W(k,:) = a*W(k,:);
    Z(k,:) = a*Z(k,:);
    z0(k)  = a*z0(k);
    q(k)   = a*q(k);

    for i=1:mq
        if i ~= k               
            a = -Z(i,k_mem)/Z(k,k_mem) ; %on doit avoir 
                                         %"Z(i,k_mem) + a*Z(k,k_mem) = 0"        
            W(i,:) = W(i,:) + a*W(k,:);
            Z(i,:) = Z(i,:) + a*Z(k,:);
            z0(i)  = z0(i)  + a*z0(k);
            q(i)   = q(i)   + a*q(k);
        end
    end
    base(k,:) = [0, k_mem];
    
    tableau = [ W Z z0 q];
    
        
    
    
    k_mem = k;
    
    %%% limiteur d'iteration %%%
    nb_iter = nb_iter + 1;
    if nb_iter > 10
        k=k0;
        disp("maximum iteration reached.")
    end
end
telapse = toc(tstart);
    

%interpretation de la matrice artificielle "base" pour determiner w et z
for i=1:mq
    x = base(i,:);
    if x(1) == 0
        z(x(2)) = q(i);
    elseif x(1) == 1
        w(x(2)) = q(i)  ;
    end
end
disp("RESULT :")
z
w

nb_iter
telapse



    
    









