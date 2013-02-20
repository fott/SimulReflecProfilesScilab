function[X] = init(dx1,dx2,l)
  // initialisation des particules ds [0,l]x[0,l/2]
  //-----------------------------------------------
  //dx1 = interatomic distance ds milieu 1
  //dx2 = interatomic distance ds milieu 2
  // l/2 = longueur carre
  // X = positions initiales des particules
  //    X(i,1) = x
  //    X(i,2) = y
  //    X(i,3) = 1/2 (milieu)
  // --------------------------------------------------
  // particules ds milieu 1 (carre [0,l/2]x[0,l/2])
  n1 = int(l/(2*dx1));
  k=1;
  for i=1:n1
    for j=1:n1
      X(k,1)=i*dx1;
      X(k,2)= j*dx1;
      X(k,3) = 1; 
      k=k+1;
    end
  end  
  // particules ds milieu 2 (carre [l/2,l]x[0,l/2])  
  n2 = int(l/(2*dx2));  
  for i=1:n2
    for j=1:n2
      X(k,1)=l/2+i*dx2;
      X(k,2)= j*dx2;
      X(k,3) = 2; 
      k=k+1;
    end
  end  
endfunction

function [y] = f(x)
  // loi du centre de spike (en x)
  y = l*x ;  // uniforme sur [0,l]
endfunction

function [y] = f_fd(x)
  // loi du centre de spike 
  // (selon la distribution de fd energie deposee)
  if((0<x) & (x<0.808)) then
    delta = 3.33e-5 - 3.46e-5*x ;
    y = (5.77e-3 - sqrt(delta))/1.73e-5 ;
  else
    if((0.808<x) & (x<1)) then   
      delta = 6.92e-6*(1-x) ;
      y = (1.845e-3 - sqrt(delta)) / 3.46e-6 ;
    else
      y=0
    end
  end
endfunction

               

function [y] = g(x)
  // loi du centre de spike (en y)
   y = l/2 * x; // uniforme sur [0,l/2]
endfunction

function [c] = center_spike(type_spike)
  // on tire un centre de spike
  // -----------------------------
  // cas uniforme 
  // en x suivant la loi f(x) ds[0,l]
  // en y suivant la loi g(x) ds [0,l/2]
  //-------------------------------------
  if(type_spike=="uniform") then
    u = rand() ;
    cx = f(u) ;
  else
    u = rand() ;
    cx = f_fd(u) ;
  end      
  u = rand();
  cy = g(u);
  c(1,1)= cx ;
  c(1,2)= cy ;
endfunction

function [y] = hr(x)
  // 
   y = -sqrt(rc)*log(1-x) ;
 endfunction
 
function [y] = hphi(x)
  // 
   y = 2*%pi*x ;
 endfunction
  

function [d] = displ(rc,n)
  // calcule les deplacements de tous les atomes 
  // suivant une loi p(r) = exp(-r/sqrt(rc))/sqrt(rc) en r 
  // uniforme en angle 
  // ---------------------------------------
  u = rand(n,2) ;
  r = hr(u(1:n,1)) ;
  phi = hphi(u(1:n,2));
  dx = r.*cos(phi) ;
  dy = r.*sin(phi) ;
  d(1:n,1)= dx ;
  d(1:n,2)= dy ;  
endfunction

function [S] = dist(Z,c) 
  // distance entre ts atomes et le centre c
  [n,nn]= size(Z) ;
  S = sqrt((Z(1:n,1)-c(1,1))^2+(Z(1:n,2)-c(1,2))^2);
endfunction


function [Y] = int_true_tab(B)
  // renvoie 1 si true 0 si false
  [n,nn]=size(B) ;
  for i=1:n
    for j=1:nn
      if (B(i,j)) then
        Y(i,j) = 1;
      else
        Y(i,j) = 0;
      end
    end
  end
endfunction

function [Y] = period(Z,l)
  // cns periodiques en y
  [n,nn]=size(Z);
  Nl = floor(2*Z(1:n,2)/l) ;
  Y(1:n,1) = Z(1:n,1) ;
  Y(1:n,2) = Z(1:n,2)-(Nl*l/2) ;
endfunction

  function [Zd] = displace(c,rc,Z)
  // deplacements de tous les atomes
  // -------------------------------
  // deplace un atome si il est 
  // ds le cercle de centre c et de rayon sqrt(rc) 
  // sinon ne le deplace pas
  // ------------------------------------------------
  [n,nn]=size(Z);
  Zd = zeros(Z) ;
  // tab des deplacements  
  [D] = displ(rc,n) ;
  // tab des distances entre un atome et le centre c
  [S] = dist(Z,c)  ;
  // tab des atomes ds le cercle de centre c et de rayon sqrt(rc) 
  // (A=1 --> on le deplace A=0 --> pas de deplacement )
  [A] = int_true_tab((S<sqrt(rc))) ;
  // deplacement des atomes
  Zd(1:n,1) = Z(1:n,1)+ D(1:n,1) .* A  ;
  Zd(1:n,2) = Z(1:n,2)+ D(1:n,2) .* A  ; 
  // cns periodiques en y
  [Zd] = period(Zd,l) ;  
endfunction



function [Y,c] = displace_spike(X,rc)
  // choisit un centre de spike et deplace tous les points
  // ds cercle centre c rayon sqrt(rc)
  // -------------------------------------
  [n,np] = size(X);
  // choisit un centre de spike
  [c] = center_spike("nonuniform") ;
  Y = zeros(X) ; 
  // deplace tous les points ds cercle centre c rayon sqrt(rc)
  [Y(:,1:2)] = displace(c,rc,X(:,1:2)) ;
  Y(:,3) = X(:,3) ; 
endfunction


function [] = plot_mixing_spike(X,c,ks)
  //  trace atomes deplaces
    [n,np]=size(X) ;
    np1 = (int(l/(2*dx1)))^2; // ts les points du milieu 1
    scf(ks);
    plot(X(1:np1,1),X(1:np1,2),'b') ; 
    plot(X(np1+1:n,1),X(np1+1:n,2),'r') ;     
    //plot cercle spike
    t=[0:0.1:2*%pi]';
    plot(c(1,1)+rc*cos(t),c(1,2)+rc*sin(t),'g') ;
    plot(c(1,1)+sqrt(rc)*cos(t),c(1,2)+sqrt(rc)*sin(t),'g') ;
  endfunction
  
  function [] = plot_mixing(X)  
    // trace atomes deplaces
    [n,np]=size(X) ;
    np1 = (int(l/(2*dx1)))^2; // ts les points du milieu 1
   scf(0);
   plot(X(1:np1,1),X(1:np1,2),'b') ; 
   plot(X(np1+1:n,1),X(np1+1:n,2),'r') ;  
 endfunction
       
function [H,C1,C2]=concentration(X)
// concentration des atomes 1
// ------------------------------
// h : pas pour l histogramme (en x )
h= 7.1 ;
np = int(l/h) ; // nb points en x 
Cc1=zeros(np,1);
Cc2=zeros(np,1);
[n,nn]=size(X) ; // n = nb total atomes
np1 = (int(l/(2*dx1)))^2; // ts les points du milieu 1
// nb atomes 1 
for i=1:np1
  ih = int(X(i,1)/h) ;
  if ((ih <= np) & (ih >=1)) then 
    Cc1(ih)=Cc1(ih)+1 ;
  end  
end
// nb atomes 2
for i=np1+1:n
  ih = int(X(i,1)/h) ;
  if ((ih <= np) & (ih >=1)) then 
    Cc2(ih)=Cc2(ih)+1 ;
  end  
end
// X 
H=(1:h:(np*h)) ;
// concentrations 
C1 = Cc1 ./ (Cc1 + Cc2) ;
C2 = Cc2 ./ (Cc1 + Cc2) ;
endfunction


       
function [H,C1,C2]= mixing(Ns,rc)
  // fonction mixing
  // --------------------------------------
  // data
  l = 100 ; // (A) lg de boite de simulation
  dx1 = 2.5 ;  // interatomic distance ds 1 (Cr)
  dx2 = 2.5 ;  // interatomic distance ds 2 (Si)
  // initialisation particules
  [X] = init(dx1,dx2,l) ;
  [n,np]=size(X) ;
  np1 = (int(l/(2*dx1)))^2; // np1 premiers points ds milieu 1
  // deplacement des atomes sur Ns simulations
  for ks = 1:Ns 
    [Y,c] = displace_spike(X,rc) ;
    X = Y ;
    // trace atomes deplaces
   //plot_mixing_spike(X,c,ks) ;
 end  
 // trace atomes deplaces a la fin des iterations
 //plot_mixing(X) ; 
 // calcul des concentrations 
[H,C1,C2]=concentration(X) ;
endfunction


function y=FF(x,p)
  // fonction fit
     y = 1- ((1/2)+(1/2)*erf((x-p(1))/p(2))) ;
endfunction

  //The criterion function moindres carres
function e=G(p,z)
  x=z(1),y=z(2);
  e=(y-FF(x,p))^2;
endfunction


function [p,err] = fit_erf(H,C1)
  // calcul des fits de la fn de mixing
  //  ---------------------------------------
  // lecture fichier concentration Cr mixing
  // dir_data = 'J:\diffusion_entre_deux_phases\mixing\' ;  // usb
    X = H ;
    Y = C1' ;       
    Z = [X;Y] ;    
    //Solve the problem
   // p0 valeurs initiales
    p0=[50;5];
    [p,err]=datafit(G,Z,p0);   
  endfunction 
  
  function [p] = fit_mixing(Ns,rc)
    // calcul de la moyenne des parametres du fit de mixing 
    // sur 10 valeurs
    // -----------------------------------------------------
    p_mean = [0;0];
    for km = 1:10 
      [H,C1,C2]= mixing(Ns,rc) ;
      [p,err] = fit_erf(H,C1) ;
      p_mean = p_mean+p ;
    end
    p = p_mean /10;
  endfunction
  
      
      