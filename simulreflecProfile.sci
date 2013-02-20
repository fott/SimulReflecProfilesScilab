function[interface,profil] = ReadML(fic_data)
// Mode d'affichage
mode(0);

// Affiche un avertissement pour une exception en virgule flottante
ieee(1);

// lecture du fichier décrivant la ML 
// format pour chaque layer :
// substrat	0	N	0	N  0 N
// profil1		linear
//
// ajout d'une ligne en fin de fichier
//unix("cat "+fic_data+" blanc.txt > data.txt");
fid = mopen(fic_data);
//
layer=mfscanf(-1,fid,' %s %e %s %e %s %e %s \n %s %s');
mclose(fid);
// on cree un tableau layer(i,j) : i numero couche 
//                                 j=1 --> 9 (1--> 7 couche nom, ....)
//                                            8 et 9 --> profil et type de profil
// 
disp(layer);
sz=size(layer);
nbInterfaces = sz(1);// nombre d'interfaces lues = nb layer 
nbProfiles = nbInterfaces-1;// nombre d''interfaces (séparant des interfaces)
//
// creation de la structure de la multicouche
// la liste des interfaces
// interface[i].name
//               .position (with respect to the substrate surface)
//               .SLD
interface = struct("name","FirstInterface","position",0,"positionFit",0,"SLD",0,"SLDfit",0);
for i = 1:nbInterfaces
  interface(i).name = layer(i,1);
  interface(i).position = layer(i,2);
  interface(i).positionFit = layer(i,3);
  interface(i).SLD = layer(i,4);
  interface(i).SLDfit = layer(i,5);
  interface(i).SLDi = layer(i,6);
  interface(i).SLDifit = layer(i,7);
end;
// la liste des profils
// profils[i].name
//           .type
profil = struct("name","FirstProfile","type","ProfileType");
disp("********************");
for i = 1:nbProfiles
  profil(i).name = layer(i,8);
  profil(i).type = layer(i,9);
end;
endfunction

function[slice] = generateSlice(interface,profil)
// Mode d'affichage
mode(0);

// Affiche un avertissement pour une exception en virgule flottante
ieee(1);

// generation de la multicouche à partir de
// interface, profil, discretizationStep
// on démarre coté substrat
discretizationStep = 1;// épaisseur des couches du profile valeur en nm
slice = struct("depth",0,"SLD",0,"SLDi",0);
// la premiere couche décrit le vide, 
// la dernière couche décrit le substrat
sliceIdx = 1;
depthCrt = 0;
nbProfiles=size(profil,1);
nbInterfaces = size(interface,1);
for i=1:nbProfiles
    SLD0=interface(i).SLD;
    SLDi0=interface(i).SLDi;
    SLD1=interface(i+1).SLD;
    SLDi1=interface(i+1).SLDi;
    thickness=interface(i+1).position-interface(i).position;
    profileType=profil(i).type;
    nbSlices=int(thickness/discretizationStep); // entier
    for j=0:nbSlices-1
        if (profileType == "flat")
            slice(sliceIdx).SLD=SLD0;
            slice(sliceIdx).SLDi=SLDi0;
            slice(sliceIdx).depth=depthCrt;
        end 
        if (profileType =="linear")
            slice(sliceIdx).SLD=SLD0+j*(SLD1-SLD0)/(nbSlices);
            slice(sliceIdx).SLDi=SLDi0+j*(SLDi1-SLDi0)/(nbSlices);
            slice(sliceIdx).depth=depthCrt;
        end
        if (profileType == "erf")
            slice(sliceIdx).SLD=SLD0+(erf(4*(j-nbSlices/2)/nbSlices)+1)/2*(SLD1-SLD0);
            slice(sliceIdx).SLDi=SLDi0+(erf(4*(j-nbSlices/2)/nbSlices)+1)/2*(SLDi1-SLDi0);
            slice(sliceIdx).depth=depthCrt;
        end
        sliceIdx=sliceIdx+1;
        depthCrt=depthCrt+discretizationStep;
    end
end
// on finit par le subtrat
slice(sliceIdx).depth=depthCrt;
slice(sliceIdx).SLD=interface(nbInterfaces).SLD;
slice(sliceIdx).SLDi=interface(nbInterfaces).SLDi;
nbLayers=sliceIdx;
endfunction

function[ ] = plotProfile(interface,slice)
// tracé du profil d'indice optique à partir de la couche Slice
// Mode d'affichage
mode(0);

// Affiche un avertissement pour une exception en virgule flottante
ieee(1);

// tracé du profil d'indice optique à partir de la couche slice
clear("depthPosition","SLD","SLDi");
clf(1);
scf(1);
nbLayers = size(slice,1);
nbInterfaces = size(interface,1);

for i=1:nbLayers
    depthPosition(i)=slice(i).depth;
    SLD(i)=slice(i).SLD;
    SLDi(i)=slice(i).SLDi;
end
plot(depthPosition,SLD,'r-',depthPosition,SLDi,'g-')

// tracé des position des interfaces
for i=1:nbInterfaces
    plot([interface(i).position,interface(i).position],[0,interface(i).SLD],'b');
    //h(i)=impoint(gca,interface(i).position,interface(i).SLD);
    //addNewPositionCallback(h(i), updateProfile);
end
xtitle("SLD profile","depth","SLD");
//h=impoint(gca,20,5);
//addNewPositionCallback(h,@(h) title(sprintf('(%1.0f,%1.0f)',h(1),h(2))));
endfunction

function [Q,Ref,Transr] = calculateReflectivity(slice)

// reflectivity calculation using the basic recursive formalism
// notations as in F. OTT, PhD 1998
// lengths in [nm]
// reminder: SLD in [10-6 A°-2]
// we work in the reciprocal space so as to be more universal
// units [nm-1] so as to make the values more readable
// ON DUPLIQUE LE CODE DELPHI
// ET ON VERIFIE LA SELF CONSISTENCE DES VALEURS au débugage
//
discretizationStep = 1;// épaisseur des couches du profile valeur en nm
nbLayers = size(slice,1); // nb slice
Q=0.01:0.01:1;
Kz0=Q./2; 
zz=discretizationStep;
nbRefPoints=size(Q,2);
Ref=zeros(nbRefPoints,1); // reflectivity
Transr=zeros(nbRefPoints,1);
// definition du tableau SliceSld
for idxL=1:nbLayers
    SliceSld(idxL) = slice(idxL).SLD;
end
// wavevectors in the different layers
for idxk0=1:nbRefPoints // balayage des vecteurs d'ondes incidents
     // kn²=k0²-4*pi*rho*b 
    Kz=-sqrt(Kz0(idxk0)^2-4*%pi*SliceSld*1E-6*100); 
    // this value can  be real or complex
    // first coefficient 1E-6 to go back to A°-2 values; 
    // second coefficient to get nm-2   
    // ***********************************
    // remplissage des matrices de transfert d'un milieu à l'autre
    // formules I.120
    // on démarre la récursion par le dessus
    // copie du code de SimulReflec
    qq1=Kz(nbLayers);
    M=eye(2,2);
    for idxL=nbLayers-1:-1:1, // boucle du vide vers le substrat 
        qq0=Kz(idxL);
        pp=(qq1+qq0)/qq1/2;
        mm=(qq1-qq0)/qq1/2;
        // matrice D de propagation
        D = [pp mm; mm pp]*[ exp(%i*qq0*zz) 0 ; 0 exp(-%i*qq0*zz)] ;
        M=M*D ;
        qq1=qq0;
    end
    auxr=(abs(M(:,1)))^2;
    Ref(idxk0)=auxr(2)/auxr(1) ;
    Transr(idxk0)=1/auxr(1) ;
end

endfunction

function [] = plotReflectivity(Q,Ref,Transr);
//  graphe reflectivite
    clf(2);
    scf(2);
    a = gca();
    a.labels_font_size=1;
    a.auto_clear= "on";
    plot2d(Q,[Ref Transr]);
    a.log_flags = "nln" ; // set Y axes to logarithmic scale
    xtitle("Reflectivity curve","Q","Ref");
endfunction


function[]= simulreflecProfile(fic_data);
    // exec du pg de simulreflec
    [interface,profil] = ReadML(fic_data);
    [slice] = generateSlice(interface,profil);
    plotProfile(interface,slice);
    [Q,Ref,Transr] = calculateReflectivity(slice) ;
    plotReflectivity(Q,Ref,Transr);
endfunction

// exemple d execution
fic_data="MyMultiLayer_NiSi.txt";
simulreflecProfile(fic_data);
