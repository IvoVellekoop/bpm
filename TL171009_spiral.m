
%% 1) Calculating the center of mass of each focus in the spiral.
%  2) convolve each focus to a intensity profile (gaussian)
%  prerequisite file: Spiral,SpiralROI, Epoint
%  Warning:  Because of the convolution, the dimension of the output matrix will be bigger.
s=fieldnames(SpiralROI);
CenterOfMass=zeros(88,2);
threshold = 0.015;
addpixels = 50;
q = 1;
for i=1:size(s,1)
    if rem(i,12)==1
    else
        p=strcat('F_',num2str(i));
        temp=Spiral.*double(getfield(SpiralROI,s{i}));
        Sx=0; Sy=0; n=0;
        for ip=1:size(temp,1)
            for jp=1:size(temp,2)
                if temp(ip,jp)> threshold;
                    temp(ip,jp)=1;
                    Sx = Sx+ip;
                    Sy = Sy+jp;
                    n = n+1;
                else
                    temp(ip,jp)=0;
                end
            end
        end
        CenterOfMass(q,1)=round(Sx*(1/n))+addpixels;
        CenterOfMass(q,2)=round(Sy*(1/n))+addpixels;
        temp=zeros(1124,1124);
        temp(CenterOfMass(q,1),CenterOfMass(q,2))=1;
        q = q+1;
        temp=conv2(temp,abs(Epoint));
        str=sprintf('SingleField.%s=temp;',p);
        eval(str);
    end
end

%% Because of a bigger spiral file, all the medium files (Nair, Nscattered, SN_rev)
% need to be reshaped to be a larger size
NairBig=ones(1144,1144,400);
NairBig(61:1084,61:1084,201)=Nscattered(:,:,201); 

NscatBig=ones(1144,1144,400);
for i=1:400
    NscatBig(61:1084,61:1084,i)=Nscattered(:,:,i);    
end

SN_revBig=ones(1144,1144,400);
for i=1:400
    SN_revBig(61:1084,61:1084,i)=SN_rev(:,:,i);    
end


%% Propagates every focus in Nair to acquire corresponding plane wave
%  then back propagate the resulting plane wave through a reversed 
%  scattering medium (SN_revBig). Pick the final layer Eshaped_in and
%  and calculate its phase conjugation as the new shaped input for 
%  the forward scattering medium (NscatBig).  All shaped inputs are stored
%  in a structure called shaped_field
tic();
layer_thick=50;
Shaped_Field = struct;
s=fieldnames(SingleField);

for i=1:size(s,1)
    disp(s{i})
    F = Field(getfield(SingleField,s{i}), (20000/8)/1024, 0.9, 'um');
    Eout = F.propagate(NairBig,layer_thick*size(NairBig,3));
    Erev_in = Eout(:,:,size(Eout,3));
    clear Eout;
    Erev_in = Field(Erev_in.data, (20000/8)/1024, 0.9, 'um');
    Erev_out = Erev_in.propagate(SN_revBig,layer_thick*size(NscatBig,3));
    Eshaped_in = Erev_out(:,:,size(Erev_out,3));
    clear Erev_out;
    Edata = Eshaped_in.data;
    Eins = abs(Edata(:,:,size(Edata,3)));
    Edata = Edata(:,:,size(Edata,3))./Eins;
    Eshaped_conj=(1./Edata).*Eins;
    clear Edata;
    clear Eins;
    Eshaped_in = Field(Eshaped_conj, (20000/8)/1024, 0.9, 'um');
    clear Eshaped_conj;
    Shaped_Field = setfield(Shaped_Field,s{i},Eshaped_in);
end

toc()

%% Adding fields of all shaped input together
%  add a mask to prevent overlaping among different foci.
    structure = Shaped_Field; 
    if ~isstruct(structure)
        error('input needs to be a structure')
    end
    w = zeros(1144,1144);
    s = fieldnames(structure);
    for i = 1:size(fieldnames(structure))    
        temp = getfield(structure,s{i});
        w = w+temp.data;       
    end
    clear temp;
    % BigMask can be created from NonOverlayMask.m
    w=w.*BigMask;
    Emasked_spiral=Field(w,(20000/8)/1024, 0.9, 'um');
    clear w;
    Eout_maskedShaped=Emasked_spiral.propagate(NscatBig,400*50);

 
%% Phase step scanning
%  Prerequisite file: BigMask.
    dim = 1;
    
    Escan = zeros(1144,1144,11);
    
    for waven = 50:5:50
        Ephase_field = genPhaseStep(waven,Shaped_Field,CenterOfMass,dim);
        Ephase_field = Ephase_field.data.*BigMask;
        Ephase_field = Field(Ephase_field,(20000/8)/1024, 0.9, 'um');
        Eout_scanned = Ephase_field.propagate(NscatBig,size(NscatBig,3)*50);
        %Escan(:,:,(waven/5)+1) = Eout_scanned.data(:,:,400);
    end
    
    Escan = SizedArray(Escan, [(20000/8)/1024, 50], 'um');

