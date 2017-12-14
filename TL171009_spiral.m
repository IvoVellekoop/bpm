%% System parameters: 
focal_length = 10000; % of the lens
field_size = (20000/1024)*1144; 
layer_thick = 100; % Thickness of every section in Z (µm)
scattering_coefficient = 0.15; % Here and the next row defines how scatter the scattering medium is. 
GaussFilterSigma = 25;
Scat_Layer_No = 50; % Layer number of the scattering medium
Additional_layers = 0; % number of additional layers of vaccum after the scattering medium

format long

%% Generate medium, lens...
% Creates field
E = Field.plane(field_size/8,[1144 1144],0.9,'um');
E = E.aperture('gaussian',512);
% make lens
Lens = E.make_lens(focal_length,layer_thick);
% make scattering medium
SN = make_medium(E,scattering_coefficient,GaussFilterSigma,Scat_Layer_No);
N = ones(size(E,1),size(E,2),(2*focal_length/layer_thick)+Additional_layers);
% Combines everything together
N(:,:,(focal_length/layer_thick)+1) = Lens;
Nscattered(:,:,1:size(N,3)-Scat_Layer_No-Additional_layers)=N(:,:,1:size(N,3)-Scat_Layer_No-Additional_layers);
Nscattered(:,:,size(N,3)-Scat_Layer_No-Additional_layers+1:size(N,3)-Additional_layers)=SN;
Nscattered(:,:,size(N,3)-Additional_layers+1:size(N,3))=ones(size(E,1),size(E,2),Additional_layers);

% Reverse the order of the whole medium for back propagation
SN_rev = ones(size(Nscattered,1),size(Nscattered,2),size(Nscattered,3));
for j=0:size(Nscattered,3)-1
    SN_rev(:,:,j+1)=Nscattered(:,:,size(Nscattered,3)-j);
end

%% 1) Calculating the center of mass of each focus in the spiral.
%  2) convolve each focus to a intensity profile (gaussian)
%  prerequisite file: Spiral,SpiralROI, Epoint (Ask Tzu-Lun)
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
                if temp(ip,jp)> threshold
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

Layer_of_Lens = (focal_length/layer_thick)+1;


NairBig=ones(size(temp,1),size(temp,2),size(Nscattered,3));
dim_offset = ((size(temp,1)-size(Nscattered,1))/2);
NairBig(dim_offset+1:dim_offset+size(Nscattered,1),dim_offset+1:dim_offset+size(Nscattered,1),Layer_of_Lens)=Nscattered(:,:,Layer_of_Lens); 

NscatBig=ones(size(temp,1),size(temp,2),size(Nscattered,3));
for i=1:size(Nscattered,3)
    NscatBig(dim_offset+1:dim_offset+size(Nscattered,1),dim_offset+1:dim_offset+size(Nscattered,1),i)=Nscattered(:,:,i);    
end

SN_revBig=ones(size(temp,1),size(temp,2),size(Nscattered,3));
for i=1:size(Nscattered,3)
    SN_revBig(dim_offset+1:dim_offset+size(Nscattered,1),dim_offset+1:dim_offset+size(Nscattered,1),i)=SN_rev(:,:,i);    
end



 
%% Propagates every focus in Nair to acquire corresponding plane wave
%  then back propagate the resulting plane wave through a reversed 
%  scattering medium (SN_revBig). Pick the final layer Eshaped_in and
%  and calculate its phase conjugation as the new shaped input for 
%  the forward scattering medium (NscatBig).  All shaped inputs are stored
%  in a structure called shaped_field
tic();
Shaped_Field = struct;
s=fieldnames(SingleField);

for i=1:size(s,1)
    disp(s{i})
    F = Field(getfield(SingleField,s{i}), (20000/8)/1024, 0.9, 'um');
    F.gpu_enabled = false;
    tic()
    Eout = F.propagate(NairBig,layer_thick*size(NairBig,3));
    Erev_in = Eout(:,:,size(Eout,3)-Additional_layers);
    clear Eout;
    Erev_in = Field(Erev_in.data, (20000/8)/1024, 0.9, 'um');
    Erev_in.gpu_enabled = false;
    Erev_out = Erev_in.propagate(SN_revBig(:,:,Additional_layers+1:size(SN_revBig,3)),layer_thick*size(NscatBig,3));
    toc()
    Eshaped_in = Erev_out(:,:,size(Erev_out,3)-Additional_layers);
    clear Erev_out;
    Edata = Eshaped_in.data;
    Eins = abs(Edata(:,:,size(Edata,3)));
    Edata = Edata(:,:,size(Edata,3))./Eins;
    Eshaped_conj=(1./Edata).*Eins;
    clear Edata;
    clear Eins;   
    Eshaped_in = Field(Eshaped_conj,(20000/8)/1024, 0.9, 'um');
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
    
    %w=w.*BigMask;
    Ein_spiral=Field(w,(20000/8)/1024, 0.9, 'um');
    clear w;
    Eout_ShapedSpiral=Ein_spiral.propagate(NscatBig,200*100);

 

    
%% Phase step scanning
function [Escan,EscanIn] = StepPhaseScanning(MaxXw,MaxYw,Sort_Shaped_Field,SingleMasks,CenterOfMass,Nscat)
%  Prerequisite file: BigMask.
    count = 1;
    
    Escan = zeros(1144,1144,169);
    EscanIn = zeros(1144,1144,169);
    
    for wavenX = MaxXw:-5:-MaxXw
        for wavenY = -MaxYw:5:MaxYw
            tic
            Ephase_field = genPhaseStep(wavenX,wavenY,Sort_Shaped_Field,SingleMasks,CenterOfMass);
            %Ephase_field = Ephase_field.data.*BigMask;
            Ephase_field = Field(Ephase_field.data,(20000/8)/1024, 0.9, 'um');
            Eout_scanned = Ephase_field.propagate(Nscat,size(Nscat,3)*100);
            EscanIn(:,:,count) = Eout_scanned.data(:,:,1);
            Escan(:,:,count) = Eout_scanned.data(:,:,200);
            toc
            count=count+1;
        end
    end
    
    Escan = SizedArray(Escan, [(20000/8)/1024,(20000/8)/1024, 100], 'um');
    EscanIn = SizedArray(EscanIn, [(20000/8)/1024,(20000/8)/1024, 100], 'um');
end
    
    

%% Scanning by phase gradient
function [Escan,EscanInT] = GradPhaseScanning(Eshaped_input,MaxXw,MaxYw,Nscat)
count = 1;    
Escan = zeros(1144,1144,121);
EscanInT = zeros(1144,1144,121);
for wavenX = MaxXw:-7:-MaxXw
        for wavenY = -MaxYw:7:MaxYw
            phaseX = linspace(-2*pi*wavenX/2,2*pi*wavenX/2,1144);
            phaseX = phaseX'.*ones(1144,1);
            phaseY = linspace(-2*pi*wavenY/2,2*pi*wavenY/2,1144);
            phaseY = ones(1144,1)'.*phaseY;
            Etilt = exp(1i*(phaseX+phaseY));

            EtiltIn=Etilt.*Eshaped_input.data;
                %Ephase_field = Ephase_field.data.*BigMask;
            EtiltIn = Field(EtiltIn,(20000/8)/1024, 0.9, 'um');
            Eout_scanned = EtiltIn.propagate(Nscat,size(Nscat,3)*100);
            EscanInT(:,:,count) = EtiltIn.data;
            count=count+1;
         end
    end
    Escan = SizedArray(Escan, [(20000/8)/1024,(20000/8)/1024, 100], 'um');  
    EscanInT = SizedArray(EscanInT, [(20000/8)/1024,(20000/8)/1024, 100], 'um');
end


function Ephase_field = genPhaseStep(wavenX,wavenY,Shaped_Field,Mask_Field,CoM)
    % waven: number of waves
    
%     if size(fieldnames(Shaped_Field),1) ~ size(fieldnames(Mask_Field),1);
%         size(fieldnames(Shaped_Field),1)
%         size(fieldnames(Mask_Field),1)
%         error('The number of mask is not identical with the number of Field')
%     end
    MX = max(CoM(:,1));
    mX = min(CoM(:,1));
    
    MY = max(CoM(:,2));
    mY = min(CoM(:,2));
    
    num = size(fieldnames(Mask_Field),1);
    if ~isstruct(Shaped_Field)
        error('Shaped_Field needs to be a struct');
    end
    s = fieldnames(Shaped_Field);
    m = fieldnames(Mask_Field);
   
    
    w = getfield(Shaped_Field,s{1});
    Ephase_field=zeros(size(w,1),size(w,2));
    
    slopeX=wavenX/(MX-mX);
    slopeY=wavenY/(MY-mY);
    
    for i=1:size(CoM,1)
            phaseX = (-(wavenX/2)+slopeX*(CoM(i,1)-mX))*2*pi;
            phaseY = (-(wavenY/2)+slopeY*(CoM(i,2)-mY))*2*pi;
            phase = phaseX+phaseY;
       
            w = getfield(Shaped_Field,s{i});
            Mask = getfield(Mask_Field,m{i});
            if isobject(w)
                w = w.data*exp(1i*phase);
            else
                w = w*exp(1i*phase);
            end
            Ephase_field = Ephase_field+ w.*Mask ; 
    end
    Ephase_field = Field(Ephase_field,(20000/8)/1024, 0.9, 'um');
   
end


function Sort_Shaped_Field=Resort_Field_order(Shaped_Field)
    Sort_Shaped_Field = struct;
    s = fieldnames(Shaped_Field);
    for i = 1:44
%    tgm = angle(EscanIn(:,:,i));
%    tgm = tgm.data;
%    imagesc(tgm-rfm, [-2*pi,2*pi])
%    colorbar;
%    imagesc(angle(EscanIn(:,:,i)))
%    imagesc(EscanComb(:,:,i),[0 0.05])
%imagesc(EscanIn(:,:,i))
        img = getfield(Shaped_Field,s{i+44});
        Sort_Shaped_Field = setfield(Sort_Shaped_Field,strcat('s_',num2str(i)),img);
%imagesc(M)
%pause(0.01)
    end
    for i = 45:88
%    tgm = angle(EscanIn(:,:,i));
%    tgm = tgm.data;
%    imagesc(tgm-rfm, [-2*pi,2*pi])
%    colorbar;
%    imagesc(angle(EscanIn(:,:,i)))
%    imagesc(EscanComb(:,:,i),[0 0.05])
%imagesc(EscanIn(:,:,i))
        img = getfield(Shaped_Field,s{i-44});
        Sort_Shaped_Field = setfield(Sort_Shaped_Field,strcat('s_',num2str(i)),img);
%imagesc(M)
%pause(0.01)
    end
end


%% Calculating the biggest mask for spiral pattern
function [BigMask,Radius] = get_Spiral_Mask()
    BigMask=zeros(1144,1144);
    SingleMasks = struct;
    r = 21;
    Radius=zeros(11);

    while max(max(BigMask))<2
        r = r+1;
        for n = 1:size(CenterOfMass,1)
            x0 = CenterOfMass(n,1);
            y0 = CenterOfMass(n,2);
            img = zeros(1144,1144);
            for x = 1:1144
                for y = 1:1144
                    if sqrt((x-x0)^2+(y-y0)^2)<r
                        img(x,y) = 1;
                    else
                        img(x,y) = 0;
                    end
                end
            end
            SingleMasks = setfield(SingleMasks,strcat('s_',num2str(n)),img);
            BigMask = SumSpiral(SingleMasks);   
        end    
        clc;
        disp(r)
    end
    r=r-1;
    Radius(1)=r;
    
    for n = 1:size(CenterOfMass,1)
        x0 = CenterOfMass(n,1);
        y0 = CenterOfMass(n,2);
        img = zeros(1144,1144);
        for x = 1:1144
            for y = 1:1144
                if sqrt((x-x0)^2+(y-y0)^2)<r
                    img(x,y) = 1;
                else
                    img(x,y) = 0;
                end
            end
        end
        SingleMasks = setfield(SingleMasks,strcat('s_',num2str(n)),img);
        BigMask = SumSpiral(SingleMasks); 
    end

    for p = 10:-1:1
        while max(max(BigMask))<2
            r = r+1;
            for n = 1:size(CenterOfMass,1)
                x0 = CenterOfMass(n,1);
                y0 = CenterOfMass(n,2);
                img = zeros(1144,1144);
                if rem(n,11) < p && rem(n,11)~0;
                    for x = 1:1144
                        for y = 1:1144
                            if sqrt((x-x0)^2+(y-y0)^2)<r
                                img(x,y) = 1;
                            else
                                img(x,y) = 0;
                            end
                        end
                    end
                    SingleMasks = setfield(SingleMasks,strcat('s_',num2str(n)),img);
                else
                end
            end        
            BigMask = SumSpiral(SingleMasks);   
            clc;
            disp(r)
        end
        R(13-p)=r-1;
        r=r-1;
        disp(p);
    
    for n = 1:size(CenterOfMass,1)
        x0 = CenterOfMass(n,1);
        y0 = CenterOfMass(n,2);
        img = zeros(1144,1144);
        if rem(n,11) < p && rem(n,11)~0;
            for x = 1:1144
                for y = 1:1144
                    if sqrt((x-x0)^2+(y-y0)^2)<r
                        img(x,y) = 1;
                    else
                        img(x,y) = 0;
                    end
                end
            end
            SingleMasks = setfield(SingleMasks,strcat('s_',num2str(n)),img);    
        else
        end
        BigMask = SumSpiral(SingleMasks);
    end
end

%% Summing all the points in the spiral as a whole spiral pattern.    
function w=SumSpiral(structure)
    if ~isstruct(structure)
        error('input needs to be a structure')
    end
    s=fieldnames(structure);
    temp = getfield(structure,s{1});
    w=zeros(size(temp,1),size(temp,2));    
    for i=1:size(fieldnames(structure))    
       temp=getfield(structure,s{i});
       if isobject(temp)
            temp=temp.data;
       end
       w=w+temp;
       clear temp;
    end
end
