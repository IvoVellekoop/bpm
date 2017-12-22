%% System parameters: 
focal_length = 10000; % of the lens
lambda = 0.9; % um
pixel_size = lambda/2;
field_size = 3000; %um
DimSize = floor(field_size/pixel_size);
layer_thick = 100; % Thickness of every section in Z (µm)
scattering_coefficient = 2.5; % Here and the next row defines how scatter the scattering medium is. 
GaussFilterSigma = 25;
Scat_Layer_No = 3; % Layer number of the scattering medium
Additional_layers = 0; % number of additional layers of vaccum after the scattering medium


format long

%% Generate medium, lens...
% Creates field
E = Field.plane(field_size,[DimSize DimSize],lambda,'um');
E = E.aperture('gaussian',512);
% make lens
Lens = E.make_lens(focal_length,layer_thick);
% make scattering medium
SN = make_medium(E,scattering_coefficient,GaussFilterSigma,Scat_Layer_No);


% Reverse the order of the whole medium for back propagation
SN_rev = ones(size(SN,1),size(SN,2),Scat_Layer_No);
for j=0:size(SN,3)-1
    SN_rev(:,:,j+1)=SN(:,:,size(SN,3)-j);
end

%% 1) Calculating the center of mass of each focus in the spiral.
%  2) convolve each focus to a intensity profile (gaussian)
%  prerequisite file: Spiral,SpiralROI, Epoint (Ask Tzu-Lun)
%  Warning:  Because of the convolution, the dimension of the output matrix will be bigger.
s=fieldnames(SpiralROI);
CenterOfMass=zeros(88,2);
threshold = 0.015;
addpixels = 0;
q = 1;
for i=1:size(s,1)
    if rem(i,12)==1
    else
        % Calculate the position of each spiral unit by a 1024x1024 spiral
        % patterns.
%       p=strcat('F_',num2str(i));
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
        Sx=round(Sx*(1/n));
        Sy=round(Sy*(1/n));
        % Scale the position of each spiral unit to the position to the
        % current field size
        Sx=round((Sx-511.5)*(DimSize/1024)+(DimSize-1)/2);
        Sy=round((Sy-511.5)*(DimSize/1024)+(DimSize-1)/2);
        CenterOfMass(q,1)=Sx+addpixels;
        CenterOfMass(q,2)=Sy+addpixels;
%       temp=zeros(DimSize+addpixels,DimSize+addpixels);
%       temp(CenterOfMass(q,1),CenterOfMass(q,2))=1;
        q = q+1;
%       temp=conv2(temp,abs(Epoint));
%        str=sprintf('SingleField.%s=temp;',p);
%       eval(str);
    end
end
clear i j Sx Sy ip jp n q;


 
%% Propagates every focus in Nair to acquire corresponding plane wave
%  then back propagate the resulting plane wave through a reversed 
%  scattering medium (SN_revBig). Pick the final layer Eshaped_in and
%  and calculate its phase conjugation as the new shaped input for 
%  the forward scattering medium (NscatBig).  All shaped inputs are stored
%  in a structure called shaped_field
tic();
Shaped_Field = struct;

for i=1:size(CenterOfMass,1)
    disp(i)
    p=strcat('F_',num2str(i));
    temp=zeros(DimSize+addpixels,DimSize+addpixels);
    temp(CenterOfMass(i,1),CenterOfMass(i,2))=1;
    temp=conv2(temp,abs(Epoint));
    temp = temp( round((size(temp,1)-DimSize)/2)+1:round((size(temp,1)+DimSize)/2),round((size(temp,2)-DimSize)/2)+1:round((size(temp,2)+DimSize)/2));
    temp=fftshift(fft2(temp));
    temp=conj(temp);
    Erev_in= Field(temp, pixel_size, lambda, 'um');
    clear temp;
    
    F.gpu_enabled = false;
    tic()
    Erev_out = Erev_in.propagate(SN_rev,layer_thick*size(SN_rev,3));
    clear Erev_in;
    
    Erev_out = Erev_out(:,:,size(Erev_out,3));
        
    Edata = Erev_out.data;

    Erev_conj=conj(Edata);   
    Erev_back = Field(Erev_conj,pixel_size, lambda, 'um');
    clear Edata Erev_conj;   

    Eft_in = Erev_back.propagate(ones(size(Erev_back,1),size(Erev_back,2),2),layer_thick*size(SN,3));
    Eshaped_in = Eft_in.data(:,:,2);
    clear Eft_in Erev_back;
    Eshaped_in = fftshift(fft2(Eshaped_in));
    Eshaped_in = Field(Eshaped_in,pixel_size,lambda,'um');
    Shaped_Field = setfield(Shaped_Field,p,Eshaped_in);
    toc()
end

toc()

