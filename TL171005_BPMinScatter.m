%%
% An experiment:
% The incoming field is defined in line 16,17
% The medium is defined in line 18-23
% 2 mediums are generated here, one with scattering layers (Nscattered),
% the other (N) is lens+vaccum
% Line 19 defines the scattering level of the scattering layers
tic()
focal_length = 10000;
field_size = 20000;
layer_thick = 50;
scattering_coefficient = 0.15;
GaussFilterSigma = 25;
Scat_Layer_No = 100;

format long

E = Field.plane(field_size/8,[1024 1024],0.9,'um');
E = E.aperture('gaussian',512);
Lens = E.make_lens(focal_length,layer_thick);
SN = make_medium(E,scattering_coefficient,GaussFilterSigma,Scat_Layer_No);
N = ones(size(E,1),size(E,2),2*focal_length/layer_thick);
N(:,:,(focal_length/layer_thick)+1) = Lens;
Nscattered(:,:,1:size(N,3)-Scat_Layer_No)=N(:,:,1:size(N,3)-Scat_Layer_No);
Nscattered(:,:,size(N,3)-Scat_Layer_No+1:size(N,3))=SN;

%Eout_air propagates the E through N (vaccum+lens)
%Eout_scattered propagates the E through Nscattered
tic;
E.gpu_enabled = false;
N=N(:,:,1:200);
Eout_air = E.propagate(N,layer_thick*size(N,3));
Eout_scattered = E.propagate(Nscattered,layer_thick*size(N,3));
toc;

%%
%SN_rev reverses the order of Nscattered
%A field with a point is picked from the last section of Eout_air
%This layer is Ein_point (line 42,43)
%Ein_point is then propagated through the reversed medium SN_rev
%The result is Erev_shaped
SN_rev = ones(size(Nscattered,1),size(Nscattered,2),size(Nscattered,3));
for j=0:size(Nscattered,3)-1
    SN_rev(:,:,j+1)=Nscattered(:,:,size(Nscattered,3)-j);
end


Edata = Eout_air.data;
Ein_point = Field(Edata(:,:,size(Edata,3)),field_size/(8*1024),0.9,'um');
Erev_shaped = Ein_point.propagate(SN_rev,layer_thick*size(SN_rev,3));
% The last frame of Erev_shaped (Edata) is taken and 
% its conjugated field (Edata_conj)is used to be a corrected field of the
% scattering medium Nscattered. Ein_shaped is then propagated to the Nscattered
% again, and the Eout_corrected should form a focus in the bottom layer of
% the scattering medium Nscattered
Edata = Erev_shaped.data;
Eins = abs(Edata(:,:,size(Edata,3)));
Edata = Edata(:,:,size(Edata,3))./Eins;
Edata_conj=(1./Edata).*Eins; 
Ein_shaped = Field(Edata_conj,field_size/(8*1024),0.9,'um');
Eout_corrected = Ein_shaped.propagate(Nscattered,layer_thick*size(N,3));
toc()
% Edata=zeros(1024,1024,400);
% Etemp=Erev_shaped.data;
% for i=0:399
%    Edata(:,:,i+1)=Etemp(:,:,400-i); 
% end
% Erevrev = SizedArray(Edata,10000/(8*1024) , 'um');
% clear Etemp;

