% unit in meter
%--------------Define propagating source (2D)---------------------%
%E=csvread('field_1.csv');
%E=bessel_test;
%npts = size(E,1); % number of data points in one dimension, field is a square 


FieldSize = 5e-3; % the field size in reality
wavelength= 900e-9;
k0=2*pi/wavelength; % wavenumber of the beam
%npts = size(E,1); % number of data points in one dimension, field is a square
npts = round((FieldSize/wavelength)/10)*2;   % number of data points in x(or y) dimension, field is a square
Xe=linspace(-FieldSize/2,FieldSize/2,npts); %space scale in 1D for input field
[X,Y]=meshgrid(Xe,Xe);
nslices = 100; % number of propagation steps
dz = 0.05/500;  % StepSize along z : 100 µm/step 

E=exp(1i*((2*pi/wavelength)*0.*X+(2*pi/wavelength)*0.*Y)); %E: input field, plane wave here
input = gausswin(1112,1)*gausswin(1112,1)';
E=E.*input;

%------------------- preparing k-space---------------------------%
kmesh=2*pi/FieldSize*[0:floor(npts/2) -1:-ceil(npts/2)+1];
%kmesh=linspace(-(2*pi/FieldSize),(2*pi/FieldSize),npts);
[Kx,Ky]=meshgrid(kmesh,kmesh);
%Kz=sqrt(k0^2-Kx.^2-Ky.^2);
Kz=(Kx.^2+Ky.^2)./(2*k0);

N = AirnLens(5e-3,abs(dz),FieldSize,npts,1.53,30);  %Creating the propagation medium air+lens+air
size(N,3)
%% Run simulation
OutputE = zeros(npts,npts);
OutputE(:,:,1)=E;
%for z=1:size(N,3)-1
for z=1:size(N,3)-1
   disp(z) 
   Nmiddle=(((N(:,:,z)+N(:,:,z+1)))/2);
   OutputE(:,:,z+1)= ifft2(fft2(OutputE(:,:,z).*exp(-1i*dz*k0.*Nmiddle)).*exp(-1i.*Kz*dz));
   %OutputE(:,:,z+1)= ifft2(fft2(OutputE(:,:,z)).*exp(-1i.*Kz.*dz));
end

Eamp=abs(OutputE);

clear Nmiddle



