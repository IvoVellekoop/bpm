classdef FBPM < SizedArray
   properties (Access=public)
       wave_length;
       k0;
       PixelSize;
       LensMedium;
       StepSize;
   end
   methods 
        function obj = FBPM(inputField, FieldSize, lambda, unit, dz)
            PixelSize = FieldSize/size(inputField,1);
            obj = obj@SizedArray(inputField, PixelSize, unit);
            obj.PixelSize = FieldSize/size(inputField,1);
            obj.wave_length = lambda;
            obj.k0 = 2*pi/lambda;
            obj.StepSize=dz;            
        end
   end
   methods 
       function Eout = BeamPropagation(obj)
           N=LensMaker(obj,1.51,20);
           obj.LensMedium=abs(N);
           dz=obj.StepSize;
           %Eout=obj;
           
           Eout = zeros(size(obj,1), size(obj,2), size(obj.LensMedium,3));
           
           boundaries = tukeywin(size(obj, 1), 0.1) * tukeywin(size(obj, 2), 0.1).';
           %boundaries=1;
           
           FieldSize=obj.PixelSize*size(obj, 1);
           
           npts=size(obj,1);
           K=(2*pi/FieldSize)*[0:ceil(npts/2)-1 ceil(-npts/2):-1];
           [Kx,Ky]=meshgrid(K,K);
           
           Eout(:,:,1)=obj.data(:,:,1).*boundaries;
           
           for z=1:size(N,3)-1 
               clc
               disp(size(N,3))
               disp(z)
               %Nmiddle=(((N(:,:,z)+N(:,:,z+1)))/2);
               navg = mean2(N(:,:,z));
               %fE = fft2(Eout(:,:,z).*exp(1i*dz*obj.k0.*(N(:,:,z)-navg)).*boundaries);
               fE = fft2(Eout(:,:,z).*exp(1i*dz*obj.k0.*(N(:,:,z)-navg)));
               %Kx = fE.coordinates(2);
               %Ky = fE.coordinates(1);
               Kz = sqrt((navg*obj.k0)^2-Kx.^2-Ky.^2);
               Eout(:,:,z+1) = ifft2(fE.*exp(1i.*Kz*dz));
               %OutputE(:,:,z+1)= ifft2(fft2(OutputE(:,:,z)).*exp(-1i.*Kz.*dz));               
 
           end
            Eout = FBPM(Eout, FieldSize, obj.wave_length,'mm',obj.StepSize);
            Eout.LensMedium = obj.LensMedium;
       end
       
       function LensMedium = LensMaker(obj,n,f)
                HalfField = size(obj.data,1)/2;
                x = [ceil(-HalfField):1:ceil(HalfField)-1].*obj.PixelSize;
                r = sqrt((x.^2)+(x'.^2));
                k = obj.k0;
                D = max(max(r));
                T=(f-sqrt(f^2-D^2))/(n-1);
                R=(D^2+T^2)/(2*T);                
                phi=(1-n)*k.*(sqrt(R^2-r.^2)-sqrt(R^2-D^2));
                z = obj.StepSize;
                Lens = 1-((1-n)*(sqrt(R^2-r.^2)-sqrt(R^2-D^2)))/z;
                LensMedium = ones(size(obj.data,1),size(obj.data,2),1+round(f/z)*3);
                LensMedium(:,:,1+round(f/z))=Lens;
       end
       
   end
   
end