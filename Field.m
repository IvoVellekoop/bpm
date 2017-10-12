classdef Field < SizedArray
    properties
        lambda  % wavelength (in micrometers)
        k0      % 2 pi / lambda
        gpu_enabled = false;
    end
    methods 
        function obj = Field(E, pixel_size, lambda, unit)
            %% Construct a new field object
            % E      the actual field
            % dxy    size of the pixels in the field matrix (in 'unit's),
            % or cell array of coordinates
            % lambda wavelength (in 'unit's)
            %
            % tip: since this function works with fft, it is fastest to
            % have the size of E a power of 2, (or a power of 2 times a
            % power of 3 and/or 5)
            validateattributes(E, {'single', 'double'}, {'2d'});
            validateattributes(lambda, {'single', 'double'}, {'scalar', 'positive'});
            obj = obj@SizedArray(E, pixel_size, unit);
            obj.lambda = lambda;
            obj.k0  = 2 * pi / obj.lambda;
        end
        function P = power(obj)
            %% Returns the total power in the field, per unit or per unit^2
            % This function simply integrates |E|.^2 over the full field
            P = obj.data(:)' * obj.data(:);
            P = P * prod(obj.pitches); %take into account the surface area of each pixel (dx dy in the integral)
        end
        function Eout = propagate(obj, n, z)
            %% Propagate the field through refractive index map n
            % Propagates the wave over a total distance z
            % 
            % This function first scatters the light (using a slice of the
            % refractive index map n) and then propagates the light
            % (using a k-vector that is based on the average refractive
            % index of the total slice. This process is repeated for each
            % slice in n.
            % 
            % If n is a scalar, this function simply propagates over
            % distance z without scattering
            %
            % Note that the field will wrap around at the edges due to the
            % fft. Proper boundaries need to be defined.
            %
            % n: 3-dimensional array (y,x,z), or scalar
            % returns: field at each propagation step (in a 3-D matrix)
            validateattributes(z, {'single', 'double'}, {'scalar'});
            
            %% Create absorbing boundaries (todo: optimize & allow tuning)

            boundaries = tukeywin(size(obj, 1), 0.1) * tukeywin(size(obj, 2), 0.1).';
            
            
            %% Setup output array and loop variables
            Nslices = size(n, 3);
            dz = z / Nslices;
            Eout = zeros(size(obj,1), size(obj,2), Nslices);
            E = obj;
            
            %% Standard beam propagation loop: 
            % diffract, Fourier tranform, propagate, transform back, repeat...
            % We base the propagation step on the average refractive index
            % for that slice, so that the inhomogeneous part of the refractive
            % index is minimized, and the accuracy is optimal.
            %
            
            if obj.gpu_enabled
                obj.data = gpuArray(obj.data);
                obj.k0 = gpuArray(obj.k0);
                dz = gpuArray(dz);
                n = gpuArray(n);
                boundaries = gpuArray(boundaries);
            end
               
            for s=1:Nslices
                slice = n(:,:,s);
                navg = mean2(slice);
                fE = fft2(E .* (exp(1.0i * dz * obj.k0 * (slice-navg)) .* boundaries)); % scatter and Fourier transform field
                kx = fE.coordinates(2);
                ky = fE.coordinates(1);
                kz = sqrt((navg * obj.k0)^2 - ky.^2.' - kx.^2);                
                E = ifft2(fE .* exp(1.0i * dz * kz)); %propagate field and inverse Fourier transform
                Eout(:, :, s) = E.data;
            end
            Eout = gather(Eout);
            Eout = SizedArray(Eout, [obj.pitches, dz], obj.unit(1)); 
        end
        
        function Eout = lens(obj, focal_length)
            %% Applies a lens function to the field
            % When starting with a plane wave, the resulting wave is
            % perfectly spherical and converges at the focal_length (which
            % may be negative)
            validateattributes(focal_length, {'single', 'double'}, {'scalar'});
            x = obj.coordinates(2);
            y = obj.coordinates(1);
            extra_path = sqrt(focal_length^2 + y.'.^2 + x.^2)-abs(focal_length);
            if (focal_length < 0)
                extra_path = -extra_path;
            end
            Eout = obj .* exp(-1.0i * obj.k0 * extra_path);
        end
        function Eout = aperture(obj, type, r, center)
            %% Applies an intensity mask to the field, mimicking e. g. a circular apertur
            % obj       field to apply the intensity mask to
            % type      type of mask, can currently be 'circular' or 'gaussian'
            % r         'radius' of the mask. Exact meaning depends on aperture type
            % center    center position of the mask, specified in the same 
            %           units as the field coordinates. Defaults to 0.0
            %
            validateattributes(r, {'single', 'double'}, {'scalar', 'positive'});
            validateattributes(type, {'char'}, {});
            if nargin < 4
                center = [0,0];
            end
            %calculate coordinates relative to the center, scaled with
            %'radius'
            x = (obj.coordinates(2) - center(2))/r;
            y = (obj.coordinates(1) - center(1))/r;
            
            switch type
            case 'circular'
                mask = x.^2 + y.'.^2 < 1.0;
            case 'gaussian'
                mask = exp(-(x.^2 + y.'.^2));
            end
            Eout = obj .* mask;
        end
        
        function Lens = make_lens(obj,focal_length,dz)
            %% This procedure produces a GRIN lens with the focal_length specified in focal_length.
            %% The output, Lens, is the gradient refractive index map required to 
            %% form focus as focal_length specified in a meterial with thickness dz 
            %% dz   physical thickness of a layer (step_size) in the whole medium           
            x = obj.coordinates(2);
            y = obj.coordinates(1);          
            extra_path = abs(focal_length)-sqrt(focal_length^2 - y.'.^2 - x.^2);
            if (max(max(abs(diff(extra_path)))) > obj.lambda / 2)
                disp(round(max(max(abs(diff(extra_path))))/obj.lambda))
                error('pixel number is not high enough to sample the lens curvature');
            end
            if (focal_length < 0)
                extra_path = -extra_path;
            end
            Lens = ((max(max(extra_path))-extra_path)/dz)+1;
        end
        function Scatter_Medium = make_medium(obj,scattering_coef,anisotropy_f,Layer_N)
             Scatter_Medium = scattering_coef*rand(size(obj,1),size(obj,2),Layer_N)+1;
             for i = 1:Layer_N
                Scatter_Medium(:,:,i) = imgaussfilt(Scatter_Medium(:,:,i),anisotropy_f);
             end              
        end
        %% Calculating the phase conjugate of the field
        function Eout=PhaseConjugate(obj)
            Edata = obj.data;
            Eins = abs(Edata(:,:,size(Edata,3)));
            Edata = Edata(:,:,size(Edata,3))./Eins;
            Eshaped_conj=(1./Edata).*Eins;
            %%%%%%%%%%%%%%%%%%%%%%%% 20000 %%%%%%%%%%%%%%%%%%%%%%%%%%%
            Eout = Field(Eshaped_conj, (20000/8)/1024, 0.9, 'um');
        end
        
        %% This step create phase step across the field according to the position
        %  of each focus.
        function Ephase_field = genPhaseStep(waven,Shaped_Field,CoM,dim)
            % waven: 
            % the tilting distance (D) in Z measured by number of wave, 
            % e.g. waven = 50, 50 waves to reach the distance as D     
            % Shaped_Field: a structure contains the field of each focus
            % CoM: Center of mass of each focus
            % dim: desired dimension to apply phase step.
            if ~isstruct(Shaped_Field)
                error('Shaped_Field needs to be a struct');
            end            
            
            M = max(CoM(:,dim));
            m = min(CoM(:,dim));
            
            %%%%%%%%%%%%%%%% 1144 %%%%%%%%%%%%%%%%%%%%%
            Ephase_field=zeros(1144,1144);    
            slope=waven/(M-m);

            
            s = fieldnames(Shaped_Field);
            
            for i=1:size(CoM,dim)
                phase = (-(waven/2)+slope*(CoM(i,dim)-m))*2*pi;
                w = getfield(Shaped_Field,s{i});
                w = w.data*exp(1i*phase);
                Ephase_field = Ephase_field + w; 
            end
            %%%%%%%%%%%%%%% 20000 %%%%%%%%%%%%%%%%%%%%
            Ephase_field = Field(Ephase_field,(20000/8)/1024, 0.9, 'um');
        end
    end
    
    methods (Static)
        function Eout = plane(dimensions, subdivs, wavelength, unit)
            %% Generates a plane wave
            %
            % The amplitude of the plane wave is normalized to have unit
            % power per unit. So, regardless of the dimensions and number
            % of subdivisions, power(Field) = 1 (per unit or per unit^2)
            %
            % dimensions: size of field (in units). Can be 1-D or 2-D
            % subdivs:  number of subdivisions in each dimension. Can be
            %           scalar (for same pixels size in both dimension)
            % wavelength: wavelenght of the light
            % unit:     unit for dimensions and wavelength
            %
            if (isscalar(subdivs)) %repeat same value if subdivs is a scalar
                subdivs = ones(1, length(dimensions)) * subdivs;
            end
            Eout = Field(ones(subdivs), dimensions./subdivs, wavelength, unit); 
            Eout = Eout / sqrt(power(Eout));
        end
        function test()
            f = 10000; %um
            D = 1000; %um
            E = Field.plane([2*D, 2*D], 512, 0.6328, 'um');
            %E = Field(ones(256,512), 0.2, 0.6328, distance_unit);
            E = lens(E, f); %focus at a distance of 1000 mu
            E = aperture(E, 'gaussian', D/2);
            figure(1); imagesc(angle(E)); %display phase of the incident light (converging wave)
            tic();
            E3D = propagate(E, ones(1,1,32), 2*f); %propagate to the focus through air in 32 steps
            toc();
            figure(2); imagesc(E3D(end/2, :, :)); %cross section in x-z plane
            figure(3); imagesc(E3D(:,:,end/2)); %cross section in focal plane (xy)
        end
    end
end