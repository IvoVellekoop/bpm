classdef Field < SizedArray
    properties
        lambda  % wavelength (in micrometers)
        k0      % 2 pi / lambda
    end
    methods
        function obj = Field(E, pixel_size, lambda, unit, varargin)
        %% Construct a new field object from existing data
        %   E             values of the electric field
        %   pixel_size    size of the pixels in the field matrix (in 'unit's),
        %                 can be a vector to indicate rectangular pixels
        %   lambda        wavelength (in 'unit's)
        %
        % tip: since this function works with fft, it is fastest to
        % have the size of E a number that is a product of powers of 2, 3 and 5.
        %
            validateattributes(E, {'single', 'double'}, {'2d'});
            validateattributes(lambda, {'single', 'double'}, {'scalar', 'positive'});
            obj = obj@SizedArray(E, pixel_size, unit, varargin{:});
            obj.lambda = lambda;
            obj.k0  = 2 * pi / obj.lambda;
        end
        
        function P = power(obj)
            %% Returns the total power in the field, per unit or per unit^2
            % This function simply integrates |E|.^2 over the full field
            P = obj.data(:)' * obj.data(:);
            P = P * prod(obj.pitches); %take into account the surface area of each pixel (dx dy in the integral)
        end
        
        function [Eout, Elayers] = propagate(obj, n, total_distance, showplot)
%% Propagate the field through refractive index map n
%   Propagates the wave over a total distance total_distance through a
%   structure with refractive index n.
%
%   n is a 3-D matrix with refractive index values. Each slice n(:,:,i)
%   corresponds to one scattering layer. For each layer, this function
%   first propagates the light by dz/2 (with dz the thickness of the layer)
%   then scatters the light, and propagates an additional dz/2.
%
%   Eout is the electric field after propagating through all layers
%   optionally, one can specify Elayers to additionally get the field
%   in the exact center of each layer (directly after scattering).
%
%   If n is a scalar, this function simply propagates over
%   distance z without scattering
%
%   Note that the field will wrap around at the edges due to the fft. To
%   reduce this effect, absorbing boundaries are added, but this will only
%   work when the step size is small with respect to divergence of the
%   simulated field.
%
%   Note: the plot can be surpressed by passing false to the showplot argument. By
%   default this is true.
            validateattributes(total_distance, {'single', 'double'}, {'scalar'});
            store_all = nargout > 1;
            
            if nargin < 4
                showplot = false;
            end

            %% Create absorbing boundaries (todo: optimize & allow tuning)
            boundaries = tukeywin(size(obj, 1), 0.1) * tukeywin(size(obj, 2), 0.1).';

            %% Setup output array and loop variables
            Nslices = size(n, 3);
            dz = total_distance / Nslices;
            if store_all
                Elayers = zeros(size(obj,1), size(obj,2), Nslices);
            end

            %% Standard beam propagation loop:
            % diffract, Fourier tranform, propagate, transform back, repeat...
            % We base the propagation step on the average refractive index
            % for that slice, so that the inhomogeneous part of the refractive
            % index is minimized, and the accuracy is optimal.
            %
            % note: we don't need to explicitly convert to gpuArray since
            % the gpuArray property is contageous: if the original data
            % passed to create the field was a gpuArray, it will remain a
            % gpuArray and all other arrays will be converted automatically.
            %

            % start with Fourier transformed field (apply boundaries)
            fE = fft2(obj .* boundaries);
            kx2 = fE.coordinates(2).^2;
            ky2 = fE.coordinates(1).'.^2;
            fE = fE.data; %workaround for performance problem
            
            % propagate 1/2 step, apply scattering, propagate next 1/2 step
            for s=1:Nslices
                slice = obj.convert_type(n(:,:,s));
                navg = mean2(slice);
               
                % propagate half way
                ekz = exp(0.5i * dz * sqrt(complex((navg * obj.k0)^2 - ky2 - kx2)));
                fE = fE .* ekz;
                
                % scatter
                E = ifft2(fE);
                E = E .* (exp(1.0i * dz * obj.k0 * (slice-navg)) .* boundaries);
                
                if store_all
                    Elayers(:, :, s) = gather(E);
                end
                if showplot
                    imagesc(real(E)); drawnow();
                end
                fE = fft2(E);

                % propagate next half
                fE = fE .* ekz;
            end
            Eout = obj.with_data(ifft2(fE));
            if Nslices ==1
                Elayers(:,:,2) = gather(Eout.data);
            end
            if store_all
                opt.EnableGPU = false;
                Elayers = SizedArray(Elayers, [obj.pitches abs(dz)], [obj.units, obj.units(1)], opt);
            end
        end

        function [E_out,E_layers] = backPropagate(E,n,total_distance)
          %  Back propagates the wave. It does the following:
          %     1) Change positive distance to negative distance.
          %     2) Reverse the order of scattering E_layers.
          %     3) Uses propagate() to propagate the E field with above inputs.
          %  Effectively performs the adjoint of propagate().
        warning("untested, please remove warning if this function works")
          total_distance = -total_distance;

          n_reverse = flip(n,3);

          [E_out,E_layers] = propagate(E,n_reverse,total_distance);


        end

        function Eout = tilt(obj, kx, ky)
            %% Applies a tilt kx,ky to the field
            % This function multiplies the field by e^(kx x + ky y)
            % kx and ky are given in radians/unit (with the unit specified
            % in the constructor, typically micrometers)
            % If you need to tilt the wavefront by a specified angle,
            % for example in the x-direction use the equation 
            % kx = k0 sin(alpha)
            %
            Eout = obj.* exp(1.0i * kx * obj.coordinates(2) + 1.0i * ky * obj.coordinates(1).');        
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
            warning("untested, please remove warning if this function works")

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
    end

    methods (Static)
        function Eout = plane(dimensions, subdivs, wavelength, unit, varargin)
            %% Generates a plane wave
            %
            % The amplitude of the plane wave is normalized to have unit
            % power per unit. So, regardless of the dimensions and number
            % of subdivisions, power(Field) = 1 (per unit or per unit^2)
            %
            % dimensions: size of field (in units). Can be 1-D or 2-D
            % subdivs:  number of subdivisions in each dimension. Can be
            %           scalar (in which case it is the same for both dimensions)
            % wavelength: wavelenght of the light
            % unit:     unit for dimensions and wavelength
            % options, see Field constructor
            Eout = Field(ones(subdivs), dimensions./subdivs, wavelength, unit, varargin{:});
            Eout = Eout/sqrt(power(Eout));
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
            [E, E3D] = propagate(E, ones(1,1,32), 2*f, true); %propagate to the focus through air in 32 steps
            toc();
            figure(2); imagesc(E3D(end/2 + 1, :, :)); %cross section in x-z plane
            figure(3); imagesc(E3D(:,:,end/2 + 1)); %cross section in focal plane (xy)
        end
    end
end
