classdef VCAO_simplified
    %% VCAO was used to simulate the idea of volume-conjugated adaptive optics.
    % In this procedure, the forward direction means from back focal plane
    % to lens to scattering medium, and the reverse direction is the other
    % way.
    % The procedure is executed as following:  
    % 1) The system parameters was defined mostly in VCAO(...) and the
    % scattering medium is generated by CreateMedium(...)
    % 2) Need a pattern on the back focal plane of the objective, such as
    % a spiral pattern. The position of each point in the spiral pattern
    % will be calculated by CalcCoM(...) and stored
    % 3) Each point in the spiral pattern will be propagated in the 
    % scattering medium individually in GetScattered(), and Fourier
    % transformed to what it should be in the back focal plane again, which
    % is the output.
    % 4) The SumField() will added the output of all spiral points into a
    % summed field as a new input.
    % 5) revPropagate() will propagate the field in the backfocal plane by
    % first: inverse fourier transformed the field, and second: propagate
    % in vaccum back to the position of the front surface of the scattering 
    % medium, and third: propagate in the forward direction through the
    % medium.
    % 6) StepPhaseScanning() performs the step-phase scanning on the summed
    % field created by step 4).
    % 7) GradPhaseScanning() put a phase gradient instead of phase step on
    % top of the summed field in step 4) to check the memory effect.
    % 8) finding_usp(): This function is used to calculate the ltr, and the
    % size of memory effect range.
    properties
        options;
        medium;
    end
    
    methods 
        function obj = VCAO_simplified(options)
        %% Prepares the medium and grid for the simulation
        % options contains the following fields:
        % options.lambda = wavelength (in micrometers)
        % options.size = size (y, x) of simulation (in pixels). (Prefer
        % powers of 2, 3 or 5 for better performance).
        % options.pixel_size = size of a single pixel
        % options.medium_thickness = thickness of medium (in micrometers)
        % options.layer_count = number of scattering layers in the medium
        % options.scattering_strength = number that determines scattering
        % properties of the medium (together with scattering_filter).
        % Higher values correspond to stronger scattering.
        % options.scattering_filter = higher values correspond to more
        % forward scattering
            obj.options = options;
            obj.medium = VCAO_simplified.generate_medium(obj.options);
        end
        
        function Eout=bfp_to_fp(obj, Ein, scattering)
            % Propagates a field from the back focal plane to the front
            % focal plane. If 'scattering' is set to false, this function
            % only performs a Fourier tranform. If 'scattering' is set to
            % true, the functions mimicks a lens with a scattering medium
            % at some distance behind it. This is done by first Fourier
            % tranforming (propagation without scattering to the focal
            % plane), then propagating back using bpm (again without
            % scattering), and finally propagating forwards through the
            % scattering medium until the focal plane is reached.
            % This procedure avoids the time-consuming process of
            % simulating the lens and the propagation between the lens and
            % the medium.
            % Currently Ein is a plane matrix without units. We need
            % to implement the fact that it may have a different size than
            % the simulation, and incorporate the effect of the focal
            % length, so Ein should be a Field (with size of back aperture).
            %
            Eout = Field(fftshift(fft2(fftshift(Ein))), obj.options.pixel_size, obj.options.lambda, 'um');
            if scattering
                Eout = conj(Eout); %reverse direction
                Eout = Eout.propagate(1, obj.options.medium_thickness); %back propagate to sample surface
                Eout = conj(Eout); %forward direction again
                Eout = Eout.propagate(obj.medium, obj.options.medium_thickness); 
            end
        end
        
        function Eout=fp_to_bfp(obj, Ein, scattering)
            % See bfp_to_fp, but other way around. 
            % As a convention, we consider all fields to be pointed in the
            % positive z direction (from bfp to fp). Since we are
            % propagating in the reverse direction here, we need several modifications:
            % first the field is conjugated to reverse the direction
            % then it is propagated through the medium in reverse order 
            % (the medium is flipped upside down), and at the bft it is 
            % phase conjugated again.
            % this way we can do: fp_to_bfp(bfp_to_fp(source)) to get back
            % (approximately) the same field.
            if scattering
                Eout = conj(Ein); %phase conjugate field, we are now propagating to negative z
                Eout = Eout.propagate(flip(obj.medium, 3), obj.options.medium_thickness); %propagate through medium in reverse
                Eout = conj(Eout); %forward direction again
                Eout = Eout.propagate(1, obj.options.medium_thickness); %propagate from sample surface to focal plane
            else
                Eout = Ein;
            end
            Eout = fftshift(ifft2(fftshift(Eout.data))); %propagate from focal plane to back focal plane
        end
        
        function shaped_fields = get_optimized_wavefronts(obj, centers, spot_size)
        %% Each position recorded in the obj.CofM (position of foci on the
        % back focal plane) will be assinged to a gaussian PSF according to
        % the SpotSize specified. Then these points are Fourier transformed 
        % to be plane waves on the front focal plane one by one. The plane 
        % wave is then reverse propagated through the scattering medium,
        % and forward propagated through the vaccum. The output will then
        % be Fourior transformed to the field appeared in the front focal
        % plane.
            Ndots = size(centers, 2);
            shaped_fields = cell(Ndots, 1);

            for i=1:Ndots
                disp([i, Ndots]);
                % generate 'dot' with Gaussian envelope at center of mass
                % (in the back focal plane)
                x = -obj.options.size(2)/2:obj.options.size(2)/2-1;
                y = (-obj.options.size(1)/2:obj.options.size(1)/2-1).';
                Ebfp = exp(-((x-centers(1,i)).^2+(y-centers(2,i)).^2) / spot_size^2);
                
                % propagate without scattering to get the plane wave
                Eplane_wave = bfp_to_fp(obj, Ebfp, false); 
                
                % propagate plane wave back through scattering medium and
                % phase conjugate to get the ideal wavefront for creating
                % that plane wave
                shaped_fields{i} = fp_to_bfp(obj, Eplane_wave, true);
            end    
        end
        
        function [Escan,EscanIn] = StepPhaseScanning(obj,MaxXw,MaxYw,Shaped_Field)
        %% This function performs step phase scannig base on the field
        % stored in obj.OriSpiral. The center of mass stored in obj.CofM is
        % required to impose phase differences on differenct foci. MaxXw
        % and MaxYw specify the maximum phase differences in X and Y
        % dimensions, and the unit is wave. The step size is 5 waves (see
        % line 232 233).
            count = 1;
            TotalThick = obj.layer_thick*obj.Scat_Layer_No;
            for i=1:44
                CenterOfMass(i,:) = obj.CofM(i+44,:);
            end
            for i=45:58
                CenterOfMass(i,:) = obj.CofM(i-44,:);
            end           
            dim = ((MaxXw*2/5)+1)*((MaxYw*2/5)+1);
            Escan = zeros(obj.DimSize,obj.DimSize,dim);
            EscanIn = zeros(obj.DimSize,obj.DimSize,dim);
    
            for wavenX = MaxYw:-5:-MaxYw
                for wavenY = -MaxXw:5:MaxXw
                    tic()
                    Ephase_field = obj.genPhaseStep(wavenX,wavenY,Shaped_Field,CenterOfMass);
                    %Ephase_field = Ephase_field.data.*BigMask;
                    %Ephase_field = Field(Ephase_field.data,obj.pixel_size, obj.lambda, obj.unit);
                    %Eout_scanned = Ephase_field.propagate(obj.medium,size(Nscat,3)*obj.layer_thick);
                    Eout_scanned = obj.revPropagate(Ephase_field.data);
                    EscanIn(:,:,count) = Ephase_field.data;
                    Escan(:,:,count) = Eout_scanned.data;
                    toc()
                    count=count+1;
                end
            end
    
            Escan = SizedArray(Escan, [obj.pixel_size,obj.pixel_size, obj.layer_thick], obj.unit);
            EscanIn = SizedArray(EscanIn, [obj.pixel_size,obj.pixel_size, obj.layer_thick], obj.unit);
        end
        
        function Ephase_field = genPhaseStep(obj,wavenX,wavenY,Shaped_Field,CoM)
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
            
%            num = size(fieldnames(Mask_Field),1);
%             if ~isstruct(Shaped_Field)
%                 error('Shaped_Field needs to be a struct');
%             end
            s = fieldnames(Shaped_Field);
%             m = fieldnames(Mask_Field);
            
            
            w = getfield(Shaped_Field,s{1});
            Ephase_field=zeros(size(w,1),size(w,2));
            
            slopeX=wavenX/(MX-mX);
            slopeY=wavenY/(MY-mY);
            
            for i=1:size(CoM,1)
                phaseX = (-(wavenX/2)+slopeX*(CoM(i,1)-mX))*2*pi;
                phaseY = (-(wavenY/2)+slopeY*(CoM(i,2)-mY))*2*pi;
                phase = phaseX+phaseY;
                
                w = getfield(Shaped_Field,s{i});
 %               Mask = getfield(Mask_Field,m{i});
                if isobject(w)
                    w = w.data*exp(1i*phase);
                else
                    w = w*exp(1i*phase);
                end
 %               Ephase_field = Ephase_field+ w.*Mask ;
                Ephase_field = Ephase_field+w;
            end
            Ephase_field = Field(Ephase_field,obj.pixel_size, obj.lambda, obj.unit);           
        end
   
        function [Escan,EscanInT] = GradPhaseScanning(obj,MaxXw,MaxYw)
            count = 1;
            dim = ((MaxXw*2/5)+1)*((MaxYw*2/5)+1)
            Escan = zeros(obj.DimSize,obj.DimSize,dim);
            EscanInT = zeros(obj.DimSize,obj.DimSize,dim);
            for wavenX = MaxXw:-5:-MaxXw
                for wavenY = -MaxYw:5:MaxYw
                    phaseX = linspace(-2*pi*wavenX/2,2*pi*wavenX/2,obj.DimSize);
                    phaseX = phaseX'.*ones(obj.DimSize,1);
                    phaseY = linspace(-2*pi*wavenY/2,2*pi*wavenY/2,obj.DimSize);
                    phaseY = ones(obj.DimSize,1)'.*phaseY;
                    Etilt = exp(1i*(phaseX+phaseY));
                    
                    EtiltIn=Etilt.*obj.OriSpiral.data;
                    %Ephase_field = Ephase_field.data.*BigMask;
                    %EtiltIn = Field(EtiltIn,obj.pixel_size, obj.lambda, obj.unit);
                    Eout_scanned = obj.revPropagate(EtiltIn);
                    EscanInT(:,:,count) = EtiltIn;
                    Escan(:,:,count) = Eout_scanned.data;
                    count=count+1
                end
            end
            Escan = SizedArray(Escan, [obj.pixel_size,obj.pixel_size, obj.layer_thick], obj.unit);
            EscanInT = SizedArray(EscanInT, [obj.pixel_size,obj.pixel_size, obj.layer_thick], obj.unit);
        end
 
        function obj=find_musp(obj)
            %DeltaX = 0;
            %DeltaY = 0;
            %Xp=linspace(-DeltaX*pi,DeltaX*pi,obj.DimSize);
            %Yp=linspace(-DeltaY*pi,DeltaY*pi,obj.DimSize);
            
            F = ones(obj.DimSize,obj.DimSize);
            %for i = 1:DimSize
            %    for j = 1:DimSize
            %        F(i,j)=Xp(i)+Yp(j);
            %    end
            %end
            
            Etilt=exp(1i*F);
            Etilt=Field(Etilt,obj.lambda/2, obj.lambda, obj.unit);
            Etilt.gpu_enabled = false;
            
            Etilt_out = Etilt.propagate(obj.medium,obj.Scat_Layer_No*obj.layer_thick,-1);
            %subplot(3,3,5);
            %imagesc(Etilt_out(:,:,Scat_Layer_No))
            
            FFT_Etout=fftshift(fft2(Etilt_out.data(:,:,1)));
            %subplot(3,3,6);
            figure;imagesc(abs(FFT_Etout))
            
            Pkb = abs(FFT_Etout((end)/2,:));
            
            
            L = obj.Scat_Layer_No*obj.layer_thick; %Thickness,�m%
            k0 = 2*pi/obj.lambda;
            
            N = obj.DimSize;
            dx = obj.lambda/2;
            x = (0:N-1)*dx;
            dk = 2*pi/(dx*N);
            k = ((0:N-1) - N/2)*dk;
            f = fit(k',Pkb','gauss1');
            c =  f.c1  ;
            
            ltr = ((k0^2)*L/(2*c^2));
            sqrt(6*ltr/((k0^2)*L))
            obj.ltr = ltr;
        end
        
        function centers = spiral_coordinates(obj,SpiralROI,Spiral)
        % This function is used to calculate the center of mass of each
        % point on the spiral pattern, and scale them to the current field
        % size. The spiral pattern SpiralROI and Spiral are prerequisite in
        % this function. For different patterns, a new function is needed.
            s=fieldnames(SpiralROI);
            CenterOfMass=zeros(88,2);
            threshold = 0.015;
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
                    Sx=round((Sx-511.5)*(obj.DimSize/1024)+(obj.DimSize-1)/2);
                    Sy=round((Sy-511.5)*(obj.DimSize/1024)+(obj.DimSize-1)/2);
                    CenterOfMass(q,1)=Sx;
                    CenterOfMass(q,2)=Sy;
                    q = q+1;
                    obj.CofM = CenterOfMass;
                end
            end
        end   
    end
    
    methods (Static)
        function medium = generate_medium(options)
%% Creates a refractive index matrix with the given scattering
% strength, anisotropy, and number of layers.
% The resulting medium has an average refractive index of n
% and the standard deviation of the refractive index fluctuations
% is given by options.scattering_strength.
% 
% The size of the fluctuations is given by options.scattering_filter.
% which is the correlation length in um.
% The two parameters can be varied independently, where
% options.scattering_strength correlates with the scattering coefficient
% and options.scattering_filter correlates with 1-g
%
% todo: return a generator object instead to save memory.
% todo: remove medium generation from all other places!
%
% Implementation: this function first generates a randomly distributed
% value for each pixel. The resulting medium is filtered with a Gaussian
% blur. The blurring will reduce the standard deviation of the values by
% a factor of sqrt(N), with N the effective size (in number of pixels) of
% the Gaussian kernel. todo: compensation is not 100% accurate

            % Step 1: calculate correlation length in pixels and calculate
            % scattering strength compensated for blur
            correlation_length = options.scattering_filter / options.pixel_size;
            kernel_size = max(pi*correlation_length, 1);
            scattering_strength = options.scattering_strength * kernel_size;

            % Step 2: Generate white noise medium
            medium = 1 + scattering_strength * randn(options.size(1), options.size(2), options.layer_count);

            % Step 3: low pass filter
            for n = 1:size(medium, 3)
                medium(:, :, n) = imgaussfilt(medium(:, :, n), correlation_length, 'FilterDomain', 'frequency', 'Padding', 'circular');
            end
        end

        function centers = ring_coordinates(options)
%% This function returns coordinates corresponding to dots on a series of rings
% the dots are positioned to have approximately even spacing.
% options.pupil_radius = radius of the ring (same units as the returned centers)
% options.ring_count = number of rings. The rings are evenly spaced. The
% number of dots per ring is adjusted to have approximately even spacint
            
            % step 1: calculate radius and area for each of the rings
            % (include an extra ring on the outside so that we can easily
            % calculate the areas)
            r = linspace(0, options.pupil_radius, options.ring_count);
            dots_per_ring = (1:options.ring_count).^2;
            
            centers = zeros(2, 0);
            offset = 0; %starting angle
            for n=1:options.ring_count
                phases = (1:dots_per_ring(n)) / dots_per_ring(n) * 2 * pi + offset;
                xcoords = r(n)*cos(phases);
                ycoords = r(n)*sin(phases);
                coords = [xcoords; ycoords];
                centers = [centers, coords]; 
                offset = phases(1)/2;
            end
        end

        function Shaped_Field = PropagateASpiral(field_size,lambda,unit,Scat_Layer_No,scat,anisotrop,SpiralROI,Spiral)
            obj = VCAO(field_size,lambda,unit);
            obj = obj.CreateMedium(Scat_Layer_No,scat,anisotrop);
            obj = obj.CalcCoM(SpiralROI,Spiral);
            Shaped_Field=GetScattered(obj,1.5);
        end    
    end
end




 




