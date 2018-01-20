%% System parameters: 
options.lambda = 0.6; %um (inside medium)
options.size = [2048, 2048];
options.pixel_size = options.lambda/2;
options.medium_thickness = 500; %um
options.layer_count = 4;
options.scattering_strength = 0.01; %standard deviation of refractive index
options.scattering_filter = 50; %feature size in um
options.pupil_radius = 0.25*options.size(1); %radius of back focal plane (in pixels)
options.ring_count = 5;
options.spot_size = 2; %size of 'dot' in back focal plane (in pixels). Larger value gives smaller field of view (but higher contrast?)

%% Generate simulation object
disp('generating medium');
v = VCAO_simplified(options);

%% Test dot ring generation
centers = VCAO_simplified.ring_coordinates(options);
figure(3);plot(centers(1,:), centers(2,:), '.');

%% Test phase conjugation of point source
disp('testing without scattering');
source = zeros(options.size);
source(end/2+1, end/2+1) = 1;
source = imgaussfilt(source, 3);
Escat = v.bfp_to_fp(source, false); % propagate to focus: plane wave
Erefocused = v.fp_to_bfp(Escat, false); % propagate back: sharp focus
figure(1); imagesc(real(Erefocused));
drawnow;

%%
disp('testing with scattering');
Escat = v.bfp_to_fp(source, true);
Erefocused = v.fp_to_bfp(Escat, true);
figure(2); imagesc(real(Erefocused))
drawnow;

%%
disp('construct set of wavefronts that (after scattering) become plane waves in the focal plane');
centers = VCAO_simplified.ring_coordinates(options);
wavefronts = v.get_optimized_wavefronts(centers, options.spot_size);

%% show generated wavefronts
%for loop=1:10
    for f=1:length(wavefronts)
        imagesc(real(wavefronts{f})); pause(0.5);
    end
%end


%% Sum all wavefronts from experiment above and see if you get a focus
% note: the contrast will be better if we have more wavefronts (a higher
% options.ring_count)
disp('interference focus');
sumfield = 0;
for w=1:length(wavefronts)
    sumfield = sumfield + wavefronts{w};
end
Eifocus = v.bfp_to_fp(sumfield, true);
imagesc(Eifocus);


%% Scan the focus using the tilt/tilt memory effect
% if the scattering strength is high enough, there should be only a very
% limited memory effect
Efoci = cell(5,1);
Efoci{1} = Eifocus;
for shift=1:6
    phase_gradient = exp(1.0i * 8 * pi * shift * (1:options.size(1))/options.size(1));
    Efoci{shift+1} = v.bfp_to_fp(sumfield .* phase_gradient, true);
    imagesc(Efoci{shift+1});
end    

%% display 'movie' (intensity squared)
for loop=1:10
    for f=1:length(Efoci)
        imagesc(abs(Efoci{f}.data(900:1100, 900:1100)).^4); pause(0.5);
    end
end

%% VCAO step scanning
Estepfoci = cell(5,1);
Estepfoci{1} = Eifocus;
for shift=1:6
    sumfield = 0;
    for w=1:length(wavefronts)
        sumfield = sumfield + wavefronts{w} * exp(1.0i * 8 * pi * shift * centers(1, w)/options.size(1));
    end
    Estepfoci{shift+1} = v.bfp_to_fp(sumfield, true);
    imagesc(Estepfoci{shift+1});
end    

%% display 'movie' (intensity squared)
for loop=1:10
    for f=1:length(Estepfoci)
        imagesc(abs(Estepfoci{f}.data(900:1100, 900:1100)).^4); pause(0.5);
    end
end
