classdef SizedArray
    %SizedArray Array that also stores information about the coordinates
    %   Detailed explanation goes here
    
    properties
        data %raw array data
    end
    properties (Access=protected)
        pitches %vector of dx,dy,dz, etc. 
        is_fft %true indicates that we are dealing with fft-shifted coordinates 
        % when true, 'coordinates' returns negative values in the last half of the vector
        % (see coordinates). And 'unit' returns inverse units.
        units %units of x,y,z,etc. coordinates
    end
    
    methods
        function obj = SizedArray(data, pitches, units)
            %% Construct a new sized array object for the given data
            % 
            % SizedArray(data, pitches, units)
            % Pitches specifies a vector of pitches (dx, dy, dz) that
            % are used to construct standard coordinates:
            % ('almost symmetric' coordinates that include the exact point 0). 
            % If a scalar is passed, it is assumed that all dimensions have
            % the same pitch.
            % As a special case, the pitch can be imaginary. In this case,
            % fft-shifted coordinates are generated
            %
            % 'units' should be a cell array of strings with the units for
            % each dimension, or it can be an array of Unit objects
            % If only one Unit or unit string is given, it is used for all
            % dimensions
            %
            % Note: arrays of SizedArray objects are not supported!
            % 
            if nargin == 0 %required 'no argument constructor' returns trivial default object
                return
            end
            
            dims = ndims(data);
            obj.is_fft = false;
            obj.data = data;
            
            % automatically repeat pitch if only one is given
            if isscalar(pitches)
                obj.pitches = ones(1, dims) * pitches;
            else
                obj.pitches = pitches;
            end
            validateattributes(obj.pitches, {'single', 'double'}, {'vector', 'numel', dims});

            % parse 'units'
            us(dims) = Unit(); %pre-allocate
            if isa(units, 'Unit')
                us(1:dims) = units; % automatically repeats unit if only one is given
            elseif ischar(units)
                us(1:dims) = Unit(units); % convert to Unit and repeat
            elseif iscell(units)
                for u=1:dims
                    us(u) = Unit(units{u});
                end
            else
                error('Invalid input for ''units''');
            end     
            obj.units = us;
        end
        function array = with_data(obj, data)
            %% returns an object with the same properties (pitch, units) as template object 'obj'
            % but with the data set to 'data'
            %if size(data) ~= size(obj.data)
            %    error('Data must be the same size as the data array in the template object');
            %end
            array = obj;
            array.data = data;
        end            
        function coords = coordinates(obj, dim)
            %% return coordinate range for given dimension
            %If the array uses standard coordinates (is_fft = false),
            %this function returns 'almost symmetrical' coordinates
            %with the 0 positioned at ceil(end/2)
            %
            %If the array uses fft coordinates (is_fft = true),
            %this function returns 'fft' coordinates that reflect
            %the spacing between the Fourier components (including 2pi
            %pre-factor, so this function returns angular frequencies)
            %Also, in this case the coordinates are ifftshift-ed to match 
            %the coorinates of the fft data. So, the 0 will be positioned
            %in the first element of the array
            %
            N = size(obj, dim);
            if ~obj.is_fft
                coords = ((0:N-1)-floor(N/2)) * obj.pitches(dim);
            else
                coords = ((0:N-1)-floor(N/2)) * 2.0 * pi / (obj.pitches(dim) * size(obj.data, dim));
                coords = ifftshift(coords);
            end
        end
        function u = unit(obj, dim)
            if dim > ndims(obj.data)
                u = Unit(); %unitless
            else
                if ~obj.is_fft
                    u = obj.units(dim);
                else
                    u = inverse(obj.units(dim));
                end
            end
        end
        function N = size(obj, varargin)
            % returns size of the data in the SizedArray
            N=size(obj.data,varargin{:});
        end
        function N = end(obj, k, n)
            N = builtin('end', obj.data, k, n);
        end
        function F = fft(obj)
            % Perform a Fourier transform of a 1-D SizedArray
            % This function throws an error if the array is not 1-dimensional
            % because usually you don't want to take the 1-D Fourier
            % transform of a higher dimensional data set.
            %
            if ~isvector(obj)
                error('SizedArray only supports 1-D fft on 1-D arrays. (see fftn)');
            end
            F = fftn(obj);
        end
        function F = ifft(obj)
            % Perform a Fourier transform of a 1-D SizedArray
            % This function throws an error if the array is not 1-dimensional
            % because usually you don't want to take the 1-D Fourier
            % transform of a higher dimensional data set.
            %
            if ~isvector(obj)
                error('SizedArray only supports 1-D ifft on 1-D arrays. (see ifftn)');
            end
            F = ifftn(obj);
        end
        function F = fft2(obj)
            % Perform a Fourier transform of a 2-D SizedArray
            % This function throws an error if the array is not 2-dimensional
            % because usually you don't want to take the 2-D Fourier
            % transform of a higher dimensional data set.
            %
            if ~ismatrix(obj.data)
                error('SizedArray only supports 2-D fft on 2-D arrays. (see fftn)');
            end
            F = fftn(obj);
        end
        function F = ifft2(obj)
            % Perform a Fourier transform of a 2-D SizedArray
            % This function throws an error if the array is not 2-dimensional
            % because usually you don't want to take the 2-D Fourier
            % transform of a higher dimensional data set.
            %
            if ~ismatrix(obj.data)
                error('SizedArray only supports 2-D ifft on 2-D arrays. (see ifftn)');
            end
            F = ifftn(obj);
        end
        function F = fftn(obj)
           % Performs a Fourier transform over all dimensions of the array
           % returns a new SizedArray object with correct Fourier
           % coordinates in 'angular frequency' (i. e. including 2?
           % prefactor)
           if obj.is_fft
               error('fftn can only transform from standard to fft coordinates. Current coordinate type = fft');
           end
           F = obj.with_data(fftn(obj.data));
           F.is_fft = true;
        end
        function F = ifftn(obj)
           % Performs a Fourier transform over all dimensions of the array
           % returns a new SizedArray object with correct Fourier
           % coordinates in 'angular frequency' (i. e. including 2?
           % prefactor)
           if ~obj.is_fft
               error('ifftn can only transform from fft to standard coordinates. Current coordinate type = standard');
           end
           F = obj.with_data(ifftn(obj.data));
           F.is_fft = false;
        end
        function C = times(A, B)
            if ~isa(B, 'SizedArray') %only A is SizedArray
                C = A.with_data(A.data .* B);
            elseif ~isa(A, 'SizedArray') %only B is SizedArray
                C = B.with_data(B.data .* A);
            else
                assert_same_size(A, B); %both are SizedArray
                C = A.with_data(B.data .* A.data);
            end
        end
        function C = mtimes(A, B)
            if ~isa(B, 'SizedArray') %only A is sizedarray
                validateattributes(B, {'single', 'double'}, {'scalar'});
                C = A.with_data(A.data * B);
                return;
            end
            if isa(A, 'SizedArray')
                % both are sizedarray
                if ~same_size(A, B) %both are sizedarray: check dimensions
                    error('Dimensions (including their units) must agree');
                else
                    C = A.with_data(A.data * B.data);
                end
            else
                % only B is sizedarray
                C = B.with_data(B.data * A);
            end
        end
        function same = same_size(A, B)
            same = isequal(A.pitches, B.pitches) && isequal(size(A), size(B)) && isequal(A.units, B.units) && isequal(A.is_fft, B.is_fft);
        end
        function B = subsref(obj, S)
            % Implement indexing operator for sized array
            % This function is called when using the indexing / slicing
            % operator: e. g. A(:,2)
            % (unfortunately, it is also called when accessing properties
            % or methods, which is why this implementation needs to be so
            % complucated).
            if strcmp(S(1).type, '()')
                % slice data
                B = obj.with_data(subsref(obj.data, S)); 

                % matlab removes trailing singleton dimensions
                % make sure that we also remove the corresponding units &
                % pitches
                newdim = ndims(B.data);
                B.pitches = B.pitches(1:newdim); 
                B.units = B.units(1:newdim);
            else
                % used for accessing properties. Unfortunately, calling
                % the 'builtin' version of subsref ignores the fact
                % that our properties are protected, so we need to
                % impose that manually
                if any(strcmp(properties(obj), S(1).subs)) || any(strcmp(methods(obj), S(1).subs))
                     B = builtin('subsref', obj, S); 
                else
                    error(['No public property ''' S(1).subs ''' for class ''SizedArray''']);
                end
            end
        end
        function A = subsasgn(A, S, B)
            % implements index assignment operator for SizedArray:
            % a(1:2,:)=b  (where b can be a sizedarray or an ordinary matrix)
            
            % we don't want to support arrays of SizedArray objects (because it gets messy and confusing)
            % therefore, we throw an error when one tries to create one.
            % the error is identical to the one you get when trying to
            % assign a matrix to a matrix element
            if prod(builtin('size', A)) > 1
                error('In an assignment  A(:) = B, the number of elements in A and B must be the same.');
            end
            if isa(B, 'SizedArray')
                d = subsasgn(A.data, S, B.data);
            else
                d = subsasgn(A.data, S, B);
            end
            A = A.with_data(d);
        end
        function B = squeeze(obj)
            % removes singleton dimensions from the data array (as in
            % the version of squeeze operating on matrices)
            % also updates the pitches and units accordingly
            %
            B = obj;
            if ~ismatrix(B.data)
                % we need to squeeze both the array data and the units
                bsize = size(B.data);
                keep = bsize > 1;
                bsize = bsize(keep); % remove singleton dimensions
                
                % remove pitches and units for singleton dimensions
                B.pitches = B.pitches(keep);
                B.units = B.units(keep);
                
                % Make sure we keep the matrix at least 2-D by adding
                % 1 or 2 trailing singleton dimensions.
                % These added dimensions will be unitless. 
                blen = length(bsize);
                if blen < 2
                    bsize((blen+1):2) = 1;
                    B.pitches((blen+1):2) = 1;
                    B.units((blen+1):2) = Unit();
                end
                B.data = reshape(B.data, bsize);
            end
        end
%        function get.
            %subsref
            %subsasgn
        function imagesc(obj, varargin)
            labels = {['x ', char(obj.unit(2), true)], ['y ', char(obj.unit(1), true)], ['z ', char(obj.unit(3), true)]};
            labels(size(obj.data)==1) = []; % remove labels for singleton dimensions
            if length(labels) < 2
                labels{2} = []; 
            end
            s = squeeze(obj);
            if isreal(s.data)
                d = s.data;
            else
                d = abs(s.data);
            end
            imagesc(limits(obj,2), limits(obj,1), d, varargin{:});
            xlabel(labels(1));
            ylabel(labels(2));
            axis image;
        end
    end
    methods (Access = private)
        function l = limits(obj, dim)
            %returns left & right bound only (for use in image/imagesc). 
            c = coordinates(obj, dim);
            l = [min(c), max(c)]; %todo: speed up
        end
        function assert_same_size(A, B)
            if ~same_size(A, B)
                error('Arrays don''t have the same size');
            end
        end
    end
end

