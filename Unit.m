classdef Unit < matlab.mixin.CustomDisplay
    %Unit unit of measure (e. g. ?m, s, J, etc.)
    properties (Access=protected)
        unit_count
    end
    methods
        function obj = Unit(str)
            %% Convert a string to a unit
            % The string must have a form like 'm' or 'm s^-1' or 'Hz^0.5'
            % used at the moment
            % todo:alternatively, the form 'm/s' can be used. 
            if nargin == 0
                obj = Unit('');
                return
            end
            elements = strsplit(str, ' ');
            % check for empty string
            if isempty(elements{1})
                N = 0;
            else
                N = length(elements);
            end
            obj.unit_count = zeros(N * 2, 1, 'uint32');
            for e=1:N
                ec = strsplit(elements{e}, '^');
                % convert unit characters to number
                chars = zeros(4, 1, 'uint8');
                chars(1:length(ec{1})) = ec{1};
                count = 1;
                if length(ec) > 1
                    count = str2double(ec{2});
                end
                obj.unit_count(e) = typecast(chars, 'uint32');
                obj.unit_count(e+N) = typecast(single(count), 'uint32');
            end
            % sort units based on type to arrive at a canonical form.
            % (This way, m s^-1 and s^-1 m are converted to identical units).
            [sorted_units, ix] = sort(obj.unit_count(1:N));
            obj.unit_count(1:N) = sorted_units;
            obj.unit_count(N+1:end) = obj.unit_count(ix+N);
        end
        function equal = eq(A, B)
            %% test for equality
            % for a scalar A and B, we could simply compare unit_count
            % however, the code here is more complicated since we can also
            % compare arrays with each other, or scalars with arrays
            if ~isequal(size(A), size(B)) && ~isscalar(A) && ~isscalar(B)
                error('Matrix dimensions must agree.');
            end;
            equal = zeros(size(A), 'logical');
            if isa(A, 'Unit') && isa(B, 'Unit') 
                for u=1:numel(A)
                    equal(u) = isequal(A(u).unit_count,B(u).unit_count); %this is the actual comparison
                end
            end
        end
        function inv = inverse(obj)
            %% Return inverse unit: 1/unit
            inv = obj;
            inv.unit_count((end/2+1):end) = typecast(-typecast(obj.unit_count((end/2+1):end), 'single'), 'uint32');
        end
        function str = char(obj, latex)
            %% Convert unit to string representation
            %  todo sort units based on power?
            out = '';
            N = length(obj.unit_count)/2;
            if N > 0
                for e=1:N
                    unit = deblank(char(typecast(obj.unit_count(e), 'uint8')));
                    if latex && unit(1) == 'u' && length(unit) > 1 %special case for mu (micrometer, etc.)
                        unit = ['\mu', unit(2:end)];
                    end

                    count = typecast(obj.unit_count(e+N), 'single');
                    if count == 1
                        out = [out, ' ', unit];
                    else
                        if latex
                        	out = [out, ' ', unit, '^{', num2str(count), '}'];
                        else
                            out = [out, ' ', unit, '^', num2str(count)];
                        end
                    end
                end
                str = ['[', out(2:end), ']']; %strip leading space
            else
                str = '[-]'; %unitless
            end
        end
        %% todo: add functions to multiply and divide units
        %% todo: add conversion to SI units?
        %% todo: add orders of magnitude? ns, mm, etc.
    end
    methods (Access=protected)
        function displayScalarObject(obj)
            %% Function to display the unit as a human readable string
            disp(char(obj));
        end
        function displayNonScalarObject(objAry)
            %(modified from Matlab example code)
            dimStr = matlab.mixin.CustomDisplay.convertDimensionsToString(objAry);
            cName = matlab.mixin.CustomDisplay.getClassNameForHeader(objAry);
            fprintf('%s %s:\n',dimStr, cName);
            for o = 1:size(objAry(:))
                disp(objAry(o));
            end;
        end
    end
end

