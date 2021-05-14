function B = imgaussFWHMfilt3(varargin)
%IMGAUSSFILT3 3-D Gaussian filtering of 3-D images
%
%   B = imgaussfilt3(A) filters 3-D image A with a 3-D Gaussian smoothing
%   kernel with standard deviation of 0.5.
%
%   B = imgaussfilt3(A,SIGMA) filters 3-D image A with a 3-D Gaussian
%   smoothing kernel with standard deviation specified by SIGMA. SIGMA can
%   be a scalar or a 3-element vector with positive values. If sigma is a
%   scalar, a cube Gaussian kernel is used.
%
%   B = imgaussfilt3(___,Name,Value,...) filters 3-D image A with a 3-D
%   Gaussian smoothing kernel with Name-Value pairs used to control aspects
%   of the filtering.
%
%   Parameters include:
%
%   'FilterSize'    -   Scalar or 3-element vector, of positive, odd
%                       integers that specifies the size of the Gaussian
%                       filter. If a scalar Q is specified, then a square
%                       Gaussian filter of size [Q Q Q] is used. 
%                       Default value is 2*ceil(2*SIGMA)+1.
%
%   'Padding'       -   String or numeric scalar that specifies padding to
%                       be used on image before filtering. 
%
%                       If a scalar (X) is specified, input image values
%                       outside the bounds of the image are implicitly
%                       assumed to have the value X. 
%
%                       If a string is specified, it can be 'replicate',
%                       'circular' or 'symmetric'. These options are
%                       analogous to the padding options provided by
%                       imfilter. 
%
%                       'replicate'
%                       Input image values outside the bounds of the image
%                       are assumed equal to the nearest image border
%                       value.
%
%                       'circular'
%                       Input image values outside the bounds of the image
%                       are computed by implicitly assuming the input image
%                       is periodic.
%
%                       'symmetric'
%                       Input image values outside the bounds of the image
%                       are computed by mirror-reflecting the array across
%                       the array border. 
%
%                       Default value is 'replicate'.
%
%   'FilterDomain'  -   String that specifies domain in which filtering is
%                       performed. It can be 'spatial', 'frequency' or
%                       'auto'. For 'spatial', convolution is performed in
%                       the spatial domain, for 'frequency', convolution is
%                       performed in the frequency domain and for 'auto',
%                       convolution may be performed in spatial or
%                       frequency domain based on internal heuristics.
%                       Default value is 'auto'.
%
%
%   Class Support
%   -------------
%   The input image A must be a real, non-sparse matrix of 3 dimension of
%   the following classes: uint8, int8, uint16, int16, uint32, int32,
%   single or double.
%
%   Notes
%   -----
%   1. If the image A contains Infs or NaNs, the behavior of imgaussfilt3
%      for frequency domain filtering is undefined. This can happen when
%      the 'FilterDomain' parameter is set to 'frequency' or when it is set
%      to 'auto' and frequency domain filtering is used internally. To
%      restrict the propagation of Infs and NaNs in the output in a manner
%      similar to imfilter, consider setting the 'FilterDomain' parameter
%      to 'spatial'.
%
%   2. When 'FilterDomain' parameter is set to 'auto', an internal
%      heuristic is used to determine whether spatial or frequency domain
%      filtering is faster. This heuristic is machine dependent and may
%      vary for different configurations. For optimal performance, consider
%      comparing 'spatial' and 'frequency' to determine the best filtering
%      domain for your image and kernel size.
%
%   3. When the 'Padding' parameter is not specified, 'replicate' padding 
%      is used as the default, which is different from the default used in
%      imfilter.
%
%
%   Example 1
%   ---------
%   This example smooths an MRI volume with a 3-D Gaussian filter.
%
%   vol = load('mri');
%   figure, montage(vol.D), title('Original image volume')
%
%   siz = vol.siz;
%   vol = squeeze(vol.D);   
%   sigma = 2;
%
%   volSmooth = imgaussfilt3(vol, sigma);
% 
%   figure, montage(reshape(volSmooth,siz(1),siz(2),1,siz(3)))
%   title('Gaussian filtered image volume')
%
%   See also imgaussfilt, imfilter.

% Copyright 2014, The MathWorks, Inc.

narginchk(1, Inf);

[A, options] = parseInputs(varargin{:});

sigma       = options.Sigma;
hsize       = options.FilterSize;
padding     = options.Padding;
domain      = options.FilterDomain;

[domain, ippFlag] = chooseFilterImplementation(A, hsize, domain);

B = spatialGaussianFilter(A, sigma, hsize, padding, ippFlag);

end

%--------------------------------------------------------------------------
% Spatial Domain Filtering
%--------------------------------------------------------------------------
function A = spatialGaussianFilter(A, sigma, hsize, padding, ippFlag)

sizeA = size(A);
dtype = class(A);
outSize = [sizeA ones(1,3-ndims(A))];

[A, padSize] = padImage(A, hsize, padding);

[hCol,hRow,hSlc] = createSeparableFWHMkernel(4);

% cast to double to preserve precision over separable filter calls.
if ~isa(A, 'double')
    A = double(A);
end

% Row filtering
outSize1D = [size(A) ones(1,3-ndims(A))];
outSize1D(2) = outSize(2);

startRow = zeros(size(padSize));
startRow(2) = padSize(2);

A = imfiltermexCaller(A, outSize1D, hRow, startRow, ippFlag);

% Column filtering
outSize1D(1) = outSize(1);

startCol = zeros(size(padSize));
startCol(1) = padSize(1);

A = imfiltermexCaller(A, outSize1D, hCol, startCol, ippFlag);

% Slice filtering
outSize1D(3) = outSize(3);

startSlc = zeros(size(padSize));
startSlc(3) = padSize(3);

% IPP only supports 2-D filtering, disable for filtering along 3rd
% dimension.
ippFlag = false;

A = imfiltermexCaller(A, outSize1D, hSlc, startSlc, ippFlag);

if ~isa(A,dtype)
    A = cast(A,dtype);
end
    
end

function A = imfiltermexCaller(A, outSize, h, start, ippFlag)

conn = [];
nz   = [];

convMode = ippFlag; % Set the convolution mode to true when IPP is enabled
                    % to prevent unnecessary (since Gaussian kernels are 
                    % symmetric) rotation of kernel. Conversly, when IPP is
                    % disabled, set convolution mode to false so home-brew 
                    % does not rotate the kernel.
sameSize = true;

if ~ippFlag
    conn = h~=0;
    nz   = h(conn);
end

A = imfilter_mex(A, outSize, h, nz, conn, start, sameSize, convMode, ippFlag);

end

function TF = useIPPL(A, outSize)

% Querying the Java preference for IPP use is computationally expensive,
% particularly in cases where IMGAUSSIANFILTER is called in a loop. As a
% performance optimization, store the preference as persistent state. This
% state is reset whenever iptsetpref('UseIPPL',___) is called.
persistent prefFlag;
needToQueryJavaPreference = isempty(prefFlag);
if needToQueryJavaPreference
    prefFlag = iptgetpref('UseIPPL');
end

TF = prefFlag;

end

function [hcol,hrow,hslc] = createSeparableFWHMkernel(fwhm)

def_x=2.2; %mm
def_y=2.2; %mm
def_z=2.8; %mm

sigma_g=fwhm/(2*sqrt(2*log(2))); % sigma gaussienne, en mm

%pour la taille du noyau on prend de façon à avoir plus de 6 sigma, en
%arrondissant au supérieur et pour un nombre de voxel impair
%de plus on veut un noyau homogene en dimension (de forme cubique), donc on
%se cale sur x ou y
hsize=ceil(6*sigma_g/def_x);


noyau_x=zeros(hsize,1);
noyau_y=zeros(hsize,1);
noyau_z=zeros(hsize,1);

%bornes inf et sup des intervalles pour le calcul des integrales
bornes_inf=zeros(hsize,1);
bornes_sup=zeros(hsize,1);

bornes_inf((hsize+3)/2:end)=0.5:1:((hsize-3)/2)+0.5;
bornes_sup((hsize+1)/2:end)=0.5:1:((hsize-1)/2)+0.5;

noyau_x=(erf(bornes_sup.*(def_x/sigma_g))-erf(bornes_inf.*(def_x/sigma_g)))./2;
noyau_y=(erf(bornes_sup.*(def_y/sigma_g))-erf(bornes_inf.*(def_y/sigma_g)))./2;
noyau_z=(erf(bornes_sup.*(def_z/sigma_g))-erf(bornes_inf.*(def_z/sigma_g)))./2;

noyau_x(1:(hsize-1)/2)=noyau_x(end:-1:(hsize+3)/2);
noyau_x((hsize+1)/2)=2*noyau_x((hsize+1)/2);
noyau_y(1:(hsize-1)/2)=noyau_y(end:-1:(hsize+3)/2);
noyau_y((hsize+1)/2)=2*noyau_y((hsize+1)/2);
noyau_z(1:(hsize-1)/2)=noyau_z(end:-1:(hsize+3)/2);
noyau_z((hsize+1)/2)=2*noyau_z((hsize+1)/2);

hcol = noyau_y;

hrow = noyau_x;
hslc = noyau_z;


hrow = reshape(hrow, 1, hsize);
hslc = reshape(hslc, 1, 1, hsize);

end

%--------------------------------------------------------------------------
% Common Functions
%--------------------------------------------------------------------------
function [domain, ippFlag] = chooseFilterImplementation(A, hsize, domain)

ippFlag = useIPPL(A, size(A));

if strcmp(domain, 'auto')
    domain = chooseFilterDomain3FWHM(A, hsize, ippFlag);
end

end

function [A, padSize] = padImage(A, hsize, padding)

padSize = computePadSize(size(A), hsize);

if ischar(padding)
    method = padding;
    padVal = [];
else
    method = 'constant';
    padVal = padding;
end
A = padarray_algoFWHM(A, padSize, method, padVal, 'both');

end

function padSize = computePadSize(sizeA, sizeH)

rankA = numel(sizeA);
rankH = numel(sizeH);

sizeH = [sizeH ones(1,rankA-rankH)];

padSize = floor(sizeH/2);

end

%--------------------------------------------------------------------------
% Input Parsing
%--------------------------------------------------------------------------
function [A, options] = parseInputs(varargin)

A = varargin{1};

supportedClasses = {'uint8','uint16','uint32','int8','int16','int32','single','double'};
supportedImageAttributes = {'real','nonsparse', '3d'};
validateattributes(A, supportedClasses, supportedImageAttributes, mfilename, 'A');

% Default options
options = struct(...
        'Sigma',        [.5 .5 .5],...
        'FilterSize',   [3 3 3],...
        'Padding',      'replicate',...
        'FilterDomain', 'auto');

beginningOfNameVal = find(cellfun(@isstr,varargin),1);

if isempty(beginningOfNameVal) && length(varargin)==1
    %imgaussfilt3(A)
    return;
elseif beginningOfNameVal==2
    %imgaussfilt3(A,'Name',Value)
elseif (isempty(beginningOfNameVal) && length(varargin)==2) || (~isempty(beginningOfNameVal) && beginningOfNameVal==3)
    %imgaussfilt3(A,sigma,'Name',Value,...)
    %imgaussfilt3(A,sigma)
    options.Sigma = validateSigma(varargin{2});
    options.FilterSize = computeFilterSizeFromSigma(options.Sigma);
else
    error(message('images:imgaussfilt:tooManyOptionalArgs'));
end

numPVArgs = length(varargin) - beginningOfNameVal + 1;
if mod(numPVArgs,2)~=0
    error(message('images:imgaussfilt:invalidNameValue'));
end

ParamNames = {'FilterSize', 'Padding', 'FilterDomain'};
ValidateFcn = {@validateFilterSize, @validatePadding, @validateFilterDomain};

for p = beginningOfNameVal : 2 : length(varargin)-1
    
    Name = varargin{p};
    Value = varargin{p+1};
    
    idx = strncmpi(Name, ParamNames, numel(Name));
    
    if ~any(idx)
        error(message('images:imgaussfilt:unknownParamName', Name));
    elseif numel(find(idx))>1
        error(message('images:imgaussfilt:ambiguousParamName', Name));
    end
    
    validate = ValidateFcn{idx};
    options.(ParamNames{idx}) = validate(Value);
    
end

end

function sigma = validateSigma(sigma)

validateattributes(sigma, {'numeric'}, {'real','positive','finite','nonempty'}, mfilename, 'Sigma');

if isscalar(sigma)
    sigma = [sigma sigma sigma];
end

if numel(sigma)~=3
    error(message('images:imgaussfilt:invalidLength3', 'Sigma'));
end

sigma = double(sigma);

end

function filterSize = computeFilterSizeFromSigma(sigma)

filterSize = 2*ceil(2*sigma) + 1;

end

function filterSize = validateFilterSize(filterSize)

validateattributes(filterSize, {'numeric'}, {'real','positive','integer','nonempty','odd'}, mfilename, 'FilterSize');

if isscalar(filterSize)
    filterSize = [filterSize filterSize filterSize];
end

if numel(filterSize)~=3
    error(message('images:imgaussfilt:invalidLength3', 'FilterSize'));
end

filterSize = double(filterSize);

end

function padding = validatePadding(padding)

if ~ischar(padding)
    validateattributes(padding, {'numeric','logical'}, {'real','scalar','nonsparse'}, mfilename, 'Padding');
else
    padding = validatestring(padding, {'replicate','circular','symmetric'}, mfilename, 'Padding');
end

end

function filterDomain = validateFilterDomain(filterDomain)

filterDomain = validatestring(filterDomain, {'spatial','frequency','auto'}, mfilename, 'FilterDomain');

end
