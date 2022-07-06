function s = movingstd2(A)

% size of the array
n = size(A);

% Improve the numerical analysis by subtracting off the array mean
% this has no effect on the standard deviation, but when the mean
% islarge, the formula used will incur numerical issues.
A = A - mean(A(:));

% scale the array to have unit variance too. will put that
% scale factor back into the result at the end
Astd = std(A(:));
A = A./Astd;

% we will need the squared elements 
A2 = A.^2;

% we also need an array of ones of the same size as A. This will let us
% count the number of elements in each truncated window near the edges.
wuns = ones(size(A));

% convolution kernel
kernel = ones(2,2);

% compute the std using:
%     std = sqrt((sum(x.^2) - (sum(x)).^2/n)/(n-1))
N = conv2(wuns,kernel,'same');
s = sqrt((conv2(A2,kernel,'same') - ((conv2(A,kernel,'same')).^2)./N)./(N-1));

% catch any complex cases that may have fallen through the cracks.
s(imag(s) ~= 0) = 0;
% restore the scale factor that was used before to normalize the data
s = s.*Astd;
