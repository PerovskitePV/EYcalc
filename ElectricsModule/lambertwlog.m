function w = lambertwlog(logx)
% LAMBERTWLOG - Lambert W function with logarithmic input
% W = LAMBERTWLOG(LOGX) solves the equation LOGX = W+log(W) for W>=0,
% equivalent to solving the Lambert W function exp(LOGX) = W*exp(W).
% Argument LOGX must be real.  Note that this restricts the domain to
% nonnegative exp(LOGX) and thus the principal branch of the function.
% 
% See also: LAMBERTW
%
% Credits to: Michael (2021). Lambert W function (logarithmic input) 
% (https://www.mathworks.com/matlabcentral/fileexchange/57239-lambert-w-function-logarithmic-input), MATLAB Central File Exchange.
% Retrieved February 8, 2021.

if ~isreal(logx) && any(imag(logx))
	error('Complex arguments not supported: use lambertw() instead.')
end
toosmall = logx<log(realmin); % the solution here will be less than realmin
small = ~toosmall & logx<0.585; % best split for fewest iterations - usually 4 for small and 4 for large
large = ~(toosmall | small);
prec = eps; % relative precision
itermax = 5; % usually no progress past 4 iterations

xsmall = exp(logx(small));
wsmall = (sqrt(1+4*xsmall)-1)/2; % initialize approximating that exp(W)==1+W
wprev = inf;
iters = 0;
while any(abs(wsmall - wprev) > wsmall*prec) && iters < itermax
	iters = iters+1;
	wprev = wsmall;
	wsmall = wsmall - (wsmall - xsmall.*exp(-wsmall))./(wsmall+1); % Newton iteration
end

logxlarge = logx(large);
wlarge = logxlarge; % initialize approximating that log(W)==0
wprev = inf;
iterl = 0;
while any(abs(wlarge - wprev) > wlarge*prec) && iterl < itermax
	iterl = iterl+1;
	wprev = wlarge;
	wlarge = wlarge - (wlarge+log(wlarge)-logxlarge) .* wlarge./(wlarge+1); % Newton iteration
end

w = zeros(size(logx)); % initialize, set toosmall's to 0
w(small) = wsmall;
w(large) = wlarge;

% fprintf('%i iterations\n',iters+iterl)
end