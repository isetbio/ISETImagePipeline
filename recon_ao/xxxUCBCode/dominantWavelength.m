function [lambda,llmStim] = dominantWavelength(proportionRed)
% use to calculate equivalent monochromatic wavelength that is metameric to
% mixture of 543 nm (green) and 680 nm (red) primaries
% proportionRed = proportion of total luminance (trolands) coming from red (680 nm) laser

% Written by AE Boehm for AOSLO color experiments. 

% load the cone fundamentals
load('StockSharpe2000'); 

% wvl = (500:780)';
% % scale so that eew is at (0.7,1)
wvl = StockSharpe2000(:,1);
lms = 10.^StockSharpe2000(:,2:4);
% l = lms(:,1) .* 0.7 ./ sum(lms(:,1));
% m = lms(:,2) .* 0.3 ./ sum(lms(:,2));
% s = lms(:,3) .* 1.0 ./ sum(lms(:,3));

% scale the cone sensitivities so that Vlambda = l+m
%  StockSharpe2000(:,2:end) = 10.^StockSharpe2000(:,2:end); 
% wvl = StockSharpe2000(:,1); 
l = lms(:,1).*Lscale;
m = lms(:,2).*Mscale;
s = lms(:,3).*Sscale;

% macleod-boynton coordinates
llm = l./(l+m);
slm = s./(l+m);

llm680 = llm(wvl==680);
slm680 = slm(wvl==680);
llm543 = llm(wvl==543);
slm543 = slm(wvl==543);

% In a space where Cartesian coordinates represent the excitations of the three cone types 
% involved in color vision, a plane of constant luminance provides a chromaticity diagram 
% in which excitation of each cone type (at constant luminance) is represented by a linear 
% scale (horizontal or vertical), and in which the center-of-gravity rule applies with weights 
% proportional to luminance (see p. 1183 of MacLeod & Boynton, 1979).
llmStim = llm680.*proportionRed + llm543.*[1-proportionRed];
slmStim = slm680.*proportionRed + slm543.*[1-proportionRed];
idx = wvl>=wvl(llm == min(llm)) & wvl<=wvl(llm == max(llm)); % removes errors from interp because of inflection in wvl vs. llm at shorter wavelengths.

lambda = interp1(llm(idx),wvl(idx),llmStim, 'linear', 'extrap');
end