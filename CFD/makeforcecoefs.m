function S = makeforcecoefs(S, U, varargin)
% function S = makeforcecoefs(S, U, varargin)
% Converts forces and impulses in the S structure to force and impulse coefficients.
%
% Units:
%  forcescale -> [g/cm^3] * [cm] * [cm/s]^2 = g cm / s^2 = dynes / cm
%  S.Faxialtotcoef -> [dynes/cm] / [dynes/cm] = dimensionless

opt.bodylen = 4*pi;
opt.rho = 1;
opt.perimeter = 25.3203;
opt.height = 0.2906;

opt = parsevarargin(opt,varargin, 3);

%a 2D "area" = length
wettedsurface = opt.perimeter;

forcescale = 0.5 * opt.rho * wettedsurface * U^2;
impulsescale = 0.5 * opt.rho * wettedsurface * opt.bodylen * U;

% S.Faxialtcoef = S.Faxialt / forcescale;
% S.Faxialncoef = S.Faxialn / forcescale;
% S.Flateraltcoef = S.Flateralt / forcescale;
% S.Flateralncoef = S.Flateraln / forcescale;
S.Faxialtotcoef = S.Faxialtot / forcescale;
S.Flateraltotcoef = S.Flateraltot / forcescale;
 
% S.Ithrusttcoef = S.Ithrustt / impulsescale;
% S.Ithrustncoef = S.Ithrustn / impulsescale;
% S.Idragtcoef = S.Idragt / impulsescale;
% S.Idragncoef = S.Idragn / impulsescale;

S.Ithrusttotcoef = sum(sum(S.Ithrustt+S.Ithrustn,1),3) / impulsescale;
S.Idragtotcoef = sum(sum(S.Idragt+S.Idragn,1),3) / impulsescale;
S.Inormcoef = sum(sum(S.Ithrustn+S.Idragn,1),3) / impulsescale;
S.Itancoef = sum(sum(S.Ithrustt+S.Idragt,1),3) / impulsescale;

S.Iaxialtotcoef = S.Iaxialtot / impulsescale;
S.Ilateraltotcoef = S.Ilateraltot / impulsescale;
