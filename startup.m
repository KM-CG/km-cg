clear all;
restoredefaultpath
me    = mfilename;                                          % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located
addpath(mydir(1:end-1))
addpath(genpath([mydir 'impl']));
addpath(genpath([mydir 'data']));
addpath(genpath([mydir 'visualization']));
addpath(genpath([mydir 'hyper_parameters']));
clear me;
clear mydir;

mpg = [0,0.4717,0.4604]; % color [0,125,122]
dre = [0.4906,0,0]; % color [130,0,0]
ora = [255,153,51] ./ 255;
blu = [0,0,0.509];
gra = 0.5 * ones(3,1);

lightmpg = [1,1,1] - 0.5 * ([1,1,1] - mpg);
lightdre = [1,1,1] - 0.5 * ([1,1,1] - dre);
lightblu = [1,1,1] - 0.5 * ([1,1,1] - blu);
lightora = [1,1,1] - 0.5 * ([1,1,1] - ora);

cya = lightmpg;

mpg2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-mpg));
dre2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-dre));
blu2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-blu));
ora2white = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,0.6,2024)').^0.5,[1,1,1]-ora));

cya2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightmpg);
red2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightdre);
blu2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightblu);
ora2black = bsxfun(@times,(linspace(0,0.6,2024)').^0.5,lightora);

mpg2whiteLin = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,1,2024)'), [1,1,1]-mpg));
dre2whiteLin = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,1,2024)'),[1,1,1]-dre));
ora2whiteLin = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,1,2024)'),[1,1,1]-ora));
blu2whiteLin = bsxfun(@minus,[1,1,1],bsxfun(@times,(linspace(0,1,2024)'),[1,1,1]-blu));
ora2dreLin = bsxfun(@minus,dre,bsxfun(@times,(linspace(0,1,2024)'),dre-ora));
blu2mpgLin = bsxfun(@minus,blu,bsxfun(@times,(linspace(0,1,2024)'),blu-mpg));
cmapMPGDre= [flipud(mpg2white); dre2white];
cmapMPGDre= cmapMPGDre(1 : 2 : end, :); % take only every 2nd element for matlab2tikz
cmapMPGDreLin= [flipud(mpg2whiteLin); dre2whiteLin];
cmapMPGDreLin= cmapMPGDreLin(1 : 2 : end, :);
cmapWhiteOraDre = [ora2dreLin; flipud(ora2whiteLin)];
cmapWhiteOraDre= cmapWhiteOraDre(1 : 2 : end, :); 
cmapWhiteMPGblu = [blu2mpgLin; flipud(mpg2whiteLin)];
cmapWhiteMPGblu= cmapWhiteMPGblu(1 : 2 : end, :);