function [] = infra_dist_azimuth()
%
% This function calculates the distances along a great circle and theoretical 
% delay times of infrasonic signals at infrasound stations of the CTBTO IMS network
% It also accounts for travel times associated with various propagation 
% channels (tropospheric, stratospheric and thermospheric)
% The user has a choice to output a map, and print time segments to be examined.
% 
% input variables: (date, time, source)
%       date [dd mm yyyy]
%       time [hh mm ss]
%       source [lat lon]
% Last updated: July 2021
% v2.7
%
% Elizabeth A. SILBER, May 2010 (c)
% e-mail: esilber@uwo.ca
% The University of Western Ontario
%
% ========================================================================
evdate = input('Enter the event date [dd mm yyyy]:  ');
time = input ('Enter the event time [hh mm ss]: ');
source = input ('Enter the coordinates of the source [lat lon]: ');

% evdate = [29 04 2010]; 
% source = [-6 107]; 
% time = [09 30 00]; 

eventTime  = date2unixsecs(evdate(3),evdate(2),evdate(1),time(1),time(2),time(3));

k = input('Do you wish to include delay time in hh:mm:ss and ETA? 0 = no, 1 = yes    ');
map = input('Do you wish to print a map? 0 = no, 1 = yes   '); 

[stn,name,mapName] = ISMstations();
n = length(stn);

% define variables
dist = zeros(n,1); rangdeg = zeros(n,1);
az = zeros(n,1);
delay = zeros(n,3); 
eta0 = zeros(n,3);
etatr = zeros(n,6); etastr = zeros(n,6); etath = zeros(n,6);
time1 = zeros(n,3); time2 = zeros(n,3); time3 = zeros(n,3);


% Produce a .txt file with results
file = fopen(['DazD-',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'.txt'],'w');
fprintf(file,('Date of the event [dd-mm-yyyy]: %02.0f-%02.0f-%4.0f\n'),evdate(1),evdate(2),evdate(3));     
fprintf(file,('Time of the event [hh:mm:ss]:   %02.0f:%02.0f:%02.0f   UTC\n'),time(1),time(2),time(3));
fprintf(file,('Location of the event [lat(N), lon(E)]:   %02.2f, %02.2f% \n'),source(1),source(2));
fprintf(file,'\n================================================================================\n');
fprintf(file,'STN \tDist  \tBack Az \t Delay(s) \t Strato Del \t Thermo Del\n');
fprintf(file,'    \t(km)  \t(deg)   \t @340 m/s \t@ 285 m/s   \t @220 m/s\n');
fprintf(file,'================================================================================\n');

for i = 1:n,
    rangdeg(i) = distance (source(1),source(2),stn(i,2),stn(i,3),'deg');
    dist(i) = deg2km(rangdeg(i)); %great circle range in km
    az(i) = azimuth(stn(i,2),stn(i,3),source(1),source(2),'deg'); %back azimuth
    
    % Tropospheric (direct) arrival
    delay(i,1) = dist(i)/0.340; %delay time in seconds
    temp = sec2hms(delay(i,1)); time1(i,1:3) = hms2mat(temp); 
    eta0(i,1) = eventTime + delay(i,1); %arrival in unix time
    [etatr(i,1),etatr(i,2),etatr(i,3),etatr(i,4),etatr(i,5),etatr(i,6)] = unixsecs2date(eta0(i,1)); % arival in actual time (date, time)
    
    % Stratospheric arrival
    delay(i,2) = dist(i)/0.285; 
    temp = sec2hms(delay(i,2)); time2(i,1:3) = hms2mat(temp);
    eta0(i,2) = eventTime + delay(i,2);
    [etastr(i,1),etastr(i,2),etastr(i,3),etastr(i,4),etastr(i,5),etastr(i,6)] = unixsecs2date(eta0(i,2)); 
    
    % Thermospheric arrival
    delay(i,3) = dist(i)/0.220; 
    temp = sec2hms(delay(i,3)); time3(i,1:3) = hms2mat(temp);
    eta0(i,3) = eventTime + delay(i,3);
    [etath(i,1),etath(i,2),etath(i,3),etath(i,4),etath(i,5),etath(i,6)] = unixsecs2date(eta0(i,3)); 
        
     
    fprintf(file,'\n%s\t%5.1f\t\t%3.1f\t\t%8.1f\t%8.1f\t%8.1f\n',name(i,1:5),dist(i),az(i),delay(i,1),delay(i,2),delay(i,3));
    if k == 1,
        fprintf(file,'      \t Delay [hh:mm:ss]  \t\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\n',time1(i,1),time1(i,2),time1(i,3),time2(i,1),time2(i,2),time2(i,3),time3(i,1),time3(i,2),time3(i,3));
        fprintf(file,'      \t ETA   [hh:mm:ss]  \t\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\n',etatr(i,4),etatr(i,5),etatr(i,6),etastr(i,4),etastr(i,5),etastr(i,6),etath(i,4),etath(i,5),etath(i,6));
    end        
end

if map == 1,   
    infra_map_world1(source, evdate);
    infra_map_world2(source, evdate);
end

close all
disp('Run completed.');

%==========================================================================
function[] = infra_map_world1(source, evdate)

figure1 = figure;
worldmap world  % to plot whole world
whos -file coast.mat
load coast      % load coastline

plotm(lat, long) % plot coastline on the map
geoshow('landareas.shp', 'FaceColor', [0.35 0.6 0.35]);
geoshow('worldlakes.shp', 'FaceColor', 'cyan');    
    
[stn,name,mapNames] = ISMstations();
stn_lat = stn(:,2);     % station latitude
stn_long = stn(:,3);    % station longitude
stn_no = stn(:,1);      % station numbers

n = length(stn_no);
% make a map

hold all;

plotm(source(1), source(2), 'ko', 'MarkerFaceColor','y','MarkerSize',20);
textm(source(1), source(2), 'Bolide','FontSize',14);

for j = 1:n,
    plotm(stn_lat(j), stn_long(j), 'dk', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.2250 0.0980]);
    textm(stn_lat(j), stn_long(j), mapNames(j,1:6),'FontSize',14); 
end
set(get(groot, 'Children'), 'WindowState', 'maximized');

title(['Event of ',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3))],'FontSize',16);
print('-djpeg','-r300',['IMSstations_WorldMap-',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'-1']);  
saveas(figure1,['IMSstations_WorldMap-',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'-1.fig']);
%

%==========================================================================
function[] = infra_map_world2(source, evdate)

figure2 = figure;
worldmap world % to plot whole world
whos -file coast.mat
load coast % load coastline

plotm(lat, long) % plot coastline on the map
geoshow('landareas.shp', 'FaceColor', [0.35 0.6 0.35]);
geoshow('worldlakes.shp', 'FaceColor', 'cyan');    
    
[stn,name,mapNames] = ISMstations();
stn_lat = stn(:,2);     % station latitude
stn_long = stn(:,3);    % station longitude
stn_no = stn(:,1);      % station numbers

n = length(stn_no);
% make a map

hold all;

plotm(source(1), source(2), 'ko', 'MarkerFaceColor','y','MarkerSize',20);

for j = 1:n,
    plotm(stn_lat(j), stn_long(j), 'dk', 'MarkerSize', 8, 'MarkerFaceColor', [0.2500 0.1250 0.0980]);
end

set(get(groot, 'Children'), 'WindowState', 'maximized');

title(['Event of ',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3))],'FontSize',16);
print('-djpeg','-r300',['IMSstations_WorldMap-',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'-2']);  
saveas(figure2,['IMSstations_WorldMap-',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'-2.fig']);

% =========================================================================

function [stations,stnNames,mapNames] = ISMstations()

% the list of IMS stations as of 2020
% 
stations = [
 2,     -54.58057,  -67.30923;
 3,		-68.58390,	78.06850;
 4,     -34.59761,  116.35669;
 5,     -42.491,    147.681;
 6,		-12.1465,	96.8203;
 7,     -19.9348,   134.3295;
 8,     -16.21523,  -68.45345;
 9,     -15.638,    -48.016;
 10,    50.2065,    -96.0117;
 11,    15.25729,   -23.18388;
 13,    -27.1273,  -109.3627;
 14,    -33.6538,  -78.7960;
 17,    6.6704,     -4.8569;
 18,    77.7476,    -69.288;
 19,	11.4740,	43.17310;
 21,    -8.8678,    -140.1501;
 22,    -22.1845,   166.8459;
 23,    -49.3458,   70.2416;
 24,    -17.7493,   -149.2958;
 26,    48.8461,    13.7179;
 27,    -70.662,    -8.321;
 30,    35.3078,    140.3138;
 31,    50.4086,    58.0343;
 32,    -1.2422,    36.8272;
 33,    -19.01086,  47.30502;
 34,    47.8015,    106.40992;
 35,    -19.191,    17.577;
 36,    -43.9166,   -176.4834;
 37,	69.0758,	18.6087;
 39,    7.53547,    134.54701;
 41,    -26.342,    -57.312;
 42,    39.0423,    -28.0055;
 43,    56.72136,   37.21759;
 44,    53.1058,    157.7139;
 45,    44.1999,    131.9773;
 46,    53.94872,   84.81891;
 47,    -28.62112,  25.23523;
 48,    35.80523,   9.32302;
 49,    -37.08995,  -12.33192;
 50,    -7.9377,    -14.3752;
 51,    32.3615,    -64.6987;
 52,    -7.3777,    72.4844;
 53,    64.875,     -147.861;
 55,    -77.731,    167.5881;
 56,    48.264,     -117.1257;
 57,    33.6058,    -116.4532;
 59,    19.592,     -155.8931];

stnNames = char('I02AR','I03AU','I04AU','I05AU','I06AU','I07AU','I08BO','I09BR','I10CA', 'I11CV','I13CL','I14CL','I17CI','I18DK','I19DJ','I21FR','I22FR','I23FR','I24FR','I26DE','I27DE','I30JP','I31KZ','I32KE','I33MG','I34MN','I35NA','I36NZ','I37NO','I39PW','I41PY','I42PT','I43RU','I44RU','I45RU','I46RU','I47ZA','I48TN','I49GB','I50GB','I51GB','I52GB','I53US','I55US','I56US','I57US','I59US');
mapNames = char(' I02AR',' I03AU',' I04AU',' I05AU',' I06AU',' I07AU',' I08BO',' I09BR',' I10CA', ' I11CV',' I13CL',' I14CL',' I17CI',' I18DK',' I19DJ',' I21FR',' I22FR',' I23FR',' I24FR',' I26DE',' I27DE',' I30JP',' I31KZ',' I32KE',' I33MG',' I34MN',' I35NA',' I36NZ',' I37NO',' I39PW',' I41PY',' I42PT',' I43RU',' I44RU',' I45RU',' I46RU',' I47ZA',' I48TN',' I49GB',' I50GB',' I51GB',' I52GB',' I53US',' I55US',' I56US',' I57US',' I59US');

%===================================================================================================
function secs = date2unixsecs(varargin)
   nargsin = nargin;
   error(nargchk(0, 6, nargsin));
   if nargsin
      argv = { 1 1 1 0 0 0 };
      argv(1:nargsin) = varargin;
   else
      argv = num2cell(clock);
   end
   [year, month, day, hour, minute, second] = deal(argv{:});
   secs = 86400 * (date2jd(year,month,day,hour,minute,second) - date2jd(1970, 1, 1));
    
%===================================================================================================   
   function [year, month, day, hour, minute, second] = unixsecs2date(secs)
   error(nargchk(1, 1, nargin));
   [year, month, day, hour, minute, second] ...
      = jd2date(secs / 86400 + date2jd(1970, 1, 1));
   
%===================================================================================================  
  function jd = date2jd(varargin)
   nargsin = nargin;
   error(nargchk(0, 6, nargsin));
   if nargsin
      argv = {1 1 1 0 0 0};
      argv(1:nargsin) = varargin;
   else
      argv = num2cell(clock);
   end
   [year, month, day, hour, minute, second] = deal(argv{:});
   a = floor((14 - month)/12);
   y = year + 4800 - a;
   m = month + 12*a - 3;
   jd = day + floor((153*m + 2)/5) ...
        + y*365 + floor(y/4) - floor(y/100) + floor(y/400) - 32045 ...
        + ( second + 60*minute + 3600*(hour - 12) )/86400;
    
%===================================================================================================    
 function [year, month, day, hour, minute, second] = jd2date(jd)
   nargsin = nargin;
   error(nargchk(1, 1, nargsin));
   ijd = floor(jd + 0.5);               % integer part
   if nargout > 3
      fjd = jd - ijd + 0.5;             % fraction part
      [hour, minute, second] = days2hms(fjd);
   end
   a = ijd + 32044;
   b = floor((4 * a + 3) / 146097);
   c = a - floor((b * 146097) / 4);
   d = floor((4 * c + 3) / 1461);
   e = c - floor((1461 * d) / 4);
   m = floor((5 * e + 2) / 153);
   day   = e - floor((153 * m + 2) / 5) + 1;
   month = m + 3 - 12 * floor(m / 10);
   year  = b * 100 + d - 4800 + floor(m / 10);

%==========================================================================   
function [hour, minute, second] = days2hms(days)
   error(nargchk(1, 1, nargin));
   second = 86400 * days;
   hour   = fix(second/3600);           % get number of hours
   second = second - 3600*hour;         % remove the hours
   minute = fix(second/60);             % get number of minutes
   second = second - 60*minute;         % remove the minutes
%
%==========================================================================
function hms = sec2hms(sec)

%SEC2HMS Converts time from seconds to hrs:min:sec vector format
%  hms = SEC2HMS(sec) converts time from seconds to hrs:min:sec
%  vector format.
%  See also HMS2SEC, SEC2HR, MAT2HMS, HMS2MAT, TIMEDIM, TIME2STR

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:48:06 $
if nargin==0;   error('Incorrect number of arguments');   end
%  Compute the time in hms by first transforming from sec to hrs.

hms = hr2hms(sec2hr(sec));
%
%==========================================================================
function hrs = sec2hr(seconds)

% SEC2HR Converts time from seconds to hours
%  hr = SEC2HR(sec) converts time from seconds to hours.%
%  See also HR2SEC, SEC2HMS, TIMEDIM, TIME2STR
%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Brown, E. Byrns
%  $Revision: 1.7 $    $Date: 1998/08/10 17:48:06 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(seconds)
     warning('Imaginary parts of complex TIME argument ignored')
     seconds = real(seconds);
end

hrs=seconds/3600;
%
%==========================================================================
function hms=hr2hms(hrs)

% HR2HMS Converts time from hours to hrs:min:sec vector format
%  hms = HR2HMS(hr) converts time from hours to hrs:min:sec
%  vector format.
%  See also HMS2HR,  HR2SEC, MAT2HMS, HMS2MAT, TIMEDIM, TIME2STR
%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%  $Revision: 1.7 $    $Date: 1998/08/10 17:47:49 $

if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(hrs)
     warning('Imaginary parts of complex TIME argument ignored')
     hrs = real(hrs);
end

%  Test for empty inputs
if isempty(hrs);     hms = [];   return;   end

%  Construct a sign vector which has +1 when hrs >= 0 and -1 when hrs < 0.
signvec = sign(hrs);
signvec = signvec + (signvec == 0);    %  Enforce +1 when hrs == 0

%  Compute the hours, minutes and seconds
hrs = abs(hrs);           %  Work in absolute value.  Signvec will set sign later
h   = fix(hrs);           %  Hours
ms  = 60*(hrs - h);       %  Minutes and seconds
m   = fix(ms);            %  Minutes
s   = 60*(ms - m);        %  Seconds

%  Determine where to store the sign of the time.  It should be
%  associated with the largest nonzero component of h:m:s.
hsign = signvec .* (h~=0);                %  Associate with hours
msign = signvec .* (h==0 & m~=0);         %  Assoicate with minutes (h = 0)
ssign = signvec .* (h==0 & m==0 & s~=0);  %  Associate with seconds (h = m = 0)

%  In the application of signs below, the ~ operator is used so that
%  the sign vector contains only +1 and -1.  Any zero occurances causes
%  data to be lost when the sign has been applied to a higher component
%  of h:m:s.
h = (~hsign + hsign).*h;      %  Apply signs to the hours
m = (~msign + msign).*m;      %  Apply signs to minutes
s = (~ssign + ssign).*s;      %  Apply signs to seconds

hms = mat2hms(h,m,s);     %  Construct the hms vector for output
%
%==========================================================================
function hmsvec = mat2hms(h,m,s)

%MAT2HMS Converts a [hrs min sec] matrix to vector format
%
%  hms = MAT2HMS(h,m,s) converts a hrs:min:sec matrix into a vector
%  format.  The vector format is hms = 100*hrs + min + sec/100.
%  This allows h,m,s triple to be compressed into a single value,
%  which can then be employed similar to a second or hour vector.
%  The inputs h, m and s must be of equal size.  Minutes and
%  second must be between 0 and 60.
%  hms = MAT2HMS(mat) assumes and input matrix of [h m s].  This is
%  useful only for single column vector for h, m and s.%
%  hms = MAT2HMS(h,m) and hms = MAT2HMS([h m]) assume that seconds
%  are zero, s = 0.%
%  hms = MAT2HMS(h,m,s,n) uses n as the accuracy of the seconds
%  calculation.  n = -2 uses accuracy in the hundredths position,
%  n = 0 uses accuracy in the units position.  Default is n = -5.
%  For further discussion of the input n, see ROUNDN.%
%  See also HMS2MAT
%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%  $Revision: 1.8 $    $Date: 1998/08/10 17:47:55 $

if nargin == 0
   error('Incorrect number of arguments')

elseif nargin==1
   if size(h,2)== 3
       s = h(:,3);   m = h(:,2);   h = h(:,1);
   elseif size(h,2)== 2
       m = h(:,2);   h = h(:,1);   s = zeros(size(h));
   elseif size(h,2) == 0
       h = [];   m = [];   s = [];
   else
       error('Single input matrices must be n-by-2 or n-by-3.');
   end
   n = -5;

elseif nargin == 2
   s = zeros(size(h));
   n = -5;

elseif nargin == 3
   n = -5;
end

%  Test for empty arguments
if isempty(h) & isempty(m) & isempty(s);  hmsvec = [];  return;  end

%  Don't let seconds be rounded beyond the tens place.
%  If you did, then 55 seconds rounds to 100, which is not good.
if n == 2;  n = 1;   end

%  Complex argument tests
if any([~isreal(h) ~isreal(m) ~isreal(s)])
    warning('Imaginary parts of complex TIME argument ignored')
	h = real(h);   m = real(m);   s = real(s);
end

%  Dimension and value tests
if  ~isequal(size(h),size(m),size(s))
    error('Inconsistent dimensions for input arguments')
elseif any(rem(h(~isnan(h)),1) ~= 0 | rem(m(~isnan(m)),1) ~= 0)
    error('Hours and minutes must be integers')
end

if any(abs(m) > 60) | any (abs(m) < 0)       %  Actually algorithm allows for
    error('Minutes must be >= 0 and < 60')   %  up to exactly 60 seconds or
                                             %  60 minutes, but the error message
elseif any(abs(s) > 60) | any(abs(s) < 0)    %  doesn't reflect this so that angst
    error('Seconds must be >= 0 and < 60')   %  is minimized in the user docs
end

%  Ensure that only one negative sign is present
if any((s<0 & m<0) | (s<0 & h<0) | (m<0 & h<0) )
    error('Multiple negative entries in a hms specification')
elseif any((s<0 & (m~=0 | h~= 0)) | (m<0 & h~=0))
    error('Incorrect negative HMS specification')
end

%  Construct a sign vector which has +1 when
%  time >= 0 and -1 when time < 0.  Note that the sign of the
%  time is associated with the largest nonzero component of h:m:s
negvec = (h<0) | (m<0) | (s<0);
signvec = ~negvec - negvec;

%  Convert to all positive numbers.  Allows for easier
%  adjusting at 60 seconds and 60 minutes
h = abs(h);     m = abs(m);    s = abs(s);

%  Truncate seconds to a specified accuracy to eliminate round-off errors
[s,msg] = roundn(s,n);
if ~isempty(msg);   error(msg);   end

%  Adjust for 60 seconds or 60 minutes. If s > 60, this can only be
%  from round-off during roundn since s > 60 is already tested above.
%  This round-off effect has happened though.
indx = find(s >= 60);
if ~isempty(indx);   m(indx) = m(indx) + 1;   s(indx) = 0;   end

%  The user can not put minutes > 60 as input.  However, the line
%  above may create minutes > 60 (since the user can put in m == 60),
%  thus, the test below includes the greater than condition.
indx = find(m >= 60);
if ~isempty(indx);   h(indx) = h(indx) + 1;   m(indx) = m(indx)-60;   end

%  Construct the hms vector format
hmsvec = signvec .* (100*h + m + s/100);
%
%==========================================================================
function [hout,mout,sout] = hms2mat(hms,n)

%HMS2MAT Converts a hms vector format to a [hrs min sec] matrix
%
%  [h,m,s] = HMS2MAT(hms) converts a hms vector format to a
%  hrs:min:sec matrix.  The vector format is hms = 100*hrs + min + sec/100.
%  This allows compressed hms data to be expanded to a h,m,s triple,
%  for easier reporting and viewing of the data.
%
%  [h,m,s] = HMS2MAT(hms,n) uses n digits in the accuracy of the
%  seconds calculation.  n = -2 uses accuracy in the hundredths position,
%  n = 0 uses accuracy in the units position.  Default is n = -5.
%  For further discussion of the input n, see ROUNDN.
%
%  mat = HMS2MAT(...) returns a single output argument of mat = [h m s].
%  This is useful only if the input hms is a single column vector.
%
%       See also MAT2HMS

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%  $Revision: 1.7 $    $Date: 1998/08/10 17:47:47 $

if nargin == 0
     error('Incorrect number of arguments')
elseif nargin == 1
     n = -5;
end

%  Test for empty arguments

if isempty(hms); hout = []; mout = []; sout = []; return; end

%  Test for complex arguments

if ~isreal(hms)
     warning('Imaginary parts of complex TIME argument ignored')
	 hms = real(hms);
end

%  Don't let seconds be rounded beyond the tens place.
%  If you did, then 55 seconds rounds to 100, which is not good.

if n == 2;  n = 1;   end

%  Construct a sign vector which has +1 when hms >= 0 and -1 when hms < 0.

signvec = sign(hms);
signvec = signvec + (signvec == 0);   %  Ensure +1 when hms = 0

%  Decompress the hms data vector

hms = abs(hms);
h = fix(hms/100);                      %  Hours
m = fix(hms) - abs(100*h);             %  Minutes
[s,msg] = roundn(100*rem(hms,1),n);    %  Seconds:  Truncate to roundoff error
if ~isempty(msg);   error(msg);   end

%  Adjust for 60 seconds or 60 minutes.
%  Test for seconds > 60 to allow for round-off from roundn,
%  Test for minutes > 60 as a ripple effect from seconds > 60


indx = find(s >= 60);
if ~isempty(indx);   m(indx) = m(indx) + 1;   s(indx) = s(indx) - 60;   end
indx = find(m >= 60);
if ~isempty(indx);   h(indx) = h(indx) + 1;   m(indx) =  m(indx) - 60;   end

%  Data consistency checks

if any(m > 59) | any (m < 0)
    error('Minutes must be >= 0 and <= 59')

elseif any(s >= 60) | any( s < 0)
    error('Seconds must be >= 0 and < 60')
end

%  Determine where to store the sign of the time.  It should be
%  associated with the largest nonzero component of h:m:s.

hsign = signvec .* (h~=0);
msign = signvec .* (h==0 & m~=0);
ssign = signvec .* (h==0 & m==0 & s~=0);

%  In the application of signs below, the ~ operator is used so that
%  the sign vector contains only +1 and -1.  Any zero occurances causes
%  data to be lost when the sign has been applied to a higher component
%  of h:m:s.  Use fix function to eliminate potential round-off errors.

h = (~hsign + hsign).*fix(h);      %  Apply signs to the hours
m = (~msign + msign).*fix(m);      %  Apply signs to minutes
s = (~ssign + ssign).*s;           %  Apply signs to seconds


%  Set the output arguments

if nargout <= 1
    hout = [h m s];
elseif nargout == 3
    hout = h;   mout = m;   sout = s;
else
    error('Invalid number of output arguments')
end

%
%==========================================================================
function [x,msg] = roundn(x,n)

%ROUNDN  Rounds input data at specified power of 10
%
%  y = ROUNDN(x) rounds the input data x to the nearest hundredth.
%
%  y = ROUNDN(x,n) rounds the input data x at the specified power
%  of tens position.  For example, n = -2 rounds the input data to
%  the 10E-2 (hundredths) position.
%
%  [y,msg] = ROUNDN(...) returns the text of any error condition
%  encountered in the output variable msg.
%
%  See also ROUND

%  Copyright 1996-1998 by Systems Planning and Analysis, Inc. and The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1.7 $    $Date: 1998/08/10 17:48:04 $

msg = [];   %  Initialize output

if nargin == 0
    error('Incorrect number of arguments')
elseif nargin == 1
    n = -2;
end

%  Test for scalar n

if max(size(n)) ~= 1
   msg = 'Scalar accuracy required';
   if nargout < 2;  error(msg);  end
   return
elseif ~isreal(n)
   warning('Imaginary part of complex N argument ignored')
   n = real(n);
end

%  Compute the exponential factors for rounding at specified
%  power of 10.  Ensure that n is an integer.

factors  = 10 ^ (fix(-n));

%  Set the significant digits for the input data

x = round(x * factors) / factors;
%
%==========================================================================