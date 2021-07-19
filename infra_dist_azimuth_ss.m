function [] = infra_dist_azimuth_ss()
%
% This function calculates the distances along a great circle and theoretical 
% delay times of infrasonic signals (source-to-receiver) at any 
% given location on Earth. The function is written for a single source and
% a single receiver. Please refer to infra_dist_azimuth.m for arrivals at
% all stations of the CTBTO IMS network.
%
% This function also accounts for travel times associated with various 
% propagation channels (tropospheric, stratospheric and thermospheric).
% The user has a choice to output a map, and print a text file with arrival
% times for all three propagation channels, distance in km and degrees, and 
% back azimuth. 
% 
% input variables: (date, time, source, stn)
%       date [dd mm yyyy]
%       time [hh mm ss]
%       source [lat lon]
%       station [lat lon]
% Last updated: July 2021
% v2.7.1
%
% Elizabeth A. SILBER (c)
% e-mail: esilber@uwo.ca
% The University of Western Ontario
%
% ========================================================================
evdate = input('Enter the event date [dd mm yyyy]:  ');
time = input ('Enter the event time [hh mm ss]: ');
source = input ('Enter the coordinates of the source [lat lon]: ');
stn = input ('Enter the coordinates of the receiver [lat lon]: ');

eventTime  = date2unixsecs(evdate(3),evdate(2),evdate(1),time(1),time(2),time(3));

k = input('Do you wish to include delay time in hh:mm:ss and ETA? 0 = no, 1 = yes    ');
map = input('Do you wish to print a map? 0 = no, 1 = yes   '); 

% Produce a .txt file with results
file = fopen(['DazD-',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'.txt'],'w');
fprintf(file,('Date of the event [dd-mm-yyyy]: %02.0f-%02.0f-%4.0f\n'),evdate(1),evdate(2),evdate(3));     
fprintf(file,('Time of the event [hh:mm:ss]:   %02.0f:%02.0f:%02.0f   UTC\n'),time(1),time(2),time(3));
fprintf(file,('Location of the event [lat(N), lon(E)]:   %02.2f, %02.2f% \n'),source(1),source(2));
fprintf(file,'\n================================================================================\n');
fprintf(file,'STN \tDist  \tBack Az \t Delay(s) \t Strato Del \t Thermo Del\n');
fprintf(file,'    \t(km)  \t(deg)   \t @340 m/s \t@ 285 m/s   \t @220 m/s\n');
fprintf(file,'================================================================================\n');


rangdeg = distance (source(1),source(2),stn(1),stn(2),'deg');
dist = deg2km(rangdeg); %great circle range in km
az = azimuth(stn(1),stn(2),source(1),source(2),'deg'); %back azimuth

% Tropospheric (direct) arrival
delay(1) = dist/0.340; %delay time in seconds
temp = sec2hms(delay(1)); time1(1:3) = hms2mat(temp);
eta0(1) = eventTime + delay(1); %arrival in unix time
[etatr(1),etatr(2),etatr(3),etatr(4),etatr(5),etatr(6)] = unixsecs2date(eta0(1)); % arival in actual time (date, time)

% Stratospheric arrival
delay(2) = dist/0.285;
temp = sec2hms(delay(2)); time2(1:3) = hms2mat(temp);
eta0(2) = eventTime + delay(2);
[etastr(1),etastr(2),etastr(3),etastr(4),etastr(5),etastr(6)] = unixsecs2date(eta0(2));

% Thermospheric arrival
delay(3) = dist/0.220;
temp = sec2hms(delay(3)); time3(1:3) = hms2mat(temp);
eta0(3) = eventTime + delay(3);
[etath(1),etath(2),etath(3),etath(4),etath(5),etath(6)] = unixsecs2date(eta0(3));


fprintf(file,'\t%5.1f\t\t%3.1f\t\t%8.1f\t%8.1f\t%8.1f\n',dist,az,delay(1),delay(2),delay(3));
if k == 1,
    fprintf(file,'      \t Delay [hh:mm:ss]  \t\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\n',time1(1),time1(2),time1(3),time2(1),time2(2),time2(3),time3(1),time3(2),time3(3));
    fprintf(file,'      \t ETA   [hh:mm:ss]  \t\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\t%02.0f:%02.0f:%02.0f\n',etatr(4),etatr(5),etatr(6),etastr(4),etastr(5),etastr(6),etath(4),etath(5),etath(6));
end


if map == 1,   
    infra_map_world(source, evdate, stn);
    
end
close all
disp('Run completed.');

%==========================================================================
function[] = infra_map_world(source, evdate, stn)

figure1 = figure;
worldmap world  % to plot whole world
whos -file coast.mat
load coast      % load coastline

plotm(lat, long) % plot coastline on the map
geoshow('landareas.shp', 'FaceColor', [0.35 0.6 0.35]);
geoshow('worldlakes.shp', 'FaceColor', 'cyan');    
    
% make a map
hold all;

plotm(source(1), source(2), 'ko', 'MarkerFaceColor','y','MarkerSize',20);
textm(source(1), source(2), 'Source','FontSize',14);

plotm(stn(1), stn(2), 'dk', 'MarkerSize', 8, 'MarkerFaceColor', [0.8500 0.2250 0.0980]);
textm(stn(1), stn(2), 'Receiver','FontSize',14); 

set(get(groot, 'Children'), 'WindowState', 'maximized');

title(['Event of ',num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3))],'FontSize',16);
print('-djpeg','-r300',[num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'-WorldMap']);  
saveas(figure1,[num2str(evdate(1)),'-',num2str(evdate(2)),'-',num2str(evdate(3)),'-WorldMap.fig']);
%
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