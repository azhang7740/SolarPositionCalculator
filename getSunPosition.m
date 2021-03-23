% Angelina Zhang 7.31.20
% Function calculates azimuth and elevation given coordinates, a date, and
% time. Calculations implemented from Astronomical Algorithms by Jean Meeus
% and the NOAA solar position calculator. 
% azimuth is the angle between the sun's position and true north (eastward
% is positive
% elevation is the angle between the sun's position and the horizontal
% plane
% sunTime is a matrix [sunlight duration, sunrise, solar noon, sunset]
% lat, lon are the latitude and longitude in degrees
% date is a matrix [year month day]
% time is a matrix [hour minute second] (military time)
function [azimuth, elevation, sunTime] = getSunPosition(lat, lon, date, time, timeZone)
dateAndTime = datetime(date(1,1), date(1,2), date(1,3), 'TimeZone', timeZone);
[utcAdd, localDayTime] = getTime(time, dateAndTime);
julianCent = getJulian(dateAndTime, localDayTime, utcAdd);

%Geometric Mean Longitude Sun/Mean Equinox date (degrees)
mlSun = mod(280.46646+julianCent*(36000.76983 + julianCent*0.0003032), 360);
%Geometric Mean Anomaly Sun (degrees)
maSun = 357.52911+julianCent*(35999.05029 - 0.0001537*julianCent);
%Eccentricity Earth's Orbit
eccEarthOrbit = 0.016708634-julianCent*(0.000042037+0.0000001267*julianCent);

%Sun Equation of Center
eqCenter = sind(maSun)*(1.914602-julianCent*(0.004817+0.000014*julianCent))...
    +sind(2*maSun)*(0.019993-0.000101*julianCent)+sind(3*maSun)*0.000289;

%Sun True Longitude
longSun = mlSun+eqCenter;
%Sun True Anomaly
anomSun = maSun + eqCenter;

% Sun Radius Vector (AUs) (distance between sun and earth centers)
sunRad = (1.000001018*(1-eccEarthOrbit*eccEarthOrbit))/(1+eccEarthOrbit*cosd(anomSun));

% Sun Apparent Longitude (degrees)
appLong = longSun-0.00569-0.00478*sind((125.04-1934.136*julianCent));

% Mean Obliquity of Ecliptic (degrees)
meanObliquity = 23+(26+((21.448-julianCent*(46.815+julianCent*(0.00059-julianCent*0.001813))))/60)/60;
% Corrected Obliquity(degrees)
corrObliq = meanObliquity + 0.00256*cosd((125.04-1934.136*julianCent));

%Sun Rt Ascen (degrees)
sunAscen = atan2d(cosd(appLong), cosd(corrObliq)*sind(appLong));
%Sun Declination (degrees)
sunDeclin = asind(sind(corrObliq)*sind(appLong));

%var y
varY = (tand(corrObliq/2))^2;

%Equation of Time (minutes)
eqTime = 4*(varY*sind(2*mlSun)-2*eccEarthOrbit*sind(maSun)+4*eccEarthOrbit*varY*sind(maSun)*...
    cosd(2*mlSun)-0.5*varY*varY*sind(4*mlSun)-1.25*eccEarthOrbit*eccEarthOrbit*sind(2*maSun))*180/pi;
%HA Sunrise (degrees)
%(acos(cosd(90.833)/(cosd(lat)*cosd(sunDeclin))-tand(lat)*tand(sunDeclin))*180/pi
haSunrise = (acos(cosd(90.833)/(cosd(lat)*cosd(sunDeclin))-tand(lat)*tand(sunDeclin)))*180/pi;
sunTime = getSunTime(lon, eqTime, utcAdd, haSunrise);
   
%True Solar Time (minutes)
tSolTime = mod(localDayTime*1440+eqTime+4*lon-60*utcAdd,1440);

if tSolTime/4 < 0 
    hourAngle = tSolTime/4 + 180;
else
    hourAngle = tSolTime/4 - 180;
end

%Solar Zenith Angle
solZenith = (acos(sind(lat)*sind(sunDeclin)+cosd(lat)*cosd(sunDeclin)*cosd(hourAngle))) * 180/pi;
%Solar Elevation Angle
solEl = 90 - solZenith;

elevation = getElevation(solEl);
azimuth = getAzimuth(hourAngle, lat, solZenith, sunDeclin);
end

% Changes from local time to UTC time, with the day variable as a decimal
% that accounts for the time in the day
function [utcAdd, localDayTime] = getTime(time, dateAndTime)
utcAdd = hours(tzoffset(dateAndTime));
localDayTime = time(1,1)/24 + time(1,2)/1440 + time(1,3)/86400;
end

function [julianCentury] = getJulian(date, localDayTime, utcAdd)
% baseDate = datetime(1900, 1, 1);
% numDays = daysact(baseDate - 2, date);
% julianDay = numDays + 2415018.5 + localDayTime - utcAdd/24
utcTime = localDayTime - utcAdd/24;
julianDay = getJulianDay(year(date), month(date), day(date)) + utcTime;
julianCentury = (julianDay - 2451545) / 36525;
end

% Function that calculates the Julian Day with the current year, month, and
% day
function julianDay = getJulianDay(Y, M, D)
if M <= 2
    Y = Y - 1;
    M = M + 12;
end
A = fix(Y/100);
B = 2 - A + fix(A/4);
julianDay = fix(365.25*(Y+4716))+fix(30.6001*(M+1))+D+B-1524.5;
end

% Returns the apparent solar elevation by applying the calculations that
% determine the approximate atmospheric refraction
function el = getElevation(solEl)
if solEl > 85
    approxAR = 0;
elseif solEl > 5
    approxAR = 58.1/tand(solEl)-0.07/(tand(solEl)^3)+0.000086/(tand(solEl)^5);
elseif solEl > -0.575
    approxAR = 1735 + solEl *(-518.2+solEl*(103.4+solEl*(-12.79+solEl*0.711)));
else
    approxAR = -20.772/tand(solEl);
end
approxAR = approxAR / 3600;
el = solEl + approxAR;
end

% Calculates and corrects the azimuth angle so it is between 0 and 360
% degrees
function az = getAzimuth(hourAngle, lat, solZenith, sunDeclin)
temp = (180/pi)*(acos(((sind(lat)*cosd(solZenith))-sind(sunDeclin))/(cosd(lat)*sind(solZenith))));
if hourAngle>0
   az = mod(temp + 180,360);
else
   az = mod(540 - temp ,360);
end
end

% Calculates micellaneous values from NOAA's solar position calculator,
% such as the duration of the sun in the sky during a day, the local time
% of sunrise and sunset, and the time of solar noon
function [sunTime] = getSunTime(lon, eqTime, utcAdd, haSunrise)
%solar noon (LST)
solarNoon = (720-4*lon-eqTime+ utcAdd*60)/1440;
%Sunrise Time (fraction of a day in local standard time)
sunriseTime = solarNoon - haSunrise*4/1440;
%Sunset Time (fraction of a day in local standard time)
sunsetTime = solarNoon + haSunrise*4/1440;
%Sunlight Duration (minutes)
durSun = 8 * haSunrise;

sunTime = [durSun, sunriseTime, solarNoon, sunsetTime];
end