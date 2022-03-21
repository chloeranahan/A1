%UPDATES to A1 for recovered marks:

%1. Fixed the avg temperature plots by changing Vx to
Vx = sqrt((V^2)-(Vy.^2));
%so that the overall velocity was vth and therefore the overall temperature
%was 300 K (see updated Q1 and Q2 temperature plots over time
    %lost 3 marks here

    
%2. Added the units for the calculated mfp (m) and tmn (s) in the report
    %lost 2 marks here
    

%3. Fixed the standard deviattion and histogram by changing std dev to
std = vth/(sqrt(2));
%see updated histogram plot in report
    %lost 2 marks here