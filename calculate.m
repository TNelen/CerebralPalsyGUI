function [Y1,Y2] = calculate(inputX,inputY,Windows)

Fs = 200;                   %Sample frequency
L1 = length(inputX);        %Data of the healthy subject
L2 = length(inputY);        %Data of the subject with cerebral palsy, has to be equal to L1
n = 2^nextpow2(L1);         %Calculate next power of 2 after the length of the input used for fft.
win = 'No window';          %Default = no windowing
switch Windows              %Choose the appropriate windowing function.
    case 1
        w=rectwin(L1);      
        win = 'Rectangular window';
    case 2
        w=bartlett(L1);
        win = 'Bartlett window';
    case 3
        w=hann(L1);
        win = 'Hann window';
    case 4
        w=hamming(L1);
        win = 'Hamming window';
    case 5
        w=blackman(L1);
        win = 'Blackman window';
    otherwise
        w=1;                %No windowing
end

Y1 = fft(w.*inputX,n);      %Calculate the fft of the input data multiplied 
Y2 = fft(w.*inputY,n);      % by the window. Zeropadded to a length of n

Y1 = abs(Y1/L1);            %Take the absolute value of the fft-data devided by the length
Y2 = abs(Y2/L2);            %Take the amplitude(magnitude) of the fft. 
                            %Normalizing the values by dividing by the
                            %sample points(length)

Y1 = Y1(1:L1/2+1);          %Compute single-sided spectrum base on the 
Y2 = Y2(1:L2/2+1);          %Double sided spectrum and the even-valued 
                            %signal length L1 or L2

Y1(2:end-1) = 2*Y1(2:end-1); %Due to the power being divided between the 
Y2(2:end-1) = 2*Y2(2:end-1); %positive and negative side of the spectrum
                             %The negative part is lost by calculating the 
                             %onesided spectrum. This can be counteracted
                             %by multiplying by 2

val = Fs*(0:(L1/2))/L1;      %Calculate the amount of frequency bins for 
                             %plotting the graphs
                             
t1 = interp1([min(val), max(val)],[min(inputX), max(inputX)], val);
t2 = interp1([min(val), max(val)],[min(inputY), max(inputY)], val);
%The fft returns frequency bins but the data is in degrees. Therefore we
%change these values by interpolating them between the maximum and minimum
%amount of degrees. And we use these values for plotting.

figure(1);

%f = Fs*(0:(L1/2))/L1;
subplot(2,2,1);
plot((1:L1),inputX);
title('Healthy');
xlabel('sample point');
ylabel('Angle (deg)');
subplot(2,2,2);
plot((1:L2),inputY);
title('Cerebral Palsy');
xlabel('sample point');
ylabel('Angle (deg)');
subplot(2,2,3);
plot(t1,Y1);
title(strcat('Healthy: FFT (',win,')'));
xlabel('Angle (deg)');
ylabel('Amplitude');
subplot(2,2,4);
plot(t2,Y2);
title(strcat('Cerebral Palsy: FFT (',win,')'));
xlabel('Angle (deg)');
ylabel('Amplitude');


end


