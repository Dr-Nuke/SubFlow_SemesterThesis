% Create signal.
       Fs = 2500;   t = 0:1/Fs:.296;
       x = D;  % A cosine of 200Hz plus noise

       % Construct a Welch spectrum object.
       h = spectrum.welch('hamming',64);

       % Create three power spectral density estimates.
       hpsd1 = psd(h,x,'Fs',Fs);
       Pxx1 = hpsd1.Data;
       W = hpsd1.Frequencies;

       h.WindowName = 'Kaiser';
       hpsd2 = psd(h,x,'Fs',Fs);
       Pxx2 = hpsd2.Data;

       h.WindowName = 'Chebyshev';
       hpsd3 = psd(h,x,'Fs',Fs);
       Pxx3 = hpsd3.Data;

       % Instantiate a PSD data object and store the three different
       % estimates since they all share the same frequency information.
       hpsd = dspdata.psd([Pxx1, Pxx2, Pxx3],W,'Fs',Fs)