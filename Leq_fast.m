function Leq=Leq_fast(sig,fs,tau,type)

%function Leq=Leq_fast(sig,fs,tau,type)
%Computes running rms-average of signal sig, sampled with sample rate sr
%tau is the integration time in seconds
%Default value for tau is 125 ms
%Kernel is either exponential decay (default, and type ='exponetial', or
%rectangular (type = 'rectangular')
%output is an array of the same size as the input. 

%Aarhus University, Jakob Tougaard and Kristian Beedholm, Jan. 2018 
%jat@bios.au.dk
%Shared under Creative Commons license CC BY-SA 4.0 (share alike)

if size(sig,1)==1    %Flip input if horizontal array
    sig=sig';
    horizontal=true;
else
    horizontal=false;
end

if nargin<3
    tau=0.125;  %default value 125 ms
end

if nargin<4    %default type is exponential
    type='e';
end

L=floor(fs*tau);    %Time constant in samples
if type(1)=='r'     %rectangular kernel
    if L>length(sig)
        error('Signal must be longer than tau')
    else
        w=[ones(L,1);zeros(length(sig)-L,1)];
    end
else                %exponential kernel
    if L*5>length(sig)
        error('Signal must be longer than five times tau')
    else
        w=exp(-(0:length(sig)-1)/(tau*fs))';  
    end
end

%convolution)
Leq=real(sqrt((ifft(fft(sig.^2).*fft(w)))/L));

if horizontal   %flip output back to horizontal if input is a horizontal array
    Leq=Leq';
end


