function s=NOAAweighted(sig,fs,filtertype)
%function s=NOAAweighted(sig,fs,filtertype)
%Filters input signal sig according to recommendations from NOAA/NMFS:
%National Marine Fisheries Service (2016). Technical guidance for assessing
%the effects of anthropogenic sound on marine mammal hearing underwater
%acoustic thresholds for onset of permanent and temporary threshold shifts.
%NOAA Technical Memorandum NMFS-OPR-55. 178 pp.
%
% The filter to be used is specified by filtertype:
%'HF'        :HF-cetacean
%'MF'        :MF-cetacean
%'LF'        :LF-cetacean
%'Otariid'   :Otariids
%'Phocid'    :Phocids
%fs is sample rate of sig (in Hz)
%
%Detailed description of the code can be found in Tougaard and Beedholm
%(2018), "A practical implementation of auditory time and frequency weighting in marine bioacoustics"
%Applied Acoustics.
%
%Aarhus University, Kristian Beedholm and Jakob Tougaard, 2018
%Shared under Creative Commons license CC BY-SA 4.0 (share alike).

if nargin<3
    error('not enough input arguments')
end

%Flip to vertical if input signal is a horizontal array
horizontal=false;
if size(sig,1)==1 
    sig=sig';
    horizontal=true; %Remember input was horizontal
end

%set filter parameters. These parameters come from table ES2 in NMFS(2016).
%Custom filters with custom parameters can be added to code, if required.
%See NMFS(2016) for details regarding interpretations of filter parameters
switch upper(filtertype(1))   
    case 'H'    %cetHF
        a=1.8; b=2; f1=12000; f2=140000; C=1.36;
    case 'M'    %cetMF
        a=1.6; b=2; f1=8800; f2=110000; C=1.3;
    case 'L'    %cetLF
        a=1; b=2; f1=200; f2=19000; C=0.13;
    case 'O'    %Otariid
        a=1; b=2; f1=1900; f2=30000; C=0.75;
    case 'P'    %Phocid
        a=2; b=2; f1=940; f2=45000; C=0.64;
    otherwise
        error('Filter type not recognized')
end

L=size(sig,1);                  %Impulse response length equal to signal length  
f=linspace(0,fs/2,ceil(L/2))';  %Define frequency axis
W=sqrt((f/f1).^(2*a)./((1+(f/f1).^2).^a.*(1+(f/f2).^2).^b)*10^(C/20));         %Create filter response 
W=[W ; flipud(W(1+mod(L,2):end))];     %two-sided spectrum with zero phase
s=real(ifft(W.*fft(sig)));                   %convolution of signal with impulse response

%If input was a horizontal array, then flip back to horizontal before returning
if horizontal  
    s=s';
end
