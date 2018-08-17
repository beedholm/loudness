function s=fileweighting(file,ip,filtertype,tau,type)
%Use weighting in time and frequency, implemented in functions 
%NOAAweighting and Leq_fast.
%Filename must be complete with extesion. The extension identifies the file
%format. Supported formats are WAV, OGG, FLAC, AU, MP3, and MPEG-4 AAC,
%It is strictly recommended to use only uncompressed or lossless compression formats.
%The lossy compression formats (MP3 and MPEG-4) will be processed correctly,
%but the lossy compression is highly likely to introduce erroneous and spurious
%artefacts in the output.
%Parameters filtertype, tau and type are passed on to NOAAweighting and
%Leq_fast
%Decimate result by factor ip to achieve a manageable array size in output

%Aarhus University, Kristian Beedholm and Jakob Tougaard, 2018
%Shared under creative commons licence CC BY-SA 4.0 (share alike)

chunkdur=10; %segment size in seconds. Overall compromize between time it takes to
%read data from file and time it takes to perform convolutions. Depending
%on factors such as speed of the computer, available RAM, and sample rate of
%signal, improved performance may be achieved by changing the segment size
%up or down.
fileinfo = audioinfo(file);
fs=fileinfo.SampleRate;
L=fileinfo.TotalSamples;
c=fileinfo.NumChannels;

if L/fs<60; %process in one piece if less than 1 minute
    chunk=audioread(file);
    s=myresample(real(Leq_fast(NOAAweighted(chunk,fs,filtertype),fs,tau,type)),ip);
else        %Read and process one chunk at a time 
    FL=chunkdur*fs; %segment size in samples
    if mod(FL,ip)>0
        error(['Downsample factor should be divisor of ',num2str(FL)])
    end
    %Create weighting windows for adding segments, based on 50% overlapping
    %Hann-windows
    win(:,2)=hanning(FL); %window for all segments, except first and last
    win(:,1)=[win(FL/2+1:end,2);zeros(FL/2,1)]; %window for first segment
    win(:,3)=flipud(win(:,1)); %window for last segment
    
    outarr=zeros(floor((L+FL)/ip),c); %initialize array for output.
    stop=floor(2*L/FL)+1;
    
    for k=1:stop
        if k==1 %First segment 
            chunk=audioread(file,[1, FL]).*repmat(win(:,1),1,c);  
            chunk=circshift(chunk,FL/2);
        elseif k==stop %last segment
            chunk=audioread(file,[(stop-3)*FL/2+1 (stop-3)*FL/2+FL]).*repmat(win(:,3),1,c);  
            chunk=circshift(chunk,FL/2);
        else %all other segments in between
            chunk=audioread(file,[(k-2)*FL/2+1 (k-2)*FL/2+FL]).*repmat(win(:,2),1,c);  
        end
        x=myresample(real(Leq_fast(NOAAweighted(chunk,fs,filtertype),fs,tau,type)),ip); %process current segment
        outarr((k-1)*FL/ip/2+1:(k-1)*FL/ip/2+length(x),:)=...
           outarr((k-1)*FL/ip/2+1:(k-1)*FL/ip/2+length(x),:)+x; %Add result to output array
    end
    outarr(floor((L+FL/2)/ip+1):end)=[]; %Truncate output to match file length
    outarr(1:FL/(ip*2))=[];
    s=outarr;
end

function [o]=myresample(s,ip)
%this function retains only every ip'th sample of the output. The simple
%resampling is prone to aliazing, but as the Leq-averaging is really a low-pass
%filter, additional anti-aliazing filtering is not required.

L=size(s,1);
idx=(1:floor(L/ip))*ip-ip+1;
o=s(idx,:);

