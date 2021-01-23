clc
clear
close all
SampFreq = 256;
addpath('D:\tfsa_5-5\windows\win64_bin');

t = 0:1/SampFreq:1-1/SampFreq;

addpath('D:\TFSA7\TFSA7');%


Sig1 = 1*exp(1i*(-0.5*pi*(1*SampFreq*t.^2))+1i*(pi*(SampFreq*t))); %300tªÚ’ﬂ150t



Sig5 = 1*exp(1i*(pi*(0.05*SampFreq*t +0.7*SampFreq*t.^3/3)));
Sig6 = 1*exp(1i*(pi*(0.3*SampFreq*t +0.7*SampFreq*t.^3/3)));


num=3;
Sig =1*Sig1 +1*Sig5+1*Sig6;
SigA=Sig;
IF_O(:,1)=0.05*(1*SampFreq)/2+0.7*SampFreq*t.^2/2;
IF_O(:,2)=0.3*(1*SampFreq)/2+0.7*SampFreq*t.^2/2;
IF_O(:,3)=-0.5*(1*SampFreq)*t+SampFreq/2;

IF_O=IF_O/(SampFreq/2);


% HADTFD BASED

iiii=0;
dis=0;
%Sig =1*Sig1+Sig3;


NS=100;
t=16:128-15;

iiii=0;



for snr=0:5:30
    iiii=iiii+1;
    
    for k1=1:NS
        Sig =SigA;
        
        Sig=awgn(Sig,snr,'measured');
        tfd = cmpt( Sig, 'ckd', 1, 0.2,0.2, SampFreq);
        
        tfd(tfd<0)=0;
        orient=orient_img_calc(tfd ,1, 1, 1);
        orient=180*(orient)/(pi);
        
        %   tfsapl(Sig,tfd);
        for kk=0:1
            if kk==0
                [fmult,out, peaks] = component_linking_new(tfd,orient,0.05,5,5);
                [fmult]= merge_IFs(fmult,orient,30,50,length(Sig)/2);%[IF,out, peaks] = component_linking_neww_direc(tfd,orient,0.2,length(Sig)/4,30,5);
                
            else
                [fmult,~, ~] = component_linking(tfd,0.1,64);
            end
            findex=fill_zeros(fmult);
            findex1=zeros(num,length(Sig));
            [aa,~]=size(findex);
            
            for ii=1:aa
                findex1(ii,1:length(findex))=findex(ii,:);
                
            end
            findex=findex1;
            msee=0.1*ones(1,num);
            dis=0;
            [aa,~]=size(findex);
            %findex=findex/(SampFreq/2);
            
            for ii=1:num
                if ii<=aa
                    
                    IF=findex(ii,:)/(length(Sig));
                    %t=t(5:end-5);
                    for i=1:num
                        c(i)=sum(abs(IF(t).'-IF_O(t,i)).^2);
                    end
                    [a1 b1]=min(c);
                    if msee(b1)>=a1(1)/length(t)
                        msee(b1)=a1(1)/length(t);
                    end
                    if dis==1
                        figure;
                        plot(t,IF(t),'-',t,IF_O(t,b1),'d');
                    end
                end
            end
            msee1(kk+1,k1)=mean(msee);
        end
    end
    mse_adtfd_new(iiii)=mean(msee1(1,:));
    mse_adtfd_old(iiii)=mean(msee1(2,:));
end



snr=0:5:30;
plot(snr,10*(log10(mse_adtfd_new)),'-.b+','linewidth',4);
hold on;
plot(snr,10*(log10(mse_adtfd_old)),'--md','linewidth',4);
xlabel('Signal to noise ratio');
ylabel('Mean square error in dB');