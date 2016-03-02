%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% cf_ampenv

clear all

% --------------------------------------------------------
% VERSION 2
% --------------------------------------------------------
v           = 2;
foi         = [0.8 2; 4 6; 7 13; 14 30; 35 48; 50 100];
fit_interv 	= [3 100];
calc_interv = [1 150];
% --------------------------------------------------------

subj = {'bla';'Nola';'Trixie'};
restoredefaultpath

addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/

idir   = '/home/istitt/chronicECoG/';
outdir   = '/home/tpfeffer/critical_ferret/proc/dfa/';
% mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/tpfeffer/critical_ferret/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m
%%

for isubj = 1 : 3
  
  indir = [idir sprintf('%s/recordingsMat/',subj{isubj})];
  d = dir([indir sprintf('*%s_sleep*',lower(subj{isubj}))]);
  
  for irec = 1 : length(d)
    for ifoi = 1 : length(foi)
      
      if ~exist(sprintf([outdir 'cf_ampcorr_dfa_s%d_r%d_f%d_v%d_processing.txt'],isubj,irec,ifoi,v))
        system(['touch ' outdir sprintf('cf_ampcorr_dfa_s%d_r%d_f%d_v%d_processing.txt',isubj,irec,ifoi,v)]);
      else
        continue
      end
      
      load(sprintf([outdir 'cf_ampcorr_s%d_r%d_f%d_v%d.mat'],isubj,irec,ifoi,v));

      siginfo = nbt_Info;
      siginfo.converted_sample_frequency = 1;
        
      for ichan = 1 : 64
        tmp  = nbt_doDFA(squeeze(powcorr(ichan,:,:))', siginfo, fit_interv,calc_interv,0.5,0,0,[]);
        ampcorr_dfa(:,ichan) = tmp.MarkerValues; clear tmp
      end
      
      save(sprintf([outdir '/cf_ampcorr_dfa_s%d_r%d_f%d_v%d.mat'],isubj,irec,ifoi,v),'ampcorr_dfa','-v7.3');
      
    end    
  end
end

error('all done!');

%% PLOT CORRELATIONS

v = 2;

for isubj = 2
  for ifoi = 1 : 6
    
    d = dir(sprintf([outdir 'cf_ampcorr_dfa_s%d_r*_f%d_v%d.mat'],isubj,ifoi,v));
    
    for irec = 1:length(d)
      
      load(sprintf([outdir 'cf_ampcorr_dfa_s%d_r%d_f%d_v%d.mat'],isubj,irec,ifoi,v));
      load(sprintf([outdir 'cf_ampcorr_s%d_r%d_f%d_v%d.mat'],isubj,irec,ifoi,v));
      
      for ichan = 1 : 64
        load(sprintf([outdir 'cf_dfa_states_s%d_r%d_c%d_f%d_v%d.mat'],isubj,irec,ichan,ifoi,v));
        dfa(ifoi,irec,ichan) = tmp.MarkerValues;
      end
%       
      ex(ifoi,irec,:)    = nanmean(ampcorr_dfa,2);
      pcorr(ifoi,irec,:) = nanmean(nanmean(powcorr,3),2);

    end
  end
  
  NFOI = 6;
  
  for irec = 1:length(d)
    
    h=figure; set(gcf,'color','white'); 

    for ifoi = 1 : 6
      
      subplot(3,6,(0)*NFOI+ifoi);   
      
      [xi yi zi] = cf_interpECOG(squeeze(dfa(ifoi,irec,:)),'cubic');
      zi(zi==0)=NaN;
      h=imagesc(zi);
      set(h,'alphadata',~isnan(zi)); box off
      clear xi yi zi
      axis off
      colorbar

      subplot(3,6,(1)*NFOI+ifoi);   
      
      [xi yi zi] = cf_interpECOG(squeeze(ex(ifoi,irec,:)),'cubic');
      zi(zi==0)=NaN;
      h=imagesc(zi);
      set(h,'alphadata',~isnan(zi)); box off
      clear xi yi zi
      axis off
      colorbar
      
      subplot(3,6,(2)*NFOI+ifoi);   
      
      [xi yi zi] = cf_interpECOG(squeeze(pcorr(ifoi,irec,:)),'cubic');
      zi(zi==0)=NaN;
      h=imagesc(zi);
      set(h,'alphadata',~isnan(zi)); box off
      clear xi yi zi
      axis off
      colorbar
       
    end
    
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Position',[50 50 1600 800])

    print(gcf,'-djpeg100',sprintf('~/critical_ferret/proc/plots/cf_overview_topos_s%d_r%d_v%d.jpg',isubj,irec,v))
%     close
  end
end

%% FC MATRIX AND CO-DFA MATRIX
    
     
v = 2;

figure; set(gcf,'color','white');

% irec = 1;
% isubj = 2;

for isubj = 1
  
  h=figure; set(gcf,'color','white'); set(gca, 'TickDir','out')
  
  for ifoi = 1 : 6
    
    d = dir(sprintf([outdir 'cf_ampcorr_dfa_s%d_r*_f%d_v%d.mat'],isubj,ifoi,v));
    
    for irec = 2%:1%length(d)
      
      load(sprintf([outdir 'cf_ampcorr_dfa_s%d_r%d_f%d_v%d.mat'],isubj,irec,ifoi,v));
      load(sprintf([outdir 'cf_ampcorr_s%d_r%d_f%d_v%d.mat'],isubj,irec,ifoi,v));
     
      ampcorr_dfa(1:65:end)=NaN;
      tmp = nanmean(powcorr,3); tmp(1:65:end)=NaN;
        
      ex(ifoi,irec,:,:)    = ampcorr_dfa;
      pcorr(ifoi,irec,:,:) = tmp;

%       end
    end
  end
  
  NFOI = 6;
  
  all_dfa   = squeeze(nanmean(ex,2));
  all_pcorr = squeeze(nanmean(pcorr,2));
  
  for irec = 2%length(d)
    
    for ifoi = 1 : 6
      
      subplot(length(d),6,(irec-1)*NFOI+ifoi);   title(sprintf('freq%d',isubj))

      h=imagesc(squeeze(all_dfa(ifoi,:,:)));
      set(h,'alphadata',~isnan(squeeze(all_dfa(ifoi,:,:)))); box off
      axis off
      colorbar
      
      subplot(length(d),6,(1+irec-1)*NFOI+ifoi);   title(sprintf('freq%d',isubj))

      h=imagesc(squeeze(all_pcorr(ifoi,:,:)));
      set(h,'alphadata',~isnan(squeeze(all_pcorr(ifoi,:,:)))); box off
      axis off
      colorbar
      
      subplot(length(d),6,(2+irec-1)*NFOI+ifoi);   title(sprintf('freq%d',isubj))

      h=imagesc(squeeze(all_pcorr(ifoi,:,:)));
      set(h,'alphadata',~isnan(squeeze(all_pcorr(ifoi,:,:)))); box off
      axis off
      colorbar
      
    end
  end
end

%% FLUC FUNCTION

%
clear tmp

irec = 1;
for subj = 1 : 1
  
  for ichan = 1 : 1
    load(sprintf([outdir 'cf_dfa_states_s%d_r%d_c%d_f%d_v%d.mat'],isubj,irec,ichan,ifoi,v));
    ex(ifoi,ichan) = tmp.MarkerValues;
  end
  %   isubj = SUBJLIST(subj);
  %   cnt = cnt + 1;
  
  %   a = squeeze(l(:,:,isubj,:)); a(a==0)=NaN;
  %   s = nanstd(a,[],3);
  %   m = nanmean(a,3);
  %   ci = [m(:,2)-s(:,2)/sqrt(size(a,3)) m(:,2)+s(:,2)/sqrt(size(a,3))];
  
  calc_small  = min(find(tmp.DFA_x>=tmp.CalcInterval(1)*tmp.Fs));		%
  calc_large  = max(find(tmp.DFA_x<=tmp.CalcInterval(2)*tmp.Fs));
  fit_small   = min(find(tmp.DFA_x>=tmp.FitInterval(1)*tmp.Fs));
  fit_large   = max(find(tmp.DFA_x<=tmp.FitInterval(2)*tmp.Fs));
  
  figure; set(gcf,'color','white'); hold on
  
  plot(log10(tmp.DFA_x(fit_small:fit_large)/tmp.Fs),log10(tmp.DFA_y{1}(fit_small:fit_large)),'ro','MarkerFace','r','MarkerSize',3)
  h=lsline; set(h,'color','b','LineWidth',2);
  
  plot(log10(tmp.DFA_x/tmp.Fs),log10(tmp.DFA_y{1}),'ro','MarkerFace','r','MarkerSize',3)
  
  %   plot(log10(m(:,1)),log10(ci(:,1)),'r--')
  %   plot(log10(m(:,1)),log10(ci(:,2)),'r--')
  
  axis([log10(min(tmp.DFA_x)/tmp.Fs-0.2) log10(max(tmp.DFA_x)/tmp.Fs+20) log10(min(tmp.DFA_y{1}(3:end))-2) log10(max(tmp.DFA_y{1})+3)])
  
  
  line([log10(tmp.FitInterval(1)) log10(tmp.FitInterval(1))],[0 2.5],'color','k','LineStyle','--')
  line([log10(tmp.FitInterval(2)) log10(tmp.FitInterval(2))],[0 2.5],'color','k','LineStyle','--')
  
  set(gca,'XTick',[log10([tmp.DFA_x/tmp.Fs])],'XTickLabel',[round(10*tmp.DFA_x/tmp.Fs)/10])
  %   xlabel('log10(window size)'); ylabel('log10(<F(t)>)');
  set(gca,'TickDir','out'); xticklabel_rotate([],45)
  
  xlabel('Window size (in minutes)')
  
  %   e = nanmean(expo(expo~=0&~isnan(expo)));
  %   e = nanmean(expo(isubj,:));
  %   s = nanstd(expo(expo~=0&~isnan(expo)))/sqrt(size(expo(expo~=0&~isnan(expo)),1));
  %   s = nanstd(expo(isubj,:));
  
  title(sprintf('S%d: dfa = %0.2f',isubj,tmp.MarkerValues))
  
  %   X = [ones(1,fit_large-fit_small+1)' log10(tmp.DFA_x(fit_small:fit_large))'];
  %   Y = log10(tmp.DFA_y{1}(fit_small:fit_large));
  %   DFA_exp = X\Y; %least-square fit
  %   DFA_exp = DFA_exp(2);
  
  clear e s X Y a s m ci
  
  saveas(gcf,sprintf('~/critical_ferret/proc/dfa/cf_dfa_ff_s%d_r%d_v%d.fig',isubj,irec,v),'fig')
  
  
end





