%% COMPUTE WHOLE BRAIN DFA IN SOURCE SPACE
% cf_dfa_states

clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
% --------------------------------------------------------

subj = {'bla';'Nola';'Trixie'};
restoredefaultpath

addpath /home/tpfeffer/Documents/MATLAB/toolboxes/NBT-ReleaseNBTRC4a/

idir   = '/home/istitt/chronicECoG/';
outdir   = '/home/tpfeffer/critical_ferret/proc/aval/';
% mridir = '/home/gnolte/neuconn/mri_data/';
plotdir = '/home/tpfeffer/critical_ferret/proc/plots/';

run ~/Documents/MATLAB/toolboxes/NBT-NBTv0.5.3-alpha/installNBT.m

hpfreq = 5;
lpfreq = 100;

%%

for isubj = 1 : 3
  
  indir = [idir sprintf('%s/recordingsMat/',subj{isubj})];
  d = dir([indir sprintf('*%s_sleep*',lower(subj{isubj}))]);
  
  for irec = 1 : length(d)
    
    if ~exist(sprintf([outdir 'cf_aval_s%d_r%d_v%d_processing.txt'],isubj,irec,v))
      system(['touch ' outdir sprintf('cf_aval_s%d_r%d_v%d_processing.txt',isubj,irec,v)]);
    else
      continue
    end
    
    for ichan = 1 : 64
      
      n = d(irec).name;
      
      indir1 = [indir n '/LFP/'];
      
      fprintf('Processing %s rec%d channel %d ...\n',subj{isubj},irec,ichan);
      
      load([indir1 sprintf('LFP_1_EU64L_%d.mat',ichan)])
      
      lfp = nbt_filter_fir(LFP,hpfreq,lpfreq,lfpFs,.8);
      thresh       = std(lfp);
      
      evts = abs(lfp)>3*thresh;
      
      all_evts(:,ichan) = zeros(size(lfp,1),1);
      
      fprintf('Isolating events ...');

      while 1
        
%         fprintf('.')
        
        idx=find(evts,1,'first');
        
        if isempty(idx)
          break
        end
        
        step = 1;
        if ~idx+step>=size(evts,1)                 
          while evts(idx+step)==1 
            step = step+1;
          end
        end
        evt_idx = idx:idx+step-1;
        
        [~,peak]=max(abs(lfp(evt_idx)));
        
        evts(evt_idx) = 0;
        all_evts(evt_idx(peak),ichan) = 1;
  
      end
      
      fprintf('Done!\n')
        
    end
    
    save(sprintf([outdir 'cf_aval_evts_s%d_r%d_v%d.mat'],isubj,irec,v),'all_evts','-v7.3');
    
    clear all_evts

  end
end
    
error('All done!')


