clear all

% --------------------------------------------------------
% VERSION 1
% --------------------------------------------------------
v         = 1;
% --------------------------------------------------------

subj = {'bla';'Nola';'Trixie'};
restoredefaultpath

outdir   = '/home/tpfeffer/critical_ferret/proc/aval/';
plotdir = '/home/tpfeffer/critical_ferret/proc/plots/';
addpath ~/pconn/matlab
%%

for isubj = 1 : 3
  
  d = dir([outdir sprintf('cf_aval_evts_s%d_r*_v%d.mat',isubj,v)]);
  
  for irec = 1 : length(d)
    
    if ~exist(sprintf([outdir 'cf_aval_ana_s%d_r%d_v%d_processing.txt'],isubj,irec,v))
      system(['touch ' outdir sprintf('cf_aval_ana_s%d_r%d_v%d_processing.txt',isubj,irec,v)]);
    else
      continue
    end
    
    load([outdir d(irec).name]);
    
    fprintf('Analyzing avalanches s%d r%d ...\n',isubj,irec)
    
    for ti = 1 : size(all_evts,1)
      
      t(ti) = any(all_evts(ti,:));
      
    end
    
    cnt = 0;
    
    while 1
      
      idx=find(t,1,'first');
      
      cnt = cnt + 1;
      
      if isempty(idx)
        break
      end
      
      step = 1; zerostep = 0; cntdown = 27; take_step  = 1;
      
      if ~[(idx+step)>=size(all_evts,1)]
        while cntdown~=0
          take_step = step + zerostep;
          if ~[(idx+take_step)>=size(all_evts,1)]
            if any(all_evts(idx+take_step,:)==1)
              step = step + zerostep;
              zerostep = 0;
              step = step + 1;
            else
              zerostep = zerostep + 1;
              cntdown = cntdown - 1;
            end
          else
            break
          end
        end
      end
      
      evt_dur(cnt) = step;
      
      evt_idx = idx:idx+step-1;
      
      evt_size(cnt)=sum(sum(all_evts(evt_idx,:)));
      
      t(evt_idx) = 0;
      
    end
    
    fprintf('Analyzing avalanches s%d r%d ... %d avalanches found! \n',isubj,irec,cnt)
    
    for i = 1 : max(evt_size)
      h(i) = length(find(evt_size==i));
    end
    
    tmp = pconn_regress(log10(unique(evt_size)),log10(h(h~=0))');
    slp = tmp(2);
    
    % collect data
    aval.slp = tmp; clear tmp slp
    aval.evt_size = evt_size;
    aval.evt_dur = evt_dur;
    aval.aval_hist = h;
    
    
    save(sprintf([outdir 'cf_aval_ana_s%d_r%d_v%d.mat'],isubj,irec,v),'aval','-v7.3');
  end
end

error('all done')

%% PLOT STUFF

figure; set(gcf,'color','white');
plot(log10(unique(aval.evt_size)),log10(aval.aval_hist(aval.aval_hist~=0)),'.');
xlabel('Avalanche size'); ylabel('Frequency');
set(gca,'tickdir','out','xtick',log10([2 4 8 16 32 64 128 256 512]),'xticklabel',[2 4 8 16 32 64 128 256 512])

%% relate to DFA

isubj = 1;
irec = 1;
ifoi = 3;
v= 2;

for irec = 1 : 8

  for ichan = 1 : 64
    load(sprintf(['/home/tpfeffer/critical_ferret/proc/dfa/' 'cf_dfa_states_s%d_r%d_c%d_f%d_v%d.mat'],isubj,irec,ichan,ifoi,v));
    dfa(ichan) = tmp.MarkerValues;
  end

  load(sprintf([outdir 'cf_aval_ana_s%d_r%d_v%d.mat'],isubj,irec,1));

  dfaaval(irec,:) = [mean(dfa) aval.slp(2)];
  
end

