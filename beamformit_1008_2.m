% no 'overall weighting'
% no 'post processing'
% no 'channels elimination'

function beamformit_1008(infilenames_cell, outfilename)

nmic = size(infilenames_cell,2);
npair = nmic - 1;

x_orig = cell(1,nmic);
sr = 16000;

for m = 1:nmic
    [x_orig{m}, sr] = audioread(infilenames_cell{m}); % estimated ref_mic:
end
nsample = size(x_orig{1},1);
x = zeros(nmic, nsample);

for m = 1:nmic
    x(m,:) = x_orig{m}';
end

%% make hamming window
%%


nwin = 16000; % 1 sec
%     hamm_val = 0.54 - 0.46*cos(6.283185307*i/(window-1));
win = zeros(1,nwin);
for i = 1:nwin
    win(i) = 0.54 - 0.46 * cos(6.283185307*(i-1)/(nwin-1));
end

% figure; plot(win);
% disp(win(1:10)');

%% 
%% testing calculating xcorr
%%

% Comparing from frame 300 to 16300
% 1st xcorr calculation iteration 1: 0.239386
% 2nd xcorr calculation iteration 1: 0.152243
% avg xcorr between 1 and 2
% 0.230185

npiece = 200;
nfft = 16384;
nbest = 2;

scroll = floor(nsample / (npiece+2)); 

stft1 = fft([x(1,scroll:(scroll+nwin-1)) .* win, zeros(1,nfft-nwin)]); 
stft2 = fft([x(2,scroll:(scroll+nwin-1)) .* win, zeros(1,nfft-nwin)]); 

numerator = stft1 .* conj(stft2);
ccorr = real(ifft(numerator ./ (abs(numerator))));
ccorr = [ccorr(1:481), ccorr(end-480+1:end)];
best2_ccorr = maxk(ccorr, 2);

% disp('scroll');
% disp(scroll);
% plot(ccorr);
% disp(best2_ccorr); 

%% calculate avg_ccorr
%%
avg_ccorr = zeros(nmic, nmic);

for i = 1:npiece
    st = i * scroll;
    ed = st + nwin - 1;
    if st + nfft >= nsample
        break;
    end

    for m1 = 1:(nmic-1)
        avg_ccorr(m1, m1) = 0;
        for m2 = (m1+1):nmic
            stft1 = fft([x(m1,st:ed) .* win, zeros(1,nfft-nwin)]); 
            stft2 = fft([x(m2,st:ed) .* win, zeros(1,nfft-nwin)]);
            numerator = stft1 .* conj(stft2);
            ccorr = real(ifft(numerator ./ (abs(numerator))));
            ccorr = [ccorr(1:481), ccorr(end-480+1:end)];
           
            avg_ccorr(m1, m2) = avg_ccorr(m1, m2) + sum(maxk(ccorr, nbest));
            avg_ccorr(m2, m1) = avg_ccorr(m1, m2); 
        end
    end
end

avg_ccorr = avg_ccorr / (nbest * npiece * (nmic-1));
% disp(avg_ccorr);
% disp(sum(avg_ccorr,1));

%% calculating scaling factor
%%
% Set the total number of segments from 10 to 1 for computing scalling factor, with segment duration 60798
% Processing channel 0
% amount_frames_read: 0.000000
% segment_duration: 0.000000
% The Median maximum energy for channel 0 is 0.352051
% Set the total number of segments from 10 to 1 for computing scalling factor, with segment duration 60798
% Processing channel 1
% amount_frames_read: 0.000000
% segment_duration: 0.000000
% The Median maximum energy for channel 1 is 0.341858
% Set the total number of segments from 10 to 1 for computing scalling factor, with segment duration 60798
% Processing channel 2
% amount_frames_read: 0.000000
% segment_duration: 0.000000
% The Median maximum energy for channel 2 is 0.499451
% Weighting calculated to adjust the signal: 0.754173

nsegment = 10;

max_val = zeros(nmic, 1);

if size(x,2) <= 160000 % 10 seconds
    for m = 1:nmic
        max_val(m) = max(abs(x(m,:)));
    end
else
    if size(x,2) < 1600000 % 100 seconds
        nsegment = ceil(size(x,2) / 16000);
    end
    scroll = floor(size(x,2) / nsegment);
    max_val_candidate = zeros(nmic, nsegment);
    
    for s = 0:(nsegment-1)
        st = s * scroll + 1;
        ed = st + scroll - 1;
        for m = 1:nmic
            max_val_candidate(m,s+1) = abs(max(x(m,st:ed)));
        end
    end
    
    for m = 1:nmic
        
        sorted = sort(max_val_candidate(m,:), 'ascend');
        max_val(m) = sorted(floor(end/2) + 1);
    end
end

overall_weight = (0.3 * nmic) / sum(max_val);

% disp(max_val);
% disp(overall_weight);
% disp((0.3 * 3) / (0.352051 + 0.341858 + 0.499451));

%% setting channel weight
%%


%% compute total number of delays
%%
% int totalNumDelays
% % = (int)((m_frames - (*m_config).windowFrames - m_biggestSkew - m_UEMGap)/((*m_config).
% rate*m_sampleRateInMs));
% sr_in_ms = 16000 / 1000; % 16

% too complicated. I should do hard coding
nframe = floor(( nsample - 8000 ) / (4000));

% disp(nframe);

%% recreating hamming window
%%
nwin = 8000; % 0.5 sec
%     hamm_val = 0.54 - 0.46*cos(6.283185307*i/(window-1));
win = zeros(1,nwin);
for i = 1:nwin
    win(i) = 0.54 - 0.46 * cos(6.283185307*(i-1)/(nwin-1));
end

% figure; plot(win);
% disp(win(1:10));
% disp(16000 * 30 / 1000);

%% compute TDOA

nbest = 4;
[dummy, ref_mic] = max(sum(avg_ccorr));
%ref_mic = 2;
micpair = zeros(nmic, npair);


gcc_nbest = zeros(npair, nframe, nbest);
tdoa_nbest = zeros(npair, nframe, nbest);
lrfilp_gcc = zeros(npair, nframe, nfft);

for m = 1:nmic
    p = 1; % pair idx
    for i = 1:nmic
        if i == m
            continue;
        end
        micpair(m,p) = i;
        p = p + 1;
    end
end


for t = 1:(nframe)
    st = (t-1) * 4000 + 1;
    ed = st + nwin - 1;
%     disp(st);
%     disp(ed);
        for p = 1:npair
            
            m = micpair(ref_mic,p);
            
            stft_ref = fft([x(ref_mic,st:ed) .* win, zeros(1,nfft-nwin)]); 
            stft_m = fft([x(m,st:ed) .* win, zeros(1,nfft-nwin)]);
            numerator = stft_m .* conj(stft_ref);
            gcc = real(ifft(numerator ./ (eps+abs(numerator))));
            gcc = [gcc(end-479:end), gcc(1:480)];
            [gcc_nbest(p,t,:), tdoa_nbest(p,t,:)] = maxk(gcc, nbest);
            tdoa_nbest(p,t,:) = tdoa_nbest(p,t,:) - (481); % index shifting
          
        end

end

% disp(squeeze(gcc_nbest(:,:,1)));
% disp(squeeze(tdoa_nbest(:,:,1)));
% disp(squeeze(tdoa_nbest(:,:,2)));

%% find noise threshold
%%
% Threshold is 0.140125 (min: 0.132648 max 0.563892)
% Thresholding noisy frames lower than 0.140125
th_idx = floor((0.1 * nframe)) + 1;

sorted = sort(sum(gcc_nbest(:,:,1),1), 'ascend');

threshold = sorted(th_idx)/npair;
%disp(sorted);
% disp(th_idx);
% disp(threshold);

%% noise filtering
%%
noise_filter = zeros(npair, nframe);

for p = 1:npair
    for t = 1:nframe
        if gcc_nbest(p,t,1) < threshold
            noise_filter(p,t) = 1;
        
            if t == 1  % it's silence
                gcc_nbest(p,t,:) = 0; % masking with discouraging values
                gcc_nbest(p,t,1) = 1; % only one path
                tdoa_nbest(p,t,:) = 480; % masking with discouraging values
                tdoa_nbest(p,t,1) = 0;
            else
                tdoa_nbest(p,t,:) = tdoa_nbest(p,t-1,:);
            end
        end
    end
end

% disp(gcc_nbest(:,:,1));
% disp(tdoa_nbest(:,:,1));
% 
% disp(gcc_nbest(:,:,2));
% disp(tdoa_nbest(:,:,2));


%% single channel viterbi - emission, trans prob
%%
emission1 = zeros(npair, nframe, nbest);
diff1 = zeros(npair, nframe, nbest, nbest); % do not using 1st idx
transition1 = zeros(npair, nframe, nbest, nbest); % do not using 1st idx

for p = 1:npair
    for t = 1:nframe
        for n = 1:nbest 
            if gcc_nbest(p,t,n) == 0
                emission1(p,t,n) = -1000;
            else
                emission1(p,t,n) = log10(gcc_nbest(p,t,n));
            end
        end
    end
end



for p = 1:npair
    for t = 2:nframe
        for n = 1:nbest
            for nprev = 1:nbest
                diff1(p,t,n,nprev) = abs(tdoa_nbest(p,t,n) - tdoa_nbest(p,t-1,nprev));
            end
        end
    end
end

maxdiff1 = max(diff1(:));

% for p = 1:npair
%     for t = 2:nframe
%         for n = 1:nbest
%             if maxdiff1 < max(diff1(p,t,n,:));
%             maxdiff1(p,t,n) = max(diff1(p,t,n,:));
%         end
%     end
% end

for p = 1:npair
    for t = 2:nframe
        for n = 1:nbest
            for nprev = 1:nbest
                % there is a computational bug.
%                 disp((2+maxdiff1(p,t,n)));
%                 disp(diff1(p,t,n,nprev));
%                 disp(maxdiff1(p,t,n));
%                 disp(log10(481/482));
%                 disp(1 + maxdiff1(p,t,n) - diff1(p,t,n,nprev));
%                 disp((2+maxdiff1(p,t,n)));
                nume = 1 + maxdiff1 - diff1(p,t,n,nprev);
                deno = (2+maxdiff1);
                transition1(p,t,n,nprev) = log10(nume / deno); 
%                 disp(transition1(p,t,n,nprev));
            end
        end
    end
end

%% single channel viterbi - searching
%%

nbest2 = 2;
score1 = zeros(npair, nframe, nbest, nbest); % tmp variable
score1_table = zeros(npair, nframe, nbest);
back1_table = zeros(npair, nframe, nbest);
bestpath1 = zeros(npair, nframe, nbest2); % state idx stored.

for p = 1:npair
    
    for n = 1:nbest
        score1(p,1,n,:) = emission1(p,1,n); % broadcasting
        score1_table(p,1,n) = emission1(p,1,n);
    end
end

for p = 1:npair
    for t = 2:nframe
        for n = 1:nbest
            for nprev = 1:nbest
                score1(p,t,n,nprev) = ...
                    score1_table(p,t-1,nprev,1) + ...
                    25*transition1(p,t,n,nprev) + ... 
                    emission1(p,t,n);    
            end
            %score1_table(p,t,n) = max(score1(p,t,n,:));
            [score1_table(p,t,n),back1_table(p,t,n)] = max(score1(p,t,n,:));
        end
    end
end

for p = 1:npair
    [dummy, bestpath1(p,end,1)] = max(score1_table(p,end,:));
end

for p = 1:npair
    for back_t = 0:(nframe-2)        
        t = nframe - back_t;
        bestpath1(p,t-1,1) = back1_table(p,t,bestpath1(p,t,1));
    end
end

% If you put more weight on transition, the back table changes.
% If I had put weight 10 on transition, the back table became intuitive.
% disp(back1_table(1,:,:)); % compare first mic pair trellis plot
% disp(score1_table(1,:,:));
% disp(bestpath1(1,:,1));
%% single channel viterbi - 2-best path
%%

emission1_bak = emission1;
for p = 1:npair
    for t = 1:nframe
        best1 = bestpath1(p,t,1);
        emission1(p,t,best1) = -1000;
    end
end

for p = 1:npair
    for n = 1:nbest
        score1(p,1,n,:) = emission1(p,1,n); % broadcasting
        score1_table(p,1,n) = emission1(p,1,n);
    end
end

for p = 1:npair
    for t = 2:nframe
        for n = 1:nbest
            for nprev = 1:nbest
                score1(p,t,n,nprev) = ...
                    score1_table(p,t-1,nprev,1) + ...
                    25*transition1(p,t,n,nprev) + ... 
                    emission1(p,t,n);    
            end
            [score1_table(p,t,n),back1_table(p,t,n)] = max(score1(p,t,n,:));
        end
    end
end

for p = 1:npair
    [dummy, bestpath1(p,end,2)] = max(score1_table(p,end,:));
end

for p = 1:npair
    for back_t = 0:(nframe-2)        
        t = nframe - back_t;
        bestpath1(p,t-1,2) = back1_table(p,t,bestpath1(p,t,2));
    end
end


% disp(squeeze(bestpath1(1,:,:))');
% disp(squeeze(bestpath1(2,:,:))');


%% multi channel viterbi - state define
%%
nbest2 = 2;
nstate = nbest2 ^ npair;
g = zeros(nstate, npair);
tmp_row = zeros(3,1);
l = 1;
for ibest = 1:nbest2
    [g, l] = fill_all_comb(1, npair, ibest, nbest2, tmp_row, g, l);
end

% disp(g);

%% multi channel viterbi - emission, trans prob

emission2 = zeros(nframe, nstate);
diff2 = zeros(nframe, nstate, nstate); % do not using 1st idx
transition2 = zeros(nframe, nstate, nstate); % do not using 1st idx

for t = 1:nframe
    for l = 1:nstate
        for m = 1:npair
            ibest = bestpath1(m, t, g(l,m));
            if gcc_nbest(m,t,ibest) > 0
                emission2(t, l) = emission2(t, l) + log10(gcc_nbest(m,t,ibest));
	    else
		emission2(t, l) = -1000;
            end
            
        end
    end
end

for t = 2:nframe
    for l = 1:nstate
        for lprev = 1:nstate
            for m = 1:npair
                ibest = bestpath1(m, t, g(l,m));
                jbest = bestpath1(m, t, g(lprev,m));
                transition2(t,l,lprev)...
                    = transition2(t,l,lprev)...
                    + transition1(m,t,ibest,jbest);
            end
        end
    end
end
%% multi channel viterbi - searching

score2 = zeros(nframe, nstate, nstate);
score2_table = zeros(nframe, nstate);
back2_table = zeros(nframe, nstate);

for l = 1:nstate
    score2_table(1,l) = emission2(1,l);
end

for t = 2:nframe
    for l = 1:nstate
        for lprev = 1:nstate
             score2(t,l,lprev) = ...
                score2_table(t-1,lprev) + ...
                25*transition2(t,l,lprev) + ... 
                1*emission2(t,l);
            
        end
        [score2_table(t,l), back2_table(t,l)] = max(score2(t,l,:));
    end
end

bestpath2 = zeros(nframe,1); % state idx stored.

[dummy, bestpath2(end)] = max(score2_table(end,:));

for back_t = 0:(nframe-2)        
    t = nframe - back_t;
    bestpath2(t-1) = back2_table(t,bestpath2(t));
end


besttdoa = zeros(npair, nframe);
bestgcc2 = zeros(npair, nframe);

for t = 1:nframe
    for p = 1:npair
        l = bestpath2(t);
        ibest = bestpath1(p, t, g(l,p));
        besttdoa(p,t) = tdoa_nbest(p,t,ibest);
        bestgcc2(p,t) = gcc_nbest(p,t,ibest);
    end
end

% disp(back2_table(:,:)'); % compare first mic pair trellis plot
% disp(bestpath2');
% disp(besttdoa);
% disp(bestgcc2);


%% compute local xcorr
%%
tmp_localxcorr = zeros(nmic, nmic, nframe);
localxcorr = zeros(nmic, nframe);

mic2refpair = zeros(nmic,1);
mic2refpair(ref_mic) = 0;

for p = 1:npair
    m = micpair(ref_mic,p);
    mic2refpair(m) = p;
end

for t = 1:nframe
    ref_st = (t-1) * 4000 + 1;
    ref_ed = min(ref_st + 8000 - 1, nsample);
    
    for m1 = 1:(nmic-1)
        for m2 = (m1+1):nmic
            
            if m1 == ref_mic
                st1 = ref_st;
                ed1 = ref_ed;
            else
                p = mic2refpair(m1);
                st1 = max(1,ref_st + besttdoa(p,t));
                ed1 = min(nsample, ref_ed + besttdoa(p,t));
            end
            
            if m2 == ref_mic
                st2 = ref_st;
                ed2 = ref_ed;
            else
                p = mic2refpair(m2);
                st2 = max(1,ref_st + besttdoa(p,t));
                ed2 = min(nsample, ref_ed + besttdoa(p,t));
            end
            
            buf1 = x(m1,st1:ed1);
            buf2 = x(m2,st2:ed2);
            
            ener1 = sum(buf1(:).^2);
            ener2 = sum(buf2(:).^2);
            
            min_ed = min(ed1-st1, ed2-st2) + 1;
            tmp_localxcorr(m1,m2,t)...
                = sum(...
                buf1(1:min_ed) .* buf2(1:min_ed)...
                / (ener1 * ener2));
        
            tmp_localxcorr(m2,m1,t) = tmp_localxcorr(m1,m2,t);
        end
    end
end

localxcorr = squeeze(sum(tmp_localxcorr,1));
% disp(localxcorr);

%% compute sum weight
%%
out_weight = ones(nmic, nframe) * ( 1 / nmic);
alpha = 0.05;


for t = 1:nframe
    
    if sum(localxcorr(:,t)) == 0
        localxcorr(:,t) = 1 / nmic;
    end
    
    localxcorr(:,t) = localxcorr(:,t) / sum(localxcorr(:,t));
    
    for m = 1:nmic
        if m == ref_mic
             out_weight(m,t) = ...
            (1-alpha) * out_weight(m,max(1,t-1)) ...
            + alpha * localxcorr(m,t);    
        
        else
            p = mic2refpair(m);
            if noise_filter(p,t) == 0
                out_weight(m,t) = ...
                    (1-alpha) * out_weight(m,max(1,t-1)) ...
                    + alpha * localxcorr(m,t);    
            end
        end
    end
    
    out_weight(:,t) = out_weight(:,t) / sum(out_weight(:,t));
    
end

% disp(mic2refpair);
% disp(out_weight);
%% 
% 
% 
% 
%% Channel sum
%%
out_x = zeros(1,nsample);

% figure;
for t = 1:nframe
    ref_st = (t-1) * 4000 + 1;
    ref_ed = min(ref_st + 8000 - 1, nsample);
    
    for m = 1:nmic   
        if m == ref_mic
            st = ref_st;
            ed = ref_ed;
        else
            p = mic2refpair(m);
            st = max(1,ref_st + besttdoa(p,t));
            ed = min(nsample, ref_ed + besttdoa(p,t));
        end
       
        triwin = triang(8000)';
%         if t == 1
%             triwin(1:4000) = 1;
%         end
%         size(out_x(ref_st:min_ed))
%         size(triwin(1:min(8000,min_ed-ref_st+1)))
        
         diff = 0;
%         if ref_ed > ed
%             diff = ref_ed - ed;
%         end

        if (ref_ed - ref_st) ~= (ed - st) % if buf is small (always)
          
                diff = ref_ed - ed;
          
        end
        out_x(ref_st:ref_ed-diff)...
            = out_x(ref_st:ref_ed-diff)...   
            + (squeeze(x(m,st:ed))...
            * out_weight(m,t)...
            .* triwin(1:min(8000,ed-st+1))...
            * overall_weight);
            
        
            %.* triwin(min(8000,min_ed-ref_st+1))...
            %* out_weight(m,t));
%         plot(out_x);
    end
end

ref_st = (t) * 4000 + 1;
ref_ed = min(ref_st + 8000 - 1, nsample);

for m = 1:nmic
    if m == ref_mic
        st = ref_st;
        ed = ref_ed;
    else
        p = mic2refpair(m);
        st = max(1,ref_st + besttdoa(p,t));
        ed = min(nsample, ref_ed + besttdoa(p,t));
    end
    buf = squeeze(x(m,st:ed));
    diff = (ref_ed - ref_st) - (ed-st);
    if diff > 0
        buf = [buf, zeros(1,diff)];
    else
        buf = buf(1:end-diff);
    end
    
    triwin = triang(8000)';
%     triwin(4001:end) = 1;
    %buf = buf .* triwin(1:size(buf,2));
    fprintf('diff: %d, ref_ed: %d, ref_st %d\n',diff, ref_ed, ref_st);
    fprintf('ed: %d, st %d\n', ed, st);
    fprintf('buf size: %d\n',size(buf,2));

    out_x(ref_st:ref_ed)...
            = out_x(ref_st:ref_ed)...   
            + (buf...
            * out_weight(m,t)...
            .* triwin(1:min(8000,ref_ed-ref_st+1))...
            * overall_weight);
%     if ref_ed < nsample
%         out_x(ref_ed+1:end) = out_x(ref_ed+1:end) + (x(ref_mic,ref_ed+1:end) * out_weight(m,t) * overall_weight);
%     end
end

%
% figure;
% subplot(4,1,1); plot(x(1,16000:16020)); % (+) delayed
% subplot(4,1,2); plot(x(2,16000:16020)); % ref_mic assumed
% subplot(4,1,3); plot(x(3,16000:16020)); % no idea
% subplot(4,1,4); plot(out_x(16000:16020));
% figure; plot(out_x);
audiowrite(outfilename,out_x,16000);
%%

end
function [table, l] = fill_all_comb(ipair, npair, ibest, nbest, tmp_row, table, l)
    tmp_row(ipair) = ibest;
    %fprintf('ipair: %d ibest: %d l: %d\n', ipair, ibest, l);
    
    if ipair == npair    
        for j = 1:npair
            table(l, j) = tmp_row(j);
        end
        l = l + 1;
    else
        for ibest = 1:nbest
            [table, l] = fill_all_comb(ipair + 1, npair, ibest, nbest, tmp_row, table, l);
        end
    end
end

function [sorted, indices] = maxk(list, k)
    [sorted, indices] = sort(list,'descend');
    sorted = sorted(1:k);
    indices = indices(1:k);
end
