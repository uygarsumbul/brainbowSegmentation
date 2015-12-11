function clusterCount = vdpgm_countEstimate(colorData)
cd ~/bb/vdpgm/
vdpgmopts                  = mkopts_avdp;
vdpgmopts.max_target_ratio = 0.99;
[~, result]                = evalc('vdpgm(transpose(colorData), vdpgmopts)');
clusterCount               = result.K;
cd ~/bb/
