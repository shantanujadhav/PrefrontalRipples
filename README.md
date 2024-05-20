Code accompanying the paper: Shin, J. D. & Jadhav, S. P. (2024). Prefrontal cortical ripples mediate top-down suppression of hippocampal 
reactivation during sleep memory consolidation.  
Current Biology.

Codes were created in MATLAB 2022a. 

FILES and FOLDERS
-----------------
  ./Figure1  
  <b>jds_CA1_PFCrips_xcorr_M.m</b>  -  CA1-PFC ripple cross-correlation. Can use to compute cross-correlation for other events as well.  
  <b>jds_ctxrip_plotExample_M.m</b>  -  Plots example ripple events.  
  <b>jds_plotRaster_sleepProject_M.m</b>  -  Plots example raster with spikes, LFP, and ripples highlighted.  
  <b>jds_rippletrig_spindlepwr_M.m</b>  -  Compares spindle power during coordinated and independent PFC ripple events.  
  <b>jds_sleepripplerate_overtime_M.m</b>  -  Calculates rate of ripples over sleep sessions.  
  <b>jds_triggered_wavelet_M.m</b>  -  Plots the PFC ripple triggered wavelet spectrogram. Uses the tetrode with highest number of ripples.  

  ./Figure2  
  <b>jds_ripMod_ripPropertiesQuartiles_M.m</b>  -  Plots the relationship between PFC ripple features (frequency, amplitude, legnth) and CA1 modulation during independent PFC ripples.  
  <b>jds_ripTrigModLargerWindow_M.m</b>  -  Plots ripple triggered modulation for CA1 and PFC. Run DFSjds_getripalignspiking_M.m to generate modulation files aligned to different ripple events.  

  ./Figure3  
  <b>jds_CA1_rippleBursting_latency_M.m</b>  -  Plots bursting and spike latency metrics for CA1 mod cells during SWRs. Also plots the relationship between CA1 modulation index and bursting/spike latency.  
  <b>jds_CA1mod_bumpPlots_M.m</b>  -  Plots example place fields sorted by peak on the maze.  
  <b>jds_CA1mod_getFieldMetricsPrePost_M.m</b>  -  Compares place field metrics for CA1 mod cells (number fields, field width, etc).  
  <b>jds_CA1mod_spatialinfo_M.m</b>  -  Calculates and plot spatial information for CA1 mod cells.  
  <b>jds_CA1mod_stability_M.m</b>  -  Calculates stability of CA1 mod cells for consecutive run sessions.  
  <b>jds_PPC_thetalockingstrength_M.m</b>  -  Calculates theta phase locking strength for CA1 mod cells (pairwise phase consistency-PPC, and kappa concentration parameter).  
  <b>jds_comparePFCripmod_SWRparticipation_M.m</b>  -  Calculates and plots NREM firing rate, firing rate gain during SWRs, etc. for CA1 mod cells.  
  <b>jds_ripTrigMod_compareRips_scatter_M.m</b>  -  Plots the relationship between modulation during coordinated ripples and independent PFC ripples.  

  ./Figure4  
  <b>jds_GLM_ripplepredictionArea_M.m</b>  -  Predicts ripple type (coordinated or independent) in CA1 or PFC using CA1 or PFC population activity, respsectively.  
  <b>jds_assemblySuppressionReinstatement_M.m</b>  -  Evaluates the relationship between assembly reinstatement and CA1 assembly suppression during independent PFC ripples.  
  <b>jds_noncoordrippletriggered_assemblystrength_M.m</b>  -  Ripple triggered assembly reactivation plots.  
  <b>jds_plotRasterReactivation_sleepProject_M.m</b>  -  Example assembly reactivation rasters for CA1 and PFC.  
  <b>jds_rippletriggered_compareStrengthAll_M.m</b>  -  Calculates and compares the reactivation strength of assemblies during different ripple events.  

  ./Figure5  
  <b>jds_CA1_SWRdimensionality_M.m</b>  -  Calculates and compares the dimensionality of independent and coordinated SWRs.  
  <b>jds_CA1replay_rippleTypeSeqMetrics_M.m</b>  -  Compares weighted correlation, jump distances, R2, % reverse replay, and % significant replay for independent and coordinated SWR replay (line fitting decoding method).  
  <b>jds_PCCmodNonmod_M.m</b>  -  Compares per cell comtribution (PCC) for CA1 mod and nonmod cells.  
  <b>jds_PCCrippleReplayDifference_M.m</b>  -  Compares PCC differences between coordinated and independent SWR replay for CA1 mod and nonmod cells.  
  <b>jds_SWR_rankordercorr_M.m</b>  -  Computes rank order correlations for independent and coordinated SWR events.  
  <b>jds_plotReplayRaster_M.m</b>  -  Plots example replay rasters.  
  <b>jds_sequenceRobustness_withawake_M.m</b>  -  Compares and plots sequence degradataion of independent, coordinated, and awake SWR replay events (single cell shuffles).  
  <b>jds_sequenceScoreIndCoord_M.m</b>  -  Compares and plots weighted correlation and jump distance for independent and coordinated events (rZ shuffle method).  

  ./Figure6  
  <b>jds_SOphaselocking_M.m</b>  -  Calculates and plots slow oscillation phase locking of spindles, pfc ripples, and CA1 SWRs.  
  <b>jds_periSO_eventProbability_M.m</b>  -  Plots the event occurences of spindles, pfc ripples, and CA1 SWRS surrounding SO troughs (up states).  
  <b>jds_periSO_eventProbability_ReactivationQuartiles_M.m</b>  -  Plots peri-SO CA1 reacvtivation strengths according to quartiles.  
  <b>jds_periSO_eventProbability_ReactivationSuppression_M.m</b>  -  Plots the relationship between peri-SO CA1 reactivation strength and independent PFC ripple triggered CA1 assembly suppression.  
  <b>jds_periSO_eventProbability_Reactivation_M.m</b>  -  Plots SO trough aligned reactivation strength for CA1 and PFC to determine reactivation timing surrounding events.  
  <b>jds_periSO_eventProbability_bootstrapTiming_M.m</b>  -  Resampling of occurence matrix to determine significance of peaks of fold change probability plots.  
  <b>jds_updown_replayProbabilityCoordNoncoord_M.m</b>  -  Calculates the probability of SWR replay events during up and down states of the slow oscillations.  

  ./SupplementalFigures  
  <b>jds_CA1PFCcoreactivationRates_M.m</b>  -  Calculates CA1-PFC co-reactivation rates during different ripple types.  
  <b>jds_CA1relativeSuppression_M.m</b>  -  Calculates suppression of CA1 assemblies during independent PFC ripples at different points in time (early, middle, late sleep).  
  <b>jds_GLM_ripplepredictionAreaLeadLag_M.m</b>  -  Prediction of leading or lagging ripple type based on CA1-PFC activity.  
  <b>jds_SWS_RepFig_M.m</b>  -  Plots representative sleep figure.  
  <b>jds_assemblyContribSpatialCorr_M.m</b>  -  Calculates the spatial correlation of cell pairs for assembly member and nonmember cells.  
  <b>jds_cellsactive_rip_M.m</b>  -  Calculates the proportion of cells active during different ripple events.  
  <b>jds_compareRipCofiring_coordNoncoord_members_M.m</b>  -  Calculates the z-scored ripple co-firing of member cell pairs during different ripple types.  
  <b>jds_compare_pcweights_contribution_M.m</b>  -  Evaluates the relationship between absolute assembly weight and assembly contribution.  
  <b>jds_leadlagrippletriggered_assemblystrength_M.m</b>  -  Plots the leading/lagging ripple triggered assembly strength for CA1 or PFC.  
  <b>jds_plotSpikeInfo_M.m</b>  -  Plots clustering information, CA1 PYR and INT separation, and spiking metrics for CA1 neurons.  
  <b>jds_relativeRippleRate_M.m</b>  -  Compares ripple rates for first and second half of sleep epochs.  
  <b>jds_ripTrigMod_modProportions.m</b>  -  Calculates and plots the proportion of EXC, INH, and modulated cells in CA1 and PFC aligned to different ripple types.  
  <b>jds_ripplecoordinationtriggered_assemblystrength_M.m</b>  -  Plots coordinated ripple triggered assembly strength for CA1 and PFC.  
  <b>jds_rippletrig_crossripplepwr_M.m</b>  -  Plots the ripple power in the opposing area aligned to either independent or coordinated ripple events.  
  <b>jds_rippletriggered_CA1PFCassemblies_M.m</b>  -  Plots ripple triggerend CA1-PFC joint reactivation strength and looks at the relationship between different ripple events.  
  <b>jds_sleepEventRates_comparerunsleep_M.m</b>  -  Plots the rates of sleep oscillations (independent and coordinated ripples, SOs, spindles, run ripples).  
  <b>jds_sleepripplelengths_overtime_M.m</b>  -  Calculates ripple lengths over epochs.  
  <b>jds_sleepripplerate_performance_M.m</b>  -  Evaluates the relationship between ripple rate and behavioral performance.  
  <b>jds_theta_phaselocking_modcells_M.m</b>  -  Plots the preferred theta phase of CA1 mod cells during run.  

  ./ProcessData  
  <b>DFSjds_getripalignspiking_M.m</b>  -  Get ripple aligned modulation.  
  <b>jds_CA1PFCassemblyReactivation_ICA_M.m</b>  -  Get joint CA1-PFC assemblies.  
  <b>jds_PCC_replay_decoding_CA1_allcells_M.m</b>  -  Replay decoding using rZ shuffle and per cell contribution calculations.  
  <b>jds_assemblyReactivationReinstatement_ICA_M.m</b>  -  Get CA1 or PFC assemblies for reinstatement analysis.  
  <b>jds_extractSlowOscillations_M.m</b>  -  Extract slow oscillations.  
  <b>jds_getSpikeInfo_M.m</b>  -  Get cell cluster information for subsequent plotting.  
  <b>jds_getripples_amp_freq_sleep_M.m</b>  -  Get frequency and amplitude or ripple events.  
  <b>jds_reactivation_contribution_M.m</b>  -  Get contribution to activation/reactivation strength for comparison with assembly weights.  
  <b>jds_replay_decoding_CA1_M.m</b>  -  Replay decoding using line fitting and monte carlo method.  
  <b>jds_replay_spikematrix_ev_singleday_M.m</b>  -  Get spike matrices during run and sleep using specified time bin sizes.  

CONTACT
=======
Justin D. Shin
jdshin@brandeis.edu
