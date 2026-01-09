#pragma once

#include "aureonoise3beta/core/state.hpp"

// Reset/report utilities
void aureonoise_reset_reports(t_aureonoise* x);

// Lifecycle helpers
void aureonoise_prepare_instance(t_aureonoise* x);
void aureonoise_prepare_dsp(t_aureonoise* x, double sr);
void aureonoise_update_sync(t_aureonoise* x, double block_dt);
void aureonoise_set_tempo(t_aureonoise* x, double bpm);
void aureonoise_update_modal(t_aureonoise* x);

// Spatial profile helpers
void aureonoise_spatial_profile_refresh(t_aureonoise* x);
void aureonoise_spatial_profile_params_changed(t_aureonoise* x);

// Max/MSP hooks
void aureonoise_clear(t_aureonoise* x);
void aureonoise_dsp64(t_aureonoise* x, t_object* dsp64, short* count, double sr, long n, long flags);
void aureonoise_perform64(t_aureonoise* x, t_object* dsp64, double** ins, long numins,
                          double** outs, long numouts, long sampleframes, long flags, void* userparam);
