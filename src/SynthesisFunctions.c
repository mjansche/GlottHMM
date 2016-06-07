/*
 * The MIT License is a permissive free software license, which permits reuse
 * within both open source and proprietary software. The software is licensed
 * as is, and no warranty is given as to fitness for purpose or absence of
 * infringement of third parties' rights, such as patents. Generally use for
 * research activities is allowed regardless of any third party patents, but
 * commercial use may be subject to a separate license.
 *
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Tuomo Raitio
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 *
 *
 *
 * <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 *                  GlottHMM Speech Synthesis
 * <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 *
 * This program reads speech parameters and a glottal pulse/
 * pulse library, and synthesizes speech from them.
 *
 * This program has been written in Aalto University,
 * Department of Signal Processign and Acoustics, Espoo, Finland
 *
 * Main author: Tuomo Raitio
 * Acknowledgements: Antti Suni, Paavo Alku, Martti Vainio
 *
 * File SynthesisFunctions.c
 * Version: 1.1
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sndfile.h> 				/* Read and write wav */
#include <gsl/gsl_vector.h>			/* GSL, Vector */
#include <gsl/gsl_matrix.h>			/* GSL, Matrix */
#include <gsl/gsl_fft_real.h>		/* GSL, FFT */
#include <gsl/gsl_fft_complex.h>	/* GSL, FFT complex */
#include <gsl/gsl_fft_halfcomplex.h>/* GSL, FFT halfcomplex */
#include <gsl/gsl_permutation.h>	/* GSL, Permutations */
#include <gsl/gsl_linalg.h>			/* GSL, Linear algebra */
#include <gsl/gsl_spline.h>			/* GSL, Interpolation */
#include <gsl/gsl_errno.h>			/* GSL, Error handling */
#include <gsl/gsl_poly.h>			/* GSL, Polynomials */
#include <gsl/gsl_sort_double.h>	/* GSL, Sort double */
#include <gsl/gsl_sort_vector.h>	/* GSL, Sort vector */
#include <gsl/gsl_complex.h>		/* GSL, Complex numbers */
#include <gsl/gsl_complex_math.h>	/* GSL, Arithmetic operations for complex numbers */
#include <gsl/gsl_rng.h>			/* GSL, Random number generation */
#include <gsl/gsl_randist.h>		/* GSL, Random number generation */
#include <libconfig.h>
#include "SynthesisFunctions.h"










/**
 * Function Check_command_line
 *
 * Checks command line formant and prints instructions
 *
 * @param argc number of input arguments
 */
int Check_command_line(int argc) {

	if (argc < 3 || argc > 4) {

		/* Incorrect command line, print instructions */
		printf("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
		printf("             GlottHMM - Speech Synthesizer (%s)\n",VERSION);
		printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n\n");
		printf("Description:\n\n");
		printf("    Synthesis of speech signal according to speech parameters.\n\n");
		printf("Usage:\n\n");
		printf("    Synthesis file_name config_default config_user\n\n");
		printf("	file_name       - File name without extensions (.wav, .lab)\n");
		printf("	config_default  - Default config file name\n");
		printf("	config_user     - User config file name (OPTIONAL)\n\n");
		printf("    Synthesized speech signal is saved to \"file_name.syn.wav\".\n\n");
		printf("Version:\n\n");
		printf("    %s (%s)\n\n",VERSION,DATE);
		return EXIT_FAILURE;
	} else {

		/* Correct command line, print program title */
		printf("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
		printf("                  Speech Synthesizer (%s)\n",VERSION);
		printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n\n");
		return EXIT_SUCCESS;
	}
}




/**
 * Function Read_config
 *
 * Read configuration file
 *
 * @param filename name of the configuration file
 * @return conf pointer to configuration structure
 */
struct config_t *Read_config(char *filename) {

	struct config_t *conf = (struct config_t *)malloc(sizeof(struct config_t));
	config_init(conf);
	if(config_read_file(conf,filename) != CONFIG_TRUE) {
		printf("\nError reading configuration file \"%s\", line %i\n",filename,config_error_line(conf));
		printf("%s\n",config_error_text(conf));
		config_destroy(conf);
		free(conf);
		return NULL;
	}
	return conf;
}







/**
 * Function Assign_config_parameters
 *
 * Assign configuration parameters to variables. Destroy and free config file.
 *
 * @param conf config file
 * @param params parameter structure
 */
int Assign_config_parameters(const char *filename, struct config_t *conf, PARAM *params, int conf_type) {

	int i,bool,multisyn_old = 0,synlistlen_old = 0;
	long int ival;
	double fval;
	const char *tmp;
	const config_setting_t *paramweights_conf;
	const config_setting_t *dnn_dim_conf;

	/* For synthesis list */
	if(conf_type == DEF_CONF) {
		multisyn_old = 0;
		synlistlen_old = 0;
	} else if(conf_type == USR_CONF) {
		multisyn_old = params->multisyn;
		synlistlen_old = params->synlistlen;
	}

	/* Assign configuration file parameters */
	if(config_lookup(conf,SAMPLING_FREQUENCY) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", SAMPLING_FREQUENCY);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,SAMPLING_FREQUENCY) != NULL) {
		config_lookup_int(conf, SAMPLING_FREQUENCY,&ival);
		params->FS = (int)ival;
	}
	if(config_lookup(conf,FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,FRAME_LENGTH) != NULL) {
		config_lookup_float(conf, FRAME_LENGTH,&fval);
		params->frame_length_ms = fval;
	}
	if(config_lookup(conf,FRAME_SHIFT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", FRAME_SHIFT);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,FRAME_SHIFT) != NULL) {
		config_lookup_float(conf, FRAME_SHIFT,&fval);
		params->shift_ms = fval;
	}
	if(config_lookup(conf,F0_FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", F0_FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,F0_FRAME_LENGTH) != NULL) {
		config_lookup_float(conf, F0_FRAME_LENGTH,&fval);
		params->f0_frame_length_ms = fval;
	}
	if(config_lookup(conf,LPC_ORDER) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",LPC_ORDER);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,LPC_ORDER) != NULL) {
		config_lookup_int(conf,LPC_ORDER,&ival);
		params->lpc_order_vt = (int)ival;
	}
	if(config_lookup(conf,LPC_ORDER_SOURCE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",LPC_ORDER_SOURCE);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,LPC_ORDER_SOURCE) != NULL) {
		config_lookup_int(conf,LPC_ORDER_SOURCE,&ival);
		params->lpc_order_gl = (int)ival;
	}
	if(config_lookup(conf,WARPING_VT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",WARPING_VT);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,WARPING_VT) != NULL) {
		config_lookup_float(conf,WARPING_VT,&fval);
		params->lambda_vt = fval;
	}
	if(config_lookup(conf,WARPING_GL) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",WARPING_GL);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,WARPING_GL) != NULL) {
		config_lookup_float(conf,WARPING_GL,&fval);
		params->lambda_gl = fval;
	}
	if(config_lookup(conf,DIFFERENTIAL_LSF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", DIFFERENTIAL_LSF);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,DIFFERENTIAL_LSF) != NULL) {
		config_lookup_bool(conf, DIFFERENTIAL_LSF,&bool);
		params->differential_lsf = bool;
	}
	if(config_lookup(conf,USE_PULSE_LIBRARY) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",USE_PULSE_LIBRARY);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_PULSE_LIBRARY) != NULL) {
		config_lookup_bool(conf,USE_PULSE_LIBRARY,&bool);
		params->use_pulselib = bool;
	}
	if(config_lookup(conf,NUMBER_OF_PULSE_CANDIDATES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",NUMBER_OF_PULSE_CANDIDATES);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NUMBER_OF_PULSE_CANDIDATES) != NULL) {
		config_lookup_int(conf,NUMBER_OF_PULSE_CANDIDATES,&ival);
		params->n_pulsecandidates = (int)ival;
	}
	if(config_lookup(conf,CONCATENATION_COST) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",CONCATENATION_COST);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,CONCATENATION_COST) != NULL) {
		config_lookup_float(conf,CONCATENATION_COST,&fval);
		params->concatenation_cost = fval;
	}
	if(config_lookup(conf,PULSE_ERROR_BIAS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",PULSE_ERROR_BIAS);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,PULSE_ERROR_BIAS) != NULL) {
		config_lookup_float(conf,PULSE_ERROR_BIAS,&fval);
		params->pulse_error_bias = fval;
	}
	if(config_lookup(conf,TARGET_COST) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",TARGET_COST);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,TARGET_COST) != NULL) {
		config_lookup_float(conf,TARGET_COST,&fval);
		params->target_cost = fval;
	}
	if(config_lookup(conf,USE_PULSE_CLUSTERING) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",USE_PULSE_CLUSTERING);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_PULSE_CLUSTERING) != NULL) {
		config_lookup_bool(conf,USE_PULSE_CLUSTERING,&bool);
		params->pulse_clustering = bool;
	}
	if(config_lookup(conf,USE_PULSE_INTERPOLATION) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",USE_PULSE_INTERPOLATION);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_PULSE_INTERPOLATION) != NULL) {
		config_lookup_bool(conf,USE_PULSE_INTERPOLATION,&bool);
		params->pulse_interpolation = bool;
	}
	if(config_lookup(conf,MAX_PULSES_IN_CLUSTER) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",MAX_PULSES_IN_CLUSTER);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,MAX_PULSES_IN_CLUSTER) != NULL) {
		config_lookup_int(conf,MAX_PULSES_IN_CLUSTER,&ival);
		params->max_pulses_in_cluster = (int)ival;
	}
	if(config_lookup(conf,NOISE_GAIN_VOICED) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",NOISE_GAIN_VOICED);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NOISE_GAIN_VOICED) != NULL) {
		config_lookup_float(conf,NOISE_GAIN_VOICED,&fval);
		params->noise_gain_voiced = fval;
	}
	if(config_lookup(conf,USE_HMM) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",USE_HMM);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_HMM) != NULL) {
		config_lookup_bool(conf,USE_HMM,&bool);
		params->use_hmm = bool;
	}
	if(config_lookup(conf,POSTFILTER_COEFFICIENT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",POSTFILTER_COEFFICIENT);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,POSTFILTER_COEFFICIENT) != NULL) {
		config_lookup_float(conf,POSTFILTER_COEFFICIENT,&fval);
		params->postfilter_alpha = fval;
	}
	if(config_lookup(conf,PITCH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",PITCH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,PITCH) != NULL) {
		config_lookup_float(conf,PITCH,&fval);
		params->pitch = fval;
	}
	if(config_lookup(conf,SPEED) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",SPEED);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,SPEED) != NULL) {
		config_lookup_float(conf,SPEED,&fval);
		params->speed = fval;
	}
	if(config_lookup(conf,GAIN_UNVOICED) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",GAIN_UNVOICED);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,GAIN_UNVOICED) != NULL) {
		config_lookup_float(conf,GAIN_UNVOICED,&fval);
		params->gain_unvoiced = fval;
	}
	if(config_lookup(conf,FILTER_UPDATE_INTERVAL_VT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",FILTER_UPDATE_INTERVAL_VT);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,FILTER_UPDATE_INTERVAL_VT) != NULL) {
		config_lookup_float(conf,FILTER_UPDATE_INTERVAL_VT,&fval);
		params->filter_update_interval_vt_ms = fval;
	}
	if(config_lookup(conf,FILTER_UPDATE_INTERVAL_GL) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",FILTER_UPDATE_INTERVAL_GL);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,FILTER_UPDATE_INTERVAL_GL) != NULL) {
		config_lookup_float(conf,FILTER_UPDATE_INTERVAL_GL,&fval);
		params->filter_update_interval_gl_ms = fval;
	}
	if(config_lookup(conf,GLFLOWSP_SMOOTH_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",GLFLOWSP_SMOOTH_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,GLFLOWSP_SMOOTH_LEN) != NULL) {
		config_lookup_int(conf,GLFLOWSP_SMOOTH_LEN,&ival);
		params->glflowsp_smooth_len = (int)ival;
	}
	if(config_lookup(conf,LSF_SMOOTH_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",LSF_SMOOTH_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,LSF_SMOOTH_LEN) != NULL) {
		config_lookup_int(conf,LSF_SMOOTH_LEN,&ival);
		params->lsf_smooth_len = (int)ival;
	}
	if(config_lookup(conf,HARMONICS_SMOOTH_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",HARMONICS_SMOOTH_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,HARMONICS_SMOOTH_LEN) != NULL) {
		config_lookup_int(conf,HARMONICS_SMOOTH_LEN,&ival);
		params->harmonics_smooth_len = (int)ival;
	}
	if(config_lookup(conf,GAIN_SMOOTH_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",GAIN_SMOOTH_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,GAIN_SMOOTH_LEN) != NULL) {
		config_lookup_int(conf,GAIN_SMOOTH_LEN,&ival);
		params->gain_smooth_len = (int)ival;
	}
	if(config_lookup(conf,NORM_GAIN_SMOOTH_V_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",NORM_GAIN_SMOOTH_V_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NORM_GAIN_SMOOTH_V_LEN) != NULL) {
		config_lookup_int(conf,NORM_GAIN_SMOOTH_V_LEN,&ival);
		params->norm_gain_smooth_v_len = (int)ival;
	}
	if(config_lookup(conf,NORM_GAIN_SMOOTH_UV_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",NORM_GAIN_SMOOTH_UV_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NORM_GAIN_SMOOTH_UV_LEN) != NULL) {
		config_lookup_int(conf,NORM_GAIN_SMOOTH_UV_LEN,&ival);
		params->norm_gain_smooth_uv_len = (int)ival;
	}
	if(config_lookup(conf,GAIN_UNVOICED_FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",GAIN_UNVOICED_FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,GAIN_UNVOICED_FRAME_LENGTH) != NULL) {
		config_lookup_float(conf,GAIN_UNVOICED_FRAME_LENGTH,&fval);
		params->gain_unvoiced_frame_length_ms = fval;
	}
	if(config_lookup(conf,GAIN_VOICED_FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",GAIN_VOICED_FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,GAIN_VOICED_FRAME_LENGTH) != NULL) {
		config_lookup_float(conf,GAIN_VOICED_FRAME_LENGTH,&fval);
		params->gain_voiced_frame_length_ms = fval;
	}
	if(config_lookup(conf,HNR_SMOOTH_LEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",HNR_SMOOTH_LEN);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,HNR_SMOOTH_LEN) != NULL) {
		config_lookup_int(conf,HNR_SMOOTH_LEN,&ival);
		params->hnr_smooth_len = (int)ival;
	}
	if(config_lookup(conf,HNR_CHANNELS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", HNR_CHANNELS);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,HNR_CHANNELS) != NULL) {
		config_lookup_int(conf, HNR_CHANNELS,&ival);
		params->hnr_channels = (int)ival;
	}
	if(config_lookup(conf,NOISE_LOW_FREQ_LIMIT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", NOISE_LOW_FREQ_LIMIT);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NOISE_LOW_FREQ_LIMIT) != NULL) {
		config_lookup_float(conf, NOISE_LOW_FREQ_LIMIT,&fval);
		params->noise_low_freq_limit = fval;
	}
	if(config_lookup(conf,ADD_NOISE_PULSELIB) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", ADD_NOISE_PULSELIB);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,ADD_NOISE_PULSELIB) != NULL) {
		config_lookup_bool(conf, ADD_NOISE_PULSELIB,&bool);
		params->add_noise_pulselib = bool;
	}
	if(config_lookup(conf,USE_HARMONIC_MODIFICATION) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_HARMONIC_MODIFICATION);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_HARMONIC_MODIFICATION) != NULL) {
		config_lookup_bool(conf, USE_HARMONIC_MODIFICATION,&bool);
		params->use_harmonic_modification = bool;
	}
	if(config_lookup(conf,NUMBER_OF_HARMONICS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", NUMBER_OF_HARMONICS);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NUMBER_OF_HARMONICS) != NULL) {
		config_lookup_int(conf, NUMBER_OF_HARMONICS,&ival);
		params->number_of_harmonics = (int)ival;
	}
	if(config_lookup(conf,HP_FILTER_F0) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", HP_FILTER_F0);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,HP_FILTER_F0) != NULL) {
		config_lookup_bool(conf, HP_FILTER_F0,&bool);
		params->hpfiltf0 = bool;
	}
	if(config_lookup(conf,WAVEFORM_SAMPLES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", WAVEFORM_SAMPLES);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,WAVEFORM_SAMPLES) != NULL) {
		config_lookup_int(conf, WAVEFORM_SAMPLES,&ival);
		params->waveform_samples = (int)ival;
	}
	if(config_lookup(conf,NOISE_ROBUST_SPEECH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", NOISE_ROBUST_SPEECH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,NOISE_ROBUST_SPEECH) != NULL) {
		config_lookup_bool(conf, NOISE_ROBUST_SPEECH,&bool);
		params->noise_robust_speech = bool;
	}
	if(config_lookup(conf,SEPARATE_VOICED_UNVOICED_SPECTRUM) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", SEPARATE_VOICED_UNVOICED_SPECTRUM);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,SEPARATE_VOICED_UNVOICED_SPECTRUM) != NULL) {
		config_lookup_bool(conf, SEPARATE_VOICED_UNVOICED_SPECTRUM,&bool);
		params->sep_vuv_spectrum = bool;
	}
	if(config_lookup(conf,SYNTHESIZE_MULTIPLE_FILES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", SYNTHESIZE_MULTIPLE_FILES);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,SYNTHESIZE_MULTIPLE_FILES) != NULL) {
		config_lookup_bool(conf, SYNTHESIZE_MULTIPLE_FILES,&bool);
		params->multisyn = bool;
	}
	if(config_lookup(conf,USE_TILT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_TILT);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_TILT) != NULL) {
		config_lookup_bool(conf, USE_TILT,&bool);
		params->use_tilt = bool;
	}
	if(config_lookup(conf,USE_HNR) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_HNR);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_HNR) != NULL) {
		config_lookup_bool(conf, USE_HNR,&bool);
		params->use_hnr = bool;
	}
	if(config_lookup(conf,USE_HARMONICS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_HARMONICS);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_HARMONICS) != NULL) {
		config_lookup_bool(conf, USE_HARMONICS,&bool);
		params->use_harmonics = bool;
	}
	if(config_lookup(conf,USE_H1H2) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_H1H2);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_H1H2) != NULL) {
		config_lookup_bool(conf, USE_H1H2,&bool);
		params->use_h1h2 = bool;
	}
	if(config_lookup(conf,USE_NAQ) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_NAQ);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_NAQ) != NULL) {
		config_lookup_bool(conf, USE_NAQ,&bool);
		params->use_naq = bool;
	}
	if(config_lookup(conf,USE_WAVEFORM) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", USE_WAVEFORM);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,USE_WAVEFORM) != NULL) {
		config_lookup_bool(conf, USE_WAVEFORM,&bool);
		params->use_waveform = bool;
	}
	if(config_lookup(conf,WRITE_EXCITATION_TO_WAV) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", WRITE_EXCITATION_TO_WAV);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,WRITE_EXCITATION_TO_WAV) != NULL) {
		config_lookup_bool(conf, WRITE_EXCITATION_TO_WAV,&bool);
		params->write_excitation_to_wav = bool;
	}
	if(config_lookup(conf,JITTER) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", JITTER);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,JITTER) != NULL) {
		config_lookup_float(conf,JITTER,&fval);
		params->jitter = fval;
	}
	if (config_lookup(conf, NOISE_REDUCTION_SYNTHESIS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", NOISE_REDUCTION_SYNTHESIS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, NOISE_REDUCTION_SYNTHESIS) != NULL) {
		config_lookup_bool(conf, NOISE_REDUCTION_SYNTHESIS,&bool);
		params->noise_reduction_synthesis = bool;
	}
	if (config_lookup(conf, NOISE_REDUCTION_LIMIT_DB) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", NOISE_REDUCTION_LIMIT_DB);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, NOISE_REDUCTION_LIMIT_DB) != NULL) {
		config_lookup_float(conf, NOISE_REDUCTION_LIMIT_DB,&fval);
		params->noise_reduction_limit_db = fval;
	}
	if (config_lookup(conf, NOISE_REDUCTION_DB) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", NOISE_REDUCTION_DB);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, NOISE_REDUCTION_DB) != NULL) {
		config_lookup_float(conf, NOISE_REDUCTION_DB,&fval);
		params->noise_reduction_db = fval;
	}
	if (config_lookup(conf, NORMALIZE_PULSELIB) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", NORMALIZE_PULSELIB);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, NORMALIZE_PULSELIB) != NULL) {
		config_lookup_bool(conf, NORMALIZE_PULSELIB,&bool);
		params->normalize_pulselib = bool;
	}
	if (config_lookup(conf, ADAPT_TO_PULSELIB) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", ADAPT_TO_PULSELIB);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, ADAPT_TO_PULSELIB) != NULL) {
		config_lookup_bool(conf, ADAPT_TO_PULSELIB,&bool);
		params->adapt_to_pulselib = bool;
	}
	if (config_lookup(conf, ADAPT_COEFF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", ADAPT_COEFF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, ADAPT_COEFF) != NULL) {
		config_lookup_float(conf, ADAPT_COEFF,&fval);
		params->adapt_coeff = fval;
	}
	if (config_lookup(conf, USE_PULSELIB_LSF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_PULSELIB_LSF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_PULSELIB_LSF) != NULL) {
		config_lookup_bool(conf, USE_PULSELIB_LSF,&bool);
		params->use_pulselib_lsf = bool;
	}
	if (config_lookup(conf, AVERAGE_N_ADJACENT_PULSES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", AVERAGE_N_ADJACENT_PULSES);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, AVERAGE_N_ADJACENT_PULSES) != NULL) {
		config_lookup_int(conf, AVERAGE_N_ADJACENT_PULSES,&ival);
		params->average_n_adjacent_pulses = (int)ival;
	}
	if (config_lookup(conf, LOG_F0) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", LOG_F0);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LOG_F0) != NULL) {
		config_lookup_bool(conf, LOG_F0,&bool);
		params->logf0 = bool;
	}
	if (config_lookup(conf, USE_DNN_PULSEGEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_DNN_PULSEGEN);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_DNN_PULSEGEN) != NULL) {
		config_lookup_bool(conf, USE_DNN_PULSEGEN,&bool);
		params->use_dnn_pulsegen = bool;
	}
	if (config_lookup(conf, USE_DNN_PULSELIB_SEL) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_DNN_PULSELIB_SEL);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_DNN_PULSELIB_SEL) != NULL) {
		config_lookup_bool(conf, USE_DNN_PULSELIB_SEL,&bool);
		params->use_dnn_pulselib_sel = bool;
	}
	if (config_lookup(conf, USE_DNN_SPECMATCH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_DNN_SPECMATCH);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_DNN_SPECMATCH) != NULL) {
		config_lookup_bool(conf, USE_DNN_SPECMATCH,&bool);
		params->use_dnn_specmatch = bool;
	}
	if (config_lookup(conf, DNN_INPUT_NORMALIZED) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", DNN_INPUT_NORMALIZED);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, DNN_INPUT_NORMALIZED) != NULL) {
		config_lookup_bool(conf, DNN_INPUT_NORMALIZED,&bool);
		params->dnn_input_normalized = bool;
	}
	if (config_lookup(conf, DNN_NUMBER_OF_STACKED_FRAMES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", DNN_NUMBER_OF_STACKED_FRAMES);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, DNN_NUMBER_OF_STACKED_FRAMES) != NULL) {
		config_lookup_int(conf, DNN_NUMBER_OF_STACKED_FRAMES,&ival);
		params->dnn_number_of_stacked_frames = (int)ival;
	}
	if (config_lookup(conf, HNR_COMPENSATION) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", HNR_COMPENSATION);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, HNR_COMPENSATION) != NULL) {
		config_lookup_bool(conf, HNR_COMPENSATION,&bool);
		params->hnr_compensation = bool;
	}
	if (config_lookup(conf, UNVOICED_PRE_EMPHASIS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", UNVOICED_PRE_EMPHASIS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, UNVOICED_PRE_EMPHASIS) != NULL) {
		config_lookup_bool(conf, UNVOICED_PRE_EMPHASIS,&bool);
		params->unvoiced_pre_emphasis = bool;
	}

	/* Pulse PCA as target cost */
	if (config_lookup(conf, USE_PULSE_PCA) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_PULSE_PCA);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_PULSE_PCA) != NULL) {
		config_lookup_bool(conf, USE_PULSE_PCA,&bool);
		params->use_pulse_pca = bool;
	}

	/* Pulse library PCA */
	if (config_lookup(conf, USE_PULSELIB_PCA) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_PULSELIB_PCA);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_PULSELIB_PCA) != NULL) {
		config_lookup_bool(conf, USE_PULSELIB_PCA,&bool);
		params->use_pulselib_pca = bool;
	}
	if (config_lookup(conf, PCA_ORDER) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", PCA_ORDER);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, PCA_ORDER) != NULL) {
		config_lookup_int(conf, PCA_ORDER,&ival);
		params->pca_order = (int)ival;
	}
	if (config_lookup(conf, PCA_ORDER_SYNTHESIS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", PCA_ORDER_SYNTHESIS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, PCA_ORDER_SYNTHESIS) != NULL) {
		config_lookup_int(conf, PCA_ORDER_SYNTHESIS,&ival);
		params->pca_order_synthesis = (int)ival;
	}
	if (config_lookup(conf, PCA_SPECTRAL_MATCHING) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", PCA_SPECTRAL_MATCHING);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, PCA_SPECTRAL_MATCHING) != NULL) {
		config_lookup_bool(conf, PCA_SPECTRAL_MATCHING,&bool);
		params->pca_spectral_matching = bool;
	}
	if (config_lookup(conf, PCA_PULSE_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", PCA_PULSE_LENGTH);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, PCA_PULSE_LENGTH) != NULL) {
		config_lookup_int(conf, PCA_PULSE_LENGTH,&ival);
		params->pca_pulse_length = (int)ival;
	}
	if (config_lookup(conf, TWO_PITCH_PERIOD_DIFF_PULSE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", TWO_PITCH_PERIOD_DIFF_PULSE);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, TWO_PITCH_PERIOD_DIFF_PULSE) != NULL) {
		config_lookup_bool(conf, TWO_PITCH_PERIOD_DIFF_PULSE,&bool);
		params->two_pitch_period_diff_pulse = bool;
	}

	/* Read data format */
	if (config_lookup(conf, DATA_FORMAT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",DATA_FORMAT);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, DATA_FORMAT) != NULL) {
		config_lookup_string(conf,DATA_FORMAT,&tmp);
		if(strcmp(tmp,DATA_FORMAT_ASCII) == 0)
			params->data_format = DATA_FORMAT_ID_ASCII;
		else if(strcmp(tmp,DATA_FORMAT_BINARY) == 0)
			params->data_format = DATA_FORMAT_ID_BINARY;
		else {
			printf("\nError: Invalid configuration value \"%s\".\n", DATA_FORMAT);
			return EXIT_FAILURE;
		}
	}

	/* Read formant enhancement method */
	if (config_lookup(conf, POSTFILTER_METHOD) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",POSTFILTER_METHOD);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, POSTFILTER_METHOD) != NULL) {
		config_lookup_string(conf,POSTFILTER_METHOD,&tmp);
		if(strcmp(tmp,POSTFILTER_LSF) == 0)
			params->postfilter_method = POSTFILTER_ID_LSF;
		else if(strcmp(tmp,POSTFILTER_LPC) == 0)
			params->postfilter_method = POSTFILTER_ID_LPC;
		else if(strcmp(tmp,POSTFILTER_NONE) == 0)
			params->postfilter_method = POSTFILTER_ID_NONE;
		else {
			printf("\nError: Invalid configuration value \"%s\".\n", POSTFILTER_METHOD);
			return EXIT_FAILURE;
		}
	}

	/* Read parameter weights */
	paramweights_conf = config_lookup(conf,PARAMETER_WEIGHTS);
	if(paramweights_conf == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", PARAMETER_WEIGHTS);
		return EXIT_FAILURE;
	} else if(paramweights_conf != NULL) {
		if(conf_type == USR_CONF)
			gsl_vector_free(params->paramweights);
		params->paramweights = gsl_vector_alloc(config_setting_length(paramweights_conf));
		for(i=0;i<config_setting_length(paramweights_conf);i++)
			gsl_vector_set(params->paramweights,i,config_setting_get_float_elem(paramweights_conf, i));
	}

	/* Read DNN weight dimensions */
	dnn_dim_conf = config_lookup(conf,DNN_WEIGHT_DIMS);
	if(dnn_dim_conf == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n", DNN_WEIGHT_DIMS);
		return EXIT_FAILURE;
	} else if(dnn_dim_conf != NULL) {
		if(conf_type == USR_CONF)
			gsl_vector_free(params->dnn_weight_dims);
		params->dnn_weight_dims = gsl_vector_alloc(config_setting_length(dnn_dim_conf));
		for(i=0;i<config_setting_length(dnn_dim_conf);i++)
			gsl_vector_set(params->dnn_weight_dims,i,config_setting_get_int_elem(dnn_dim_conf, i));
	}

	/* Read pulse filename */
	if(config_lookup(conf,GLOTTAL_PULSE_NAME) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",GLOTTAL_PULSE_NAME);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,GLOTTAL_PULSE_NAME) != NULL) {
		if(conf_type == USR_CONF)
			free(params->pulse_filename);
		config_lookup_string(conf,GLOTTAL_PULSE_NAME,&tmp);
		params->pulse_filename = (char *)malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(params->pulse_filename,tmp);
	}

	/* Read pulse library filename */
	if(config_lookup(conf,PULSE_LIBRARY_NAME) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",PULSE_LIBRARY_NAME);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,PULSE_LIBRARY_NAME) != NULL) {
		if(conf_type == USR_CONF)
			free(params->pulselibrary_filename);
		config_lookup_string(conf,PULSE_LIBRARY_NAME,&tmp);
		params->pulselibrary_filename = (char *)malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(params->pulselibrary_filename,tmp);
	}

	/* Read synthesis list filename */
	if(config_lookup(conf,SYNTHESIS_LIST) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",SYNTHESIS_LIST);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,SYNTHESIS_LIST) != NULL) {
		if(conf_type == USR_CONF)
			free(params->synlist_filename);
		config_lookup_string(conf,SYNTHESIS_LIST,&tmp);
		params->synlist_filename = (char *)malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(params->synlist_filename,tmp);
	}

	/* Read DNN path */
	if(config_lookup(conf,DNN_WEIGHT_PATH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read configuration \"%s\".\n",DNN_WEIGHT_PATH);
		return EXIT_FAILURE;
	} else if(config_lookup(conf,DNN_WEIGHT_PATH) != NULL) {
		if(conf_type == USR_CONF)
			free(params->dnnpath);
		config_lookup_string(conf,DNN_WEIGHT_PATH,&tmp);
		params->dnnpath = (char *)malloc((strlen(tmp)+1)*sizeof(char));
		strcpy(params->dnnpath,tmp);
	}

	/* Free memory if previously a synthesis list was allocated */
	if(conf_type == USR_CONF) {
		for(i=0;i<synlistlen_old;i++)
			free(params->synlist[i]);
		free(params->synlist);
	}

	/* Create synthesis list */
	if(params->multisyn == 1) {
		if(ReadSynthesisList(params) == EXIT_FAILURE)
			return EXIT_FAILURE;
	} else {
		params->synlistlen = 1;
		params->synlist = (char **)malloc(sizeof(char *));
		params->synlist[0] = (char *)malloc((strlen(filename)+1)*sizeof(char));
		strcpy(params->synlist[0],filename);
	}

	/* Convert milliseconds to samples */
	params->frame_length = rint(params->FS*params->frame_length_ms/1000);
	params->shift = rint(params->FS*params->shift_ms/1000);
	params->f0_frame_length = rint(params->FS*params->f0_frame_length_ms/1000);
	params->gain_voiced_frame_length = rint(params->FS*params->gain_voiced_frame_length_ms/1000);
	params->gain_unvoiced_frame_length = rint(params->FS*params->gain_unvoiced_frame_length_ms/1000);
	params->filter_update_interval_vt = GSL_MAX(rint(params->FS*params->filter_update_interval_vt_ms/1000),1);
	params->filter_update_interval_gl = GSL_MAX(rint(params->FS*params->filter_update_interval_gl_ms/1000),1);

	/* Set time */
	if(conf_type == DEF_CONF)
		params->time_temp = (double)clock();

	/* Free memory */
	config_destroy(conf);
	free(conf);

	return EXIT_SUCCESS;
}


























/**
 * Function ReadSynthesisList
 *
 * Read list of synthesis file names
 *
 * @param name filename
 * @return number of files in the list
 */
int ReadSynthesisList(PARAM *params) {

	FILE *file;
	char s[DEF_STRING_LEN];
	int ind;

	/* Open file */
	file = fopen(params->synlist_filename, "r");
	if(!file) {
		printf("Error opening synthesis list file \"%s\": %s\n", params->synlist_filename, strerror(errno));
		return EXIT_FAILURE;
	}

	/* Read lines until EOF */
	ind = 0;
	while(fscanf(file,"%s",s) != EOF)
		ind++;

	/* Allocate memory for strings */
	params->synlistlen = ind;
	char **filenames = (char **)malloc(ind*sizeof(char *));

	/* Read file names to array */
	fseek(file, 0, SEEK_SET);
	ind = 0;
	char *fn;
	while(fscanf(file,"%s",s) != EOF) {
		fn = (char *)malloc((strlen(s)+1)*sizeof(char));
		strcpy(fn,s);
		filenames[ind] = fn;
		ind++;
	}
	fclose(file);

	/* Set list to params */
	params->synlist = filenames;
	return EXIT_SUCCESS;
}




/**
 * Read_DNN_weights
 *
 * Read DNN pulse generation weights from file, check validity
 *
 */
int Read_DNN_weights(PARAM *params, gsl_matrix **DNN_W) {

	/* Do not read if DNNs are not used */
	if(params->use_dnn_pulsegen == 0)
		return EXIT_SUCCESS;

	/* Allocate weight matrices and read data from files */
	int i,numlayers = params->dnn_weight_dims->size/2;
	char temp[DEF_STRING_LEN];
	char num[20];
	FILE *wfile;
	for(i=0;i<numlayers;i++) {
		DNN_W[i] = gsl_matrix_alloc(gsl_vector_get(params->dnn_weight_dims,2*i),gsl_vector_get(params->dnn_weight_dims,2*i+1));
		strcpy(temp,params->dnnpath);
		strcat(temp,FILENAME_DNN_W);
		sprintf(num,"%d",i+1);
		strcat(temp,num);
		wfile = fopen(temp, "r");
		if(wfile == NULL) {
			printf("Error opening file \"%s\": %s\n",temp,strerror(errno));
			return EXIT_FAILURE;
		}
		gsl_matrix_fscanf(wfile,DNN_W[i]);
		fclose(wfile);
	}

	/* Return */
	return EXIT_SUCCESS;
}


/**
 * Read_input_minmax
 *
 * Read input data minimum and maximum for normalizing the input data for DNN pulse generation
 *
 */
int Read_input_minmax(PARAM *params, gsl_vector **input_minmax) {

	/* Do not read if DNNs are not used */
	if(params->use_dnn_pulsegen == 0 || params->dnn_input_normalized == 0)
		return EXIT_SUCCESS;

	/* Allocate weight matrices and read data from files */
	char temp[DEF_STRING_LEN];
	char num[20];
	FILE *datafile;
	input_minmax[0] = gsl_vector_alloc(2*(gsl_vector_get(params->dnn_weight_dims,0)-1)/params->dnn_number_of_stacked_frames);
	strcpy(temp,params->dnnpath);
	strcat(temp,FILENAME_INPUT_MINMAX);
	datafile = fopen(temp, "r");
	if(datafile == NULL) {
		printf("Error opening file \"%s\": %s\n",temp,strerror(errno));
		return EXIT_FAILURE;
	}
	gsl_vector_fscanf(datafile,input_minmax[0]);
	fclose(datafile);

	/* Return */
	return EXIT_SUCCESS;
}






/**
 * Read_pulse_library
 *
 * Read pulse library from file, check validity
 *
 */
int Read_pulse_library(PARAM *params,gsl_matrix **pulses,gsl_matrix **pulses_rs,gsl_matrix **plsf,gsl_matrix **ptilt,gsl_matrix **pharm,
		gsl_matrix **phnr,gsl_matrix **pwaveform, gsl_matrix **pca_pc,gsl_matrix **pca_w_lib,gsl_vector **stoch_env,gsl_vector **stoch_sp,gsl_vector **pgain, gsl_vector **ph1h2,
		gsl_vector **pnaq, gsl_vector **pca_mean, gsl_vector **pulse_lengths) {

	/* Do not read if pulse library is not used */
	if(params->use_pulselib == 0) {
		params->n_pulsecandidates = 1;
		params->pulsemaxlen = 1;
		params->rspulsemaxlen = 1;
		return EXIT_SUCCESS;
	}

	/* Initialize */
	printf("	- Loading pulse library...\n");
	double time1 = (double)clock();
	char temp[DEF_STRING_LEN];

	/* Load params file for pulse library, */
	strcpy(temp, params->pulselibrary_filename);
	FILE *plparams_file = fopen(strcat(temp,FILENAME_ENDING_INFO), "r");
	if(plparams_file == NULL) {
		printf("	- Pulse library \"%s\" was not found\n",params->pulselibrary_filename);
		return EXIT_FAILURE;
	}

	/* Read pulse library parameterers */
	gsl_vector *plparams = gsl_vector_alloc(NPARAMS);
	gsl_vector_fscanf(plparams_file,plparams);
	params->number_of_pulses = gsl_vector_get(plparams,9);
	params->pulsemaxlen_ms = gsl_vector_get(plparams,10);
	params->rspulsemaxlen_ms = gsl_vector_get(plparams,11);

	/* Convert milliseconds to samples */
	params->pulsemaxlen = rint((double)params->FS*(double)params->pulsemaxlen_ms/1000.0);
	params->rspulsemaxlen = rint(params->FS*params->rspulsemaxlen_ms/1000);

	/* Check compatibility with configuration parameters */
	if(params->lpc_order_vt != (int)gsl_vector_get(plparams,3)) {printf("\nPulse library error: Different vocal tract LPC order\n\n"); return EXIT_FAILURE;}
	if(params->lpc_order_gl != (int)gsl_vector_get(plparams,4)) {printf("\nPulse library error: Different source LPC order\n\n"); return EXIT_FAILURE;}
	if(params->hnr_channels != gsl_vector_get(plparams,7)) {printf("\nPulse library error: Different number of HNR channels\n\n"); return EXIT_FAILURE;}
	if(params->number_of_harmonics != gsl_vector_get(plparams,8)) {printf("\nPulse library error: Different number of harmonics\n\n"); return EXIT_FAILURE;}
	if(params->waveform_samples != gsl_vector_get(plparams,12)) {printf("\nPulse library error: Different number of samples in Waveform\n\n"); return EXIT_FAILURE;}
	if(params->FS != gsl_vector_get(plparams,13)) {printf("\nPulse library error: Different sampling frequency\n\n"); return EXIT_FAILURE;}
	if(params->data_format != (int)gsl_vector_get(plparams,14)) {printf("\nPulse library error: Different data format (ascii/binary)\n\n"); return EXIT_FAILURE;}
	if(params->number_of_pulses < params->n_pulsecandidates) {
		printf("	- Warning: Number of pulse candidates is greater than the total number of pulses.\n");
		printf("	       Number of pulse candidates changed to: %i.\n",params->number_of_pulses);
		params->n_pulsecandidates = params->number_of_pulses;
	}
	gsl_vector_free(plparams);
	fclose(plparams_file);

	/* Allocate space for pulse library data */
	*(pulses) = gsl_matrix_alloc(params->number_of_pulses,params->pulsemaxlen);
	*(pulses_rs) = gsl_matrix_alloc(params->number_of_pulses,params->rspulsemaxlen);
	*(pulse_lengths) = gsl_vector_alloc(params->number_of_pulses);
	*(plsf) = gsl_matrix_calloc(params->number_of_pulses,params->lpc_order_vt);
	*(pgain) = gsl_vector_calloc(params->number_of_pulses);
	if(params->use_tilt == 1) *(ptilt) = gsl_matrix_calloc(params->number_of_pulses,params->lpc_order_gl);
	else *(ptilt) = NULL;
	if(params->use_harmonics == 1) *(pharm) = gsl_matrix_calloc(params->number_of_pulses,params->number_of_harmonics);
	else *(pharm) = NULL;
	if(params->use_hnr == 1) *(phnr) = gsl_matrix_calloc(params->number_of_pulses,params->hnr_channels);
	else *(phnr) = NULL;
	if(params->use_waveform == 1) *(pwaveform) = gsl_matrix_calloc(params->number_of_pulses,params->waveform_samples);
	else *(pwaveform) = NULL;
	if(params->use_h1h2 == 1) *(ph1h2) = gsl_vector_calloc(params->number_of_pulses);
	else *(ph1h2) = NULL;
	if(params->use_naq == 1) *(pnaq) = gsl_vector_calloc(params->number_of_pulses);
	else *(pnaq) = NULL;

	/* Load pulse data */
	strcpy(temp, params->pulselibrary_filename);
	FILE *pulses_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_PULSES), "r");
	if(pulses_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	strcpy(temp, params->pulselibrary_filename);FILE *pulses_rs_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_RSPULSES), "r");
	if(pulses_rs_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	strcpy(temp, params->pulselibrary_filename);FILE *pulse_lengths_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_PULSELENGTHS), "r");
	if(pulse_lengths_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII) {
		gsl_matrix_fscanf(pulses_file,*(pulses));
		gsl_matrix_fscanf(pulses_rs_file,*(pulses_rs));
		gsl_vector_fscanf(pulse_lengths_file,*(pulse_lengths));
	} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
		gsl_matrix_fread(pulses_file,*(pulses));
		gsl_matrix_fread(pulses_rs_file,*(pulses_rs));
		gsl_vector_fread(pulse_lengths_file,*(pulse_lengths));
	}
	fclose(pulses_file);
	fclose(pulses_rs_file);
	fclose(pulse_lengths_file);

	/* Load parameter data */
	strcpy(temp, params->pulselibrary_filename);FILE *plsf_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_LSF), "r");
	if(plsf_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_matrix_fscanf(plsf_file,*(plsf));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(plsf_file,*(plsf));
	fclose(plsf_file);

	strcpy(temp, params->pulselibrary_filename);FILE *pgain_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_GAIN), "r");
	if(pgain_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_vector_fscanf(pgain_file,*(pgain));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_vector_fread(pgain_file,*(pgain));
	fclose(pgain_file);

	if(params->use_tilt == 1) {strcpy(temp, params->pulselibrary_filename);FILE *ptilt_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_TILT), "r");
	if(ptilt_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII)	gsl_matrix_fscanf(ptilt_file,*(ptilt));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(ptilt_file,*(ptilt));
	fclose(ptilt_file);}

	if(params->use_harmonics == 1) {strcpy(temp, params->pulselibrary_filename);FILE *pharm_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_HARM), "r");
	if(pharm_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII)	gsl_matrix_fscanf(pharm_file,*(pharm));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(pharm_file,*(pharm));
	fclose(pharm_file);}

	if(params->use_hnr == 1) {strcpy(temp, params->pulselibrary_filename);FILE *phnr_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_HNR), "r");
	if(phnr_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII)	gsl_matrix_fscanf(phnr_file,*(phnr));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(phnr_file,*(phnr));
	fclose(phnr_file);}

	if(params->use_waveform == 1) {strcpy(temp, params->pulselibrary_filename);FILE *pwaveform_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_WAVEFORM), "r");
	if(pwaveform_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII)	gsl_matrix_fscanf(pwaveform_file,*(pwaveform));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(pwaveform_file,*(pwaveform));
	fclose(pwaveform_file);}

	if(params->use_h1h2 == 1) {strcpy(temp, params->pulselibrary_filename);FILE *ph1h2_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_H1H2), "r");
	if(ph1h2_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII)	gsl_vector_fscanf(ph1h2_file,*(ph1h2));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_vector_fread(ph1h2_file,*(ph1h2));
	fclose(ph1h2_file);}

	if(params->use_naq == 1) {strcpy(temp, params->pulselibrary_filename);FILE *pnaq_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_NAQ), "r");
	if(pnaq_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->data_format == DATA_FORMAT_ID_ASCII)	gsl_vector_fscanf(pnaq_file,*(pnaq));
	else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_vector_fread(pnaq_file,*(pnaq));
	fclose(pnaq_file);}

	/* Normalize the sum of paramweights to one */
	double sum = 0;
	int i;
	for(i=0;i<params->paramweights->size;i++)
		sum += gsl_vector_get(params->paramweights,i);
	for(i=0;i<params->paramweights->size;i++)
		gsl_vector_set(params->paramweights,i,gsl_vector_get(params->paramweights,i)/sum);

	/* Allocate space for pulse library PCA parameters */
	if(params->use_pulselib_pca == 1) {
		*(pca_pc) = gsl_matrix_calloc(params->pca_pulse_length,params->pca_order);
		*(pca_mean) = gsl_vector_calloc(params->pca_pulse_length);
		*(stoch_env) = gsl_vector_calloc(params->rspulsemaxlen);
		*(stoch_sp) = gsl_vector_calloc(STOCH_SP_LPC_ORDER);

		strcpy(temp, params->pulselibrary_filename);FILE *pca_pc_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_PCA_PC), "r");
		if(pca_pc_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
		if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_matrix_fscanf(pca_pc_file,*(pca_pc));
		else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(pca_pc_file,*(pca_pc));
		fclose(pca_pc_file);

		strcpy(temp, params->pulselibrary_filename);FILE *pca_mean_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_PCA_MEAN), "r");
		if(pca_mean_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
		if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_vector_fscanf(pca_mean_file,*(pca_mean));
		else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_vector_fread(pca_mean_file,*(pca_mean));
		fclose(pca_mean_file);

		//strcpy(temp, params->pulselibrary_filename);FILE *stoch_env_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_STOCH_ENV), "r");
		//if(stoch_env_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
		//if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_vector_fscanf(stoch_env_file,*(stoch_env));
		//else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_vector_fread(stoch_env_file,*(stoch_env));
		//fclose(stoch_env_file);

		//strcpy(temp, params->pulselibrary_filename);FILE *stoch_sp_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_STOCH_SP), "r");
		//if(stoch_sp_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
		//if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_vector_fscanf(stoch_sp_file,*(stoch_sp));
		//else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_vector_fread(stoch_sp_file,*(stoch_sp));
		//fclose(stoch_sp_file);
	} else {
		*(pca_pc) = NULL;
		*(pca_mean) = NULL;
		*(stoch_env) = NULL;
		*(stoch_sp) = NULL;
	}

	/* Allocate space for pulse PCA parameters */
	if(params->use_pulse_pca == 1) {
		*(pca_w_lib) = gsl_matrix_calloc(params->number_of_pulses,params->pca_order);
		strcpy(temp, params->pulselibrary_filename);FILE *pca_w_lib_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_PCA_W), "r");
		if(pca_w_lib_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
		if(params->data_format == DATA_FORMAT_ID_ASCII) gsl_matrix_fscanf(pca_w_lib_file,*(pca_w_lib));
		else if(params->data_format == DATA_FORMAT_ID_BINARY) gsl_matrix_fread(pca_w_lib_file,*(pca_w_lib));
		fclose(pca_w_lib_file);
	} else {
		*(pca_w_lib) = NULL;
	}

	/* Print elapsed time */
	printf("	  (%1.2lf s)\n",((double)clock()-time1)/(double)CLOCKS_PER_SEC);
	params->time_temp = (double)clock();
	return EXIT_SUCCESS;
}







/**
 * Function ReadPulseFile
 *
 * Read pulse values from file
 *
 * @param name filename
 * @return vector containing the values
 */
gsl_vector *ReadPulseFile(PARAM *params) {

	FILE *file;
	double values[1000];
	char s[50];
	int i,ind;

	/* Open file */
	file = fopen(params->pulse_filename, "r");
	if(!file) {
		printf("Error opening pulse file \"%s\": %s\n", params->pulse_filename, strerror(errno));
		return NULL;
	}

	/* Read lines until EOF */
	ind = 0;
	while(fscanf(file,"%s",s) != EOF) {
		values[ind] = atof(s);
		ind++;
	}

	/* Copy values to vector and return vector pointer */
	gsl_vector *pulse = gsl_vector_alloc(ind);
	for(i=0;i<ind;i++) {
		gsl_vector_set(pulse,i,values[i]);
	}
	fclose(file);
	return pulse;
}




/**
 * Function Initialize_params
 *
 * Initialize params file
 *
 * @param params
 * @return synfilenumber
 */
int Initialize_params(PARAM *params, int synfilenumber) {

	/* Initialize parameters for synthesizing new file */
	params->synfilenumber = synfilenumber;
	params->resynth = 0;
	params->compcoeff = 1.0;
	params->hnr_reestimated = 0;
	params->pulse_tilt_decrease_coeff = 1;

	/* Modification for noise robust speech */
	Noise_robust_speech1(params);

	/* Count the number of frames */
	char temp[DEF_STRING_LEN];
	strcpy(temp, params->synlist[params->synfilenumber]);
	params->n_frames = EvalFileLength(strcat(temp,FILENAME_ENDING_F0),params);

	/* Evaluate signal length */
	params->signal_length = rint(params->n_frames*params->shift/params->speed);

	// For replicating the exact analysis-synthesis file length
	//int empty_frames = rint((params->frame_length/(double)params->shift - 1)/2);
	//int frames_orig = params->n_frames - 2*empty_frames;
	//params->signal_length = (frames_orig + params->frame_length/(double)params->shift/params->speed - 1.0)*(double)params->shift/params->speed;

	/* Return */
	if(params->n_frames == -1) {
		return EXIT_FAILURE;
	}
	if(params->n_frames == 0) {
		printf("\nError: Zero length F0 vector\n\n");
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}








/**
 * Function EvalFileLength
 *
 * If file is in ASCII mode, read file and count the number of lines.
 * If file is in BINARY mode, read file size.
 *
 * @param name filename
 * @return number of parameters
 */
int EvalFileLength(const char *name, PARAM *params) {

	FILE *file;
	char s[50];
	int fileSize = 0;

	/* Open file */
	file = fopen(name, "r");
	if(!file) {
		printf("Error opening file \"%s\": %s\n", name, strerror(errno));
		return -1;
	}

	/* Read lines until EOF */
	if(params->data_format == DATA_FORMAT_ID_ASCII) {
		while(fscanf(file,"%s",s) != EOF)
			fileSize++;
	} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
		fseek(file, 0, SEEK_END);
		fileSize = ftell(file)/sizeof(double);
	}
	fclose(file);

	return fileSize;
}






/**
 * Function Print_synthesis_settings_start
 *
 * Print synthesis settings in the beginning
 *
 * @param params
 */
void Print_synthesis_settings_start(PARAM *params) {

	printf("Synthesis of %s\n",params->synlist[0]);
}



/**
 * Function Print_synthesis_settings_middle
 *
 * Print synthesis settings in progress
 *
 * @param params
 */
void Print_synthesis_settings_middle(PARAM *params) {

	if(params->synfilenumber > 0) {
		printf("Synthesis of %s\n",params->synlist[params->synfilenumber]);
	}
}




/**
 * Function Print_synthesis_settings_end
 *
 * Print synthesis settings in the end
 *
 * @param filename
 * @param params

 */
void Print_synthesis_settings_end(PARAM *params, double time1) {

	/* Speech file name */
	char temp[DEF_STRING_LEN];
	strcpy(temp, params->synlist[params->synfilenumber]);
	strcat(temp,FILENAME_ENDING_SYNTHESIS);

	/* Print at the end of synthesis of one file */
	if(params->multisyn == 1) {
		printf("	  (%1.2lf s)\n",((double)clock()-params->time_temp)/(double)CLOCKS_PER_SEC);
		printf("	- Finished synthesis of %s\n\n",temp);
	}

	/* Print at the of the program */
	if(params->synfilenumber == params->synlistlen-1) {
		if(params->multisyn == 0) {
			printf("	  (%1.2lf s)\n",((double)clock()-params->time_temp)/(double)CLOCKS_PER_SEC);
			printf("	- Finished synthesis.\n\n");
			printf("Elapsed time: %1.2lf seconds.\n",((double)clock()-time1)/(double)CLOCKS_PER_SEC);
			printf("Synthesized speech saved to file \"%s\".\n\n",temp);
		} else {
			printf("Finished synthesis of files:\n\n");
			int i;
			for(i=0;i<params->synlistlen;i++)
				printf("   %s\n",params->synlist[i]);
			double t = ((double)clock()-time1)/(double)CLOCKS_PER_SEC;
			printf("\nTotal elapsed time: %1.2lf seconds (average %1.2lf seconds per file).\n\n",t,t/(double)params->synlistlen);
		}
	}
}





/**
 * Function Compatibility_check
 *
 * Check the compatibility with configuration and pulse library parameters
 *
 * @param params
 */
int Compatibility_check(PARAM *params) {

	/* Initialize */
	char temp[DEF_STRING_LEN];
	FILE *parameter_file;

	/* Open parameter file */
	strcpy(temp, params->synlist[params->synfilenumber]);
	parameter_file = fopen(strcat(temp,FILENAME_ENDING_INFO), "r");

	/* If parameter file is found, check compatibility */
	if(parameter_file == NULL) {
		if(params->use_hmm == 0)
			printf("\"%s\" was not found, proceed without compatibility check.\n\n",temp);
		if(parameter_file != NULL)
			fclose(parameter_file);
	} else {
		gsl_vector *parameters = gsl_vector_alloc(NPARAMS);
		gsl_vector_fscanf(parameter_file,parameters);
		if(params->frame_length_ms != gsl_vector_get(parameters,0)) {printf("\nWarning: Different frame length in analysis and synthesis\n\n");}
		if(params->shift_ms != gsl_vector_get(parameters,1)) {printf("\nWarning: Different shift length in analysis and synthesis\n\n");}
		if(params->n_frames != (int)gsl_vector_get(parameters,2)) {printf("\nWarning: infofile indicates a different number of frames compared to actual files\n\n");}
		if(params->n_frames < 1) {printf("\nError: Zero-length input parameters\n\n");return EXIT_FAILURE;}
		if(params->lpc_order_vt != (int)gsl_vector_get(parameters,3)) {printf("\nError: Different vocal tract LPC order in analysis and synthesis\n\n"); return EXIT_FAILURE;}
		if(params->lpc_order_gl != (int)gsl_vector_get(parameters,4)) {printf("\nError: Different source LPC order in analysis and synthesis\n\n"); return EXIT_FAILURE;}
		if(params->lambda_gl != gsl_vector_get(parameters,6)) {printf("\nError: Different source warping coefficient in analysis and synthesis\n\n"); return EXIT_FAILURE;}
		if(params->hnr_channels != gsl_vector_get(parameters,7)) {printf("\nError: Different number of HNR channels in analysis and synthesis\n\n"); return EXIT_FAILURE;}
		if(params->number_of_harmonics != gsl_vector_get(parameters,8)) {printf("\nError: Different number of harmonics in analysis and synthesis\n\n"); return EXIT_FAILURE;}
		if(params->lambda_vt != gsl_vector_get(parameters,5)) {printf("\nWarning: Different vocal tract warping coefficient in analysis and synthesis\n\n");}
		if(params->FS != gsl_vector_get(parameters,13)) {printf("\nWarning: Different sampling frequency in analysis and synthesis\n\n");}
		if(params->gain_voiced_frame_length_ms > params->frame_length_ms || params->gain_unvoiced_frame_length_ms > params->frame_length_ms) {printf("\nError: Voiced/unvoiced frame length must be less or equal to frame length\n\n"); return EXIT_FAILURE;}
		if(params->data_format != (int)gsl_vector_get(parameters,14)) {printf("\nError: Different data format (ascii/binary)\n\n"); return EXIT_FAILURE;}
		gsl_vector_free(parameters);
		fclose(parameter_file);
	}

	/* Check compatibility with used parameters and synthesis method */
	if(params->use_pulselib == 0) {
		if(params->use_tilt == 0) {printf("\nError: Spectral tilt (LSFsource) must be used with single pulse technique.\n\n");return EXIT_FAILURE;}
		if(params->use_hnr == 0) {printf("\nError: Harmonic to noise ratio (HNR) must be used with single pulse technique.\n\n");return EXIT_FAILURE;}
		if(params->use_harmonic_modification == 1 && params->use_harmonics == 0) {printf("\nError: Harmonics must be used with harmonic modification of the pulse.\n\n");return EXIT_FAILURE;}
	}

	/* Check compatibility with PCA/ICA */
	if(params->use_pulselib_pca == 1 && params->use_pulselib == 0) {printf("\nError: Pulse library must be used with PCA/ICA based pulse reconstruction.\n\n");return EXIT_FAILURE;}

	return EXIT_SUCCESS;
}





/**
 * Function Allocate_params
 *
 * Allocate synthesis parameters
 *
 * @param params
 */
void Allocate_params(gsl_vector **excitation_voiced, gsl_vector **excitation_unvoiced, gsl_vector **resynthesis_pulse_index,
		gsl_vector **gain_new, gsl_matrix **glflowsp_new, gsl_matrix **hnr_new, PARAM *params) {

	/* Allocate */
	*(excitation_voiced) = gsl_vector_calloc(params->signal_length);
	*(excitation_unvoiced) = gsl_vector_calloc(params->signal_length);
	*(gain_new) = gsl_vector_calloc(params->n_frames);
	*(hnr_new) = gsl_matrix_calloc(params->n_frames,params->hnr_channels);
	*(glflowsp_new) = gsl_matrix_calloc(params->n_frames,params->lpc_order_gl);
	*(resynthesis_pulse_index) = gsl_vector_calloc(params->n_frames);
}









/**
 * Function Read_synthesis_parameters
 *
 * Allocate and read synthesis parameters
 *
 * @param params
 */
int Read_synthesis_parameters(gsl_vector **gain, gsl_vector **fundf, gsl_matrix **LSF, gsl_matrix **LSF2, gsl_matrix **glflowsp,
		gsl_matrix **hnr, gsl_matrix **harmonics, gsl_matrix **waveform, gsl_vector **h1h2, gsl_vector **naq, gsl_matrix **pca_w, PARAM *params) {

	/* Initialize */
	char temp[DEF_STRING_LEN];

	/* Define variables */
	FILE *LSF_file = NULL,*LSF2_file = NULL,*Gain_file = NULL,*F0_file = NULL,*LSFsource_file = NULL,*hnr_file = NULL;
	FILE *harmonics_file = NULL,*waveform_file = NULL,*h1h2_file = NULL,*naq_file = NULL;

	/* Allocate memory for variables (affected by differential LSFs processing) */
	if(params->differential_lsf == 1)
		*(LSF) = gsl_matrix_alloc(params->n_frames,params->lpc_order_vt+1);
	else
		*(LSF) = gsl_matrix_alloc(params->n_frames,params->lpc_order_vt);
	if(params->use_tilt == 1)
		if(params->differential_lsf == 1)
			*(glflowsp) = gsl_matrix_alloc(params->n_frames,params->lpc_order_gl+1);
		else
			*(glflowsp) = gsl_matrix_alloc(params->n_frames,params->lpc_order_gl);
	else
		*(glflowsp) = NULL;
	if(params->sep_vuv_spectrum == 1)
		if(params->differential_lsf == 1)
			*(LSF2) = gsl_matrix_alloc(params->n_frames,params->lpc_order_vt+1);
		else
			*(LSF2) = gsl_matrix_alloc(params->n_frames,params->lpc_order_vt);
	else
		*(LSF2) = NULL;

	/* Allocate memory for variables */
	*(gain) = gsl_vector_alloc(params->n_frames);
	*(fundf) = gsl_vector_alloc(params->n_frames);
	if(params->use_hnr == 1) *(hnr) = gsl_matrix_alloc(params->n_frames,params->hnr_channels);
	else *(hnr) = NULL;
	if(params->use_harmonics == 1) *(harmonics) = gsl_matrix_alloc(params->n_frames,params->number_of_harmonics);
	else *(harmonics) = NULL;
	if(params->use_waveform == 1) *(waveform) = gsl_matrix_alloc(params->n_frames,params->waveform_samples);
	else *(waveform) = NULL;
	if(params->use_h1h2 == 1) *(h1h2) = gsl_vector_alloc(params->n_frames);
	else *(h1h2) = NULL;
	if(params->use_naq == 1) *(naq) = gsl_vector_alloc(params->n_frames);
	else *(naq) = NULL;


	/* Open files */
	strcpy(temp, params->synlist[params->synfilenumber]);LSF_file = fopen(strcat(temp,FILENAME_ENDING_LSF), "r");
	if(LSF_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	strcpy(temp, params->synlist[params->synfilenumber]);Gain_file = fopen(strcat(temp,FILENAME_ENDING_GAIN), "r");
	if(Gain_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	strcpy(temp, params->synlist[params->synfilenumber]);F0_file = fopen(strcat(temp,FILENAME_ENDING_F0), "r");
	if(F0_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
	if(params->use_tilt == 1) {strcpy(temp, params->synlist[params->synfilenumber]);LSFsource_file = fopen(strcat(temp,FILENAME_ENDING_LSFSOURCE), "r");
		if(LSFsource_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}}
	if(params->use_hnr == 1) {strcpy(temp, params->synlist[params->synfilenumber]);hnr_file = fopen(strcat(temp,FILENAME_ENDING_HNR), "r");
		if(hnr_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}}
	if(params->use_harmonics == 1) {strcpy(temp, params->synlist[params->synfilenumber]);harmonics_file = fopen(strcat(temp,FILENAME_ENDING_HARMONICS), "r");
		if(harmonics_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}}
	if(params->use_waveform == 1) {strcpy(temp, params->synlist[params->synfilenumber]);waveform_file = fopen(strcat(temp,FILENAME_ENDING_WAVEFORM), "r");
		if(waveform_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}}
	if(params->use_h1h2 == 1) {strcpy(temp, params->synlist[params->synfilenumber]);h1h2_file = fopen(strcat(temp,FILENAME_ENDING_H1H2), "r");
		if(h1h2_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}}
	if(params->use_naq == 1) {strcpy(temp, params->synlist[params->synfilenumber]);naq_file = fopen(strcat(temp,FILENAME_ENDING_NAQ), "r");
		if(naq_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}}
	if(params->sep_vuv_spectrum == 1) {strcpy(temp, params->synlist[params->synfilenumber]);LSF2_file = fopen(strcat(temp,FILENAME_ENDING_LSF2), "r");
		if(LSF2_file==NULL){printf("Error opening file \"%s\": %s.\nSet configuration parameter SEPARATE_VOICED_UNVOICED_SPECTRUM to zero if single spectrum is used.\n",temp,strerror(errno));return EXIT_FAILURE;}}

	/* Read parameters from files */
	if(params->data_format == DATA_FORMAT_ID_ASCII) {
		gsl_vector_fscanf(Gain_file,*(gain));
		gsl_vector_fscanf(F0_file,*(fundf));
		gsl_matrix_fscanf(LSF_file,*(LSF));
		if(params->use_tilt == 1) gsl_matrix_fscanf(LSFsource_file,*(glflowsp));
		if(params->use_hnr == 1) gsl_matrix_fscanf(hnr_file,*(hnr));
		if(params->use_harmonics == 1) gsl_matrix_fscanf(harmonics_file,*(harmonics));
		if(params->use_waveform == 1) gsl_matrix_fscanf(waveform_file,*(waveform));
		if(params->use_h1h2 == 1) gsl_vector_fscanf(h1h2_file,*(h1h2));
		if(params->use_naq == 1) gsl_vector_fscanf(naq_file,*(naq));
		if(params->sep_vuv_spectrum == 1) gsl_matrix_fscanf(LSF2_file,*(LSF2));
	} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
		gsl_vector_fread(Gain_file,*(gain));
		gsl_vector_fread(F0_file,*(fundf));
		gsl_matrix_fread(LSF_file,*(LSF));
		if(params->use_tilt == 1) gsl_matrix_fread(LSFsource_file,*(glflowsp));
		if(params->use_hnr == 1) gsl_matrix_fread(hnr_file,*(hnr));
		if(params->use_harmonics == 1) gsl_matrix_fread(harmonics_file,*(harmonics));
		if(params->use_waveform == 1) gsl_matrix_fread(waveform_file,*(waveform));
		if(params->use_h1h2 == 1) gsl_vector_fread(h1h2_file,*(h1h2));
		if(params->use_naq == 1) gsl_vector_fread(naq_file,*(naq));
		if(params->sep_vuv_spectrum == 1) gsl_matrix_fread(LSF2_file,*(LSF2));
	}

	/* Close files */
	fclose(F0_file);
	fclose(LSF_file);
	fclose(Gain_file);
	if(params->use_tilt == 1) fclose(LSFsource_file);
	if(params->use_hnr == 1) fclose(hnr_file);
	if(params->use_harmonics == 1)fclose(harmonics_file);
	if(params->use_waveform == 1) fclose(waveform_file);
	if(params->use_h1h2 == 1) fclose(h1h2_file);
	if(params->use_naq == 1) fclose(naq_file);
	if(params->sep_vuv_spectrum == 1) fclose(LSF2_file);

	/* Read pulse library PCA weights */
	if(params->use_pulselib_pca == 1 || params->use_pulse_pca == 1) {

		/* Initialize */
		FILE *pca_w_file;
		*(pca_w) = gsl_matrix_calloc(params->n_frames,params->pca_order);

		/* Load only if PC weights are used in synthesis of the mean pulse */
		if(params->use_pulselib_pca == 1 && params->pca_order_synthesis > 0 && params->use_dnn_pulsegen == 0) {
			strcpy(temp, params->synlist[params->synfilenumber]);
			pca_w_file = fopen(strcat(temp,FILENAME_ENDING_PULSELIB_PCA_W), "r");
			if(pca_w_file==NULL) {
				printf("Error opening file \"%s\": %s\n",temp,strerror(errno));
				return EXIT_FAILURE;
			}
			if(params->data_format == DATA_FORMAT_ID_ASCII) {
				gsl_matrix_fscanf(pca_w_file,*(pca_w));
			} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
				gsl_matrix_fread(pca_w_file,*(pca_w));
			}
			fclose(pca_w_file);
		}
	}
	return EXIT_SUCCESS;
}






/**
 * Function Print_elapsed_time
 *
 * Print elapsed time
 *
 * @param
 */
void Print_elapsed_time(PARAM *params) {

	printf("	  (%1.2lf s)\n",((double)clock()-params->time_temp)/(double)CLOCKS_PER_SEC);
	params->time_temp = (double)clock();
}


/**
 * Function Pulse_clustering
 *
 * Read file and count the number of lines
 *
 * @param pulse_clus_id pointer to pulse cluster ids
 * @param pulse_clusters pointer to pulse clusters
 * @param params paramteres structure
 * @param synfilenumber
 * @return
 */
int Pulse_clustering(gsl_vector **pulse_clus_id_POINTER, gsl_matrix **pulse_clusters_POINTER, PARAM *params) {

	if(params->pulse_clustering == 1) {

		/* Initialize */
		gsl_vector *pulse_clus_id;
		gsl_matrix *pulse_clusters;
		FILE *pulse_clus_id_file;
		char temp[DEF_STRING_LEN];
		int i,j;

		/* Read cluster ids per frame */
		strcpy(temp, params->synlist[params->synfilenumber]);pulse_clus_id_file = fopen(strcat(temp,FILENAME_ENDING_PULSECLUSTER), "r");
		if(pulse_clus_id_file==NULL){printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}
		pulse_clus_id = gsl_vector_alloc(params->n_frames);
		if(params->data_format == DATA_FORMAT_ID_ASCII)
			gsl_vector_fscanf(pulse_clus_id_file, pulse_clus_id);
		else if(params->data_format == DATA_FORMAT_ID_BINARY)
			gsl_vector_fread(pulse_clus_id_file, pulse_clus_id);
		fclose(pulse_clus_id_file);

		/* Read clusters that are used in utterance */
		pulse_clusters = gsl_matrix_alloc((int)(gsl_vector_max(pulse_clus_id)+1), params->max_pulses_in_cluster);
		gsl_matrix_set_all(pulse_clusters, -1);
		for(i=0;i<pulse_clus_id->size;i++) {

			sprintf(temp, "%s_clusters/%d", params->pulselibrary_filename, (int)gsl_vector_get(pulse_clus_id, i));
			int n_p = EvalFileLength(temp,params);
			FILE *cluster_file = fopen(temp, "r");
			if(cluster_file==NULL) {printf("Error opening file \"%s\": %s\n",temp,strerror(errno));return EXIT_FAILURE;}

			/* Allocate correctly sized vector */
			gsl_vector *temp_v;
			if(n_p > 0) {
				temp_v = gsl_vector_alloc(n_p);
				cluster_file = fopen(temp, "r");
				if(params->data_format == DATA_FORMAT_ID_ASCII)
					gsl_vector_fscanf(cluster_file, temp_v);
				else if(params->data_format == DATA_FORMAT_ID_BINARY)
					gsl_vector_fread(cluster_file, temp_v);
				fclose(cluster_file);

				/* Add clusters to matrix (terribly sparse) */
				for(j=0;j<temp_v->size;j++) {
					if(j == params->max_pulses_in_cluster) break;
					gsl_matrix_set(pulse_clusters, gsl_vector_get(pulse_clus_id, i), j,  gsl_vector_get(temp_v, j));
				}
				gsl_vector_free(temp_v);
			}
		}

		/* Set pointers */
		*(pulse_clus_id_POINTER) = pulse_clus_id;
		*(pulse_clusters_POINTER) = pulse_clusters;
		return EXIT_SUCCESS;
	}
	return EXIT_SUCCESS;
}





/**
 * Function LSF_fix_vector
 *
 * Check the validity of LSF and fix found errors (vector)
 *
 * @param lsf vector
 */
void LSF_fix_vector(gsl_vector *lsf) {

	int i,ok = 0;
	double mean;
	int flag_neg = 0;
	int flag_zero = 0;
	int flag_pi = 0;
	int flag_piplus = 0;
	int flag_nan = 0;
	int flag_dec = 0;
	int flag_close = 0;

	/* Repeat until LSF is fixed */
	while(ok == 0) {

		/* Set ok */
		ok = 1;

		/* Check and correct values less than zero or greater than pi */
		for(i=0;i<lsf->size;i++) {
			if(gsl_vector_get(lsf,i) < 0) {
				gsl_vector_set(lsf,i,LSF_EPSILON);
				flag_neg++;
				ok = 0;
			} else if(gsl_vector_get(lsf,i) < LSF_EPSILON) {
				gsl_vector_set(lsf,i,LSF_EPSILON);
				flag_zero++;
				ok = 0;
			} else if(gsl_vector_get(lsf,i) > M_PI) {
				gsl_vector_set(lsf,i,M_PI-LSF_EPSILON);
				flag_piplus++;
				ok = 0;
			} else if(gsl_vector_get(lsf,i) > M_PI-LSF_EPSILON) {
				gsl_vector_set(lsf,i,M_PI-LSF_EPSILON);
				flag_pi++;
				ok = 0;
			}
			if(gsl_isnan(gsl_vector_get(lsf,i))) {
				if(i == 0)
					gsl_vector_set(lsf,i,LSF_EPSILON);
				else if(i == lsf->size-1)
					gsl_vector_set(lsf,i,M_PI-LSF_EPSILON);
				else
					gsl_vector_set(lsf,i,(gsl_vector_get(lsf,i-1)+gsl_vector_get(lsf,i+1))/2.0);
				flag_nan++;
				ok = 0;
			}
		}

		/* Check and correct non-increasing values or coefficients too close */
		for(i=0;i<lsf->size-1;i++) {
			if(gsl_vector_get(lsf,i) > gsl_vector_get(lsf,i+1)) {
				mean = (gsl_vector_get(lsf,i)+gsl_vector_get(lsf,i+1))/2.0;
				gsl_vector_set(lsf,i,mean - LSF_EPSILON/2.0);
				gsl_vector_set(lsf,i+1,mean + LSF_EPSILON/2.0);
				flag_dec++;
				ok = 0;
			} else if(gsl_vector_get(lsf,i) > gsl_vector_get(lsf,i+1)-LSF_EPSILON/10.0) {
				mean = (gsl_vector_get(lsf,i)+gsl_vector_get(lsf,i+1))/2.0;
				gsl_vector_set(lsf,i,mean - LSF_EPSILON/2.0);
				gsl_vector_set(lsf,i+1,mean + LSF_EPSILON/2.0);
				flag_close++;
				ok = 0;
			}
		}
	}

	/* Report */
	if(flag_dec > 0)
		printf("Warning: %i decreasing LSFs -> Fixed!\n",flag_dec);
	if(flag_neg > 0)
		printf("Warning: %i negative LSFs -> Fixed!\n",flag_neg);
	if(flag_piplus > 0)
		printf("Warning: %i LSFs greater than Pi -> Fixed!\n",flag_piplus);
	if(flag_zero > 0)
		printf("Warning: %i LSFs too close to 0 -> Fixed!\n",flag_zero);
	if(flag_piplus > 0)
		printf("Warning: %i LSFs greater than Pi -> Fixed!\n",flag_pi);
	if(flag_close > 0)
		printf("Warning: %i LSFs too close to each other -> Fixed!\n",flag_close);
	if(flag_nan > 0)
		printf("Warning: %i NaN LSF value(s) -> Fixed!\n",flag_nan);
}





/**
 * Function LSF_fix_matrix
 *
 * Check the validity of LSF and fix found errors
 *
 * @param lsf matrix
 */
void LSF_fix_matrix(gsl_matrix *lsf) {

	int n,i,ok = 0;
	double mean;
	int flag_neg = 0;
	int flag_zero = 0;
	int flag_pi = 0;
	int flag_piplus = 0;
	int flag_nan = 0;
	int flag_dec = 0;
	int flag_close = 0;

	/* Repeat until LSF is fixed */
	while(ok == 0) {

		/* Set ok */
		ok = 1;

		/* Check and correct values less than zero or greater than pi */
		for(n=0;n<lsf->size1;n++) {
			for(i=0;i<lsf->size2;i++) {
				if(gsl_matrix_get(lsf,n,i) < 0) {
					gsl_matrix_set(lsf,n,i,LSF_EPSILON);
					flag_neg++;
					ok = 0;
				} else if(gsl_matrix_get(lsf,n,i) < LSF_EPSILON) {
					gsl_matrix_set(lsf,n,i,LSF_EPSILON);
					flag_zero++;
					ok = 0;
				} else if(gsl_matrix_get(lsf,n,i) > M_PI) {
					gsl_matrix_set(lsf,n,i,M_PI-LSF_EPSILON);
					flag_piplus++;
					ok = 0;
				} else if(gsl_matrix_get(lsf,n,i) > M_PI-LSF_EPSILON) {
					gsl_matrix_set(lsf,n,i,M_PI-LSF_EPSILON);
					flag_pi++;
					ok = 0;
				}
				if(gsl_isnan(gsl_matrix_get(lsf,n,i))) {
					if(i == 0)
						gsl_matrix_set(lsf,n,i,LSF_EPSILON);
					else if(i == lsf->size2-1)
						gsl_matrix_set(lsf,n,i,M_PI-LSF_EPSILON);
					else
						gsl_matrix_set(lsf,n,i,(gsl_matrix_get(lsf,n,i-1)+gsl_matrix_get(lsf,n,i+1))/2.0);
					flag_nan++;
					ok = 0;
				}
			}

			/* Check and correct non-increasing values or coefficients too close */
			for(i=0;i<lsf->size2-1;i++) {
				if(gsl_matrix_get(lsf,n,i) > gsl_matrix_get(lsf,n,i+1)) {
					mean = (gsl_matrix_get(lsf,n,i)+gsl_matrix_get(lsf,n,i+1))/2.0;
					gsl_matrix_set(lsf,n,i,mean - LSF_EPSILON/2.0);
					gsl_matrix_set(lsf,n,i+1,mean + LSF_EPSILON/2.0);
					flag_dec++;
					ok = 0;
				} else if(gsl_matrix_get(lsf,n,i) > gsl_matrix_get(lsf,n,i+1)-LSF_EPSILON/10.0) {
					mean = (gsl_matrix_get(lsf,n,i)+gsl_matrix_get(lsf,n,i+1))/2.0;
					gsl_matrix_set(lsf,n,i,mean - LSF_EPSILON/2.0);
					gsl_matrix_set(lsf,n,i+1,mean + LSF_EPSILON/2.0);
					flag_close++;
					ok = 0;
				}
			}
		}
	}

	/* Report */
	if(flag_dec > 0)
		printf("Warning: %i decreasing LSFs -> Fixed!\n",flag_dec);
	if(flag_neg > 0)
		printf("Warning: %i negative LSFs -> Fixed!\n",flag_neg);
	if(flag_piplus > 0)
		printf("Warning: %i LSFs greater than Pi -> Fixed!\n",flag_piplus);
	if(flag_zero > 0)
		printf("Warning: %i LSFs too close to 0 -> Fixed!\n",flag_zero);
	if(flag_pi > 0)
		printf("Warning: %i LSFs too close to Pi -> Fixed!\n",flag_pi);
	if(flag_close > 0)
		printf("Warning: %i LSFs too close to each other -> Fixed!\n",flag_close);
	if(flag_nan > 0)
		printf("Warning: %i NaN LSF value(s) -> Fixed!\n",flag_nan);
}





/**
 * Function Merge_voiced_unvoiced_spectra
 *
 * Merge voiced and unvoiced spectra, i.e., replace unvoiced spectrum of LSF with LSF2 if two spectra is used
 *
 * @param LSF
 * @param LSF2
 * @param params
 */
void Merge_voiced_unvoiced_spectra(gsl_matrix *LSF, gsl_matrix *LSF2, gsl_vector *fundf, PARAM *params) {

	int i,j;
	if(params->sep_vuv_spectrum == 1)
		for(i=0; i<params->n_frames; i++)
			if(gsl_vector_get(fundf, i) == 0)
				for(j=0; j<LSF->size2; j++)
					gsl_matrix_set(LSF, i, j, gsl_matrix_get(LSF2, i, j));
}




/**
 * Function Integrate_LSFs
 *
 * Integrate line spectral frequency (LSF) based parameters if differential LSFs are used.
 * First raise to the power of 2, integrate, scale according to the distance to PI.
 *
 * @param LSF
 * @param LSF2
 * @param params
 */
void Integrate_LSFs(gsl_matrix **LSFd, PARAM *params) {

	if(params->differential_lsf == 0)
		return;

	/* Initialize */
	int i,k;
	double len_orig,len_new,c;
	gsl_matrix *LSF = gsl_matrix_calloc((*LSFd)->size1,(*LSFd)->size2-1);

	/* Remove sqrt */
	for(i=0;i<LSF->size1;i++) {
		for(k=0;k<LSF->size2;k++)
			gsl_matrix_set(*LSFd,i,k,pow(gsl_matrix_get(*LSFd,i,k),2));

		/* Set initial LSF */
		gsl_matrix_set(LSF,i,0,gsl_matrix_get(*LSFd,i,0));
	}

	/* Integrate */
	for(i=0;i<LSF->size1;i++)
		for(k=1;k<LSF->size2;k++)
			gsl_matrix_set(LSF,i,k,gsl_matrix_get(*LSFd,i,k) + gsl_matrix_get(LSF,i,k-1));

	/* Scale LSFs according to the last coefficient so that the distance to PI matches */
	for(i=0;i<LSF->size1;i++) {
		len_orig = M_PI - pow(gsl_matrix_get(*LSFd,i,(*LSFd)->size2-1),2);
		len_new = gsl_matrix_get(LSF,i,LSF->size2-1);
		c = len_orig/len_new;
		for(k=0;k<LSF->size2;k++)
			gsl_matrix_set(LSF,i,k,gsl_matrix_get(LSF,i,k)*c);
	}

	/* Free old matrix */
	gsl_matrix_free(*(LSFd));

	/* Set pointer to new matrix */
	*(LSFd) = LSF;
}



/**
 * Function Noise_robust_speech1
 *
 * Modification for noise robust speech synthesis (1st stage)
 *
 * @param pulse_clus_id pointer to pulse cluster ids
 * @param pulse_clusters pointer to pulse clusters
 * @param params paramteres structure
 * @param synfilenumber
 * @return
 */
void Noise_robust_speech1(PARAM *params) {

	if(params->noise_robust_speech == 0)
		return;

	/* Modifications */
	params->hpfiltf0 = 1;
	params->postfilter_alpha = 0.15;        // 0<c<1, 1:OFF
	params->compcoeff = 0.2;                // 0<c<1, 1:OFF

	/* Pitch and speed */
	params->pitch = 1.4;
	params->speed = 0.8;
}



/**
 * Function Noise_robust_speech2
 *
 * Modification for noise robust speech synthesis (2nd stage)
 *
 * @param pulse_clus_id pointer to pulse cluster ids
 * @param pulse_clusters pointer to pulse clusters
 * @param params paramteres structure
 * @param synfilenumber
 * @return
 */
void Noise_robust_speech2(gsl_vector *gain, gsl_matrix *harmonics, PARAM *params) {

	if(params->noise_robust_speech == 0)
		return;

	/* Modify harmonics */
	if(params->use_harmonics == 1)
		Modify_harmonics(harmonics,0.5); 			// 0<c<1, 1:OFF
	else {
		params->pulse_tilt_decrease_coeff = 0.5; 	// 0<c<1, 1:OFF
		params->use_harmonic_modification = 1;
	}

	/* Compression of gain */
	if(params->use_pulselib == 1) {
		int i;
		double compression = 0.99; // The lower the number, the more compression will take effect
		double addition = 9;       // dB, rises the baseline gain
		double mingain = gsl_vector_min(gain);
		for(i=0;i<gain->size;i++)
			gsl_vector_set(gain,i,pow(gsl_vector_get(gain,i) - mingain, compression) + addition + mingain);
	}
}


















/**
 * Function CreateExcitation
 *
 * Creates excitation...
 *
 * @param ...
 *
 */
void CreateExcitation(PARAM *params,gsl_vector *excitation_voiced,gsl_vector *excitation_unvoiced,gsl_vector *fundf,gsl_vector *gain,
	      gsl_matrix *lsf,gsl_matrix *glflowsp,gsl_matrix *hnr,gsl_matrix *harmonics,gsl_matrix *waveform,gsl_vector *h1h2,
	      gsl_vector *naq,gsl_vector *original_pulse,gsl_matrix *pulses,gsl_matrix *pulses_rs,gsl_vector *pulse_lengths,
	      gsl_vector *pgain,gsl_matrix *plsf,gsl_matrix *ptilt,gsl_matrix *phnr,gsl_matrix *pharm,gsl_matrix *pwaveform,
	      gsl_vector *ph1h2, gsl_vector *pnaq,gsl_vector *resynthesis_pulse_index,gsl_vector *pulse_clus_id,
	      gsl_matrix *pulse_clusters,gsl_vector *oldgain,gsl_matrix *glflowsp_new,gsl_matrix *hnr_new,
	      gsl_vector *pca_mean, gsl_matrix *pca_pc, gsl_matrix *pca_w, gsl_matrix *pca_w_lib, gsl_vector *stoch_env, gsl_vector *stoch_sp,
	      gsl_matrix **DNN_W, gsl_vector **input_minmax, gsl_vector *dnnpulseindices, gsl_vector *dnnpulses) {

	/* Clear parameters (if resynth) */
	gsl_vector_set_zero(excitation_voiced);
	gsl_vector_set_zero(excitation_unvoiced);
	gsl_matrix_set_zero(glflowsp_new);

	/* Allocate memory */
	int sample_index = 0,frame_index_old = 0;
	int frame_index_center,frame_index,i,j,N,NN,pulseind = 0;
	double ngain_uv,gain_v,sum,sum_d,modg;
	int end_flag = 0;
	gsl_vector *noise;
	gsl_vector *pulse_train;
	gsl_vector *inds = gsl_vector_alloc(params->n_pulsecandidates);
	gsl_vector *inds_next = gsl_vector_alloc(params->n_pulsecandidates);
	gsl_vector *eparam = gsl_vector_alloc(params->n_pulsecandidates);
	gsl_vector *eparam_next = gsl_vector_alloc(params->n_pulsecandidates);
	gsl_vector *epulse = gsl_vector_alloc(params->n_pulsecandidates);
	gsl_vector *efinal = gsl_vector_alloc(params->n_pulsecandidates);
	gsl_vector *prevpulse = gsl_vector_alloc(params->rspulsemaxlen);

	/* Modulation */
	int mod_pulse_index = 0;
	double mod_pulse_gain = 0.0;

	/* DNN */
	int dnnpind = 0;
	int dnnpulseind = 0;


	/***********************************************************/
	/* Create excitation: loop until signal length is reached  */
	/***********************************************************/

	while(sample_index < params->signal_length) {

		/* Define current time index */
		frame_index = floor(params->n_frames*(sample_index/(double)params->signal_length));
		if(frame_index > params->n_frames-1) {frame_index = params->n_frames-1;}


		/****************************************************/
		/* 1. Segment is voiced                             */
		/****************************************************/

		if(gsl_vector_get(fundf,frame_index) != 0) {


			/****************************************************/
			/* Interpolate one pulse (original implementation)  */
			/****************************************************/

			if(params->use_pulselib == 0 && params->use_pulselib_pca == 0 && params->use_dnn_pulsegen == 0) {

				/* Interpolate the original pulse according to T0 (N), use cubic spline interpolation */
				N = rint(params->FS/gsl_vector_get(fundf,frame_index)/params->pitch);
				gsl_vector *pulse;
				if(params->two_pitch_period_diff_pulse == 0)
					pulse = gsl_vector_calloc(N);
				else
					pulse = gsl_vector_calloc(2*N);
				Interpolate(original_pulse,pulse);

				/* Create synthetic pulse train */
				if(params->two_pitch_period_diff_pulse == 0)
					pulse_train = Create_pulse_train(pulse,original_pulse,hnr,fundf,harmonics,frame_index,sample_index,params);
				else {
					pulse_train = Create_pulse_train_diff(pulse,original_pulse,hnr,fundf,harmonics,frame_index,sample_index,params);
					Integrate(pulse_train,LEAK);
				}

				/* Re-estimate HNR */
				if(params->hnr_reestimated == 0) {
					Upper_lower_envelope(pulse_train, hnr_new, gsl_vector_get(fundf,frame_index), frame_index, params);
					FillHNRValues(hnr_new,frame_index_old,frame_index);
					frame_index_old = frame_index;
					gsl_vector_free(pulse_train);
				} else { /* Add noise and analyse pulse train spectrum */
					if(params->two_pitch_period_diff_pulse == 1)
						Integrate(pulse,LEAK);
					Phase_manipulation(pulse,hnr,harmonics,frame_index,params);
					if(params->two_pitch_period_diff_pulse == 1) {
						Differentiate(pulse,LEAK);
						gsl_vector_set(pulse,0,0); // Fix DC error
					}
					Analyse_pulse_train_spectrum(pulse_train,glflowsp_new,N,sample_index,frame_index,params);
					frame_index_old = frame_index;
					gsl_vector_free(pulse_train);
				}

				/* Add jitter */
				int N_orig = N;
				if(params->jitter > 0) {
					N = N + rint(RAND()*N*params->jitter);
					gsl_vector *pulse_jitter;
					if(params->two_pitch_period_diff_pulse == 0)
						pulse_jitter = gsl_vector_alloc(N);
					else
						pulse_jitter = gsl_vector_alloc(2*N);
					Interpolate(pulse,pulse_jitter);
					gsl_vector_free(pulse);
					pulse = pulse_jitter;
				}

				/* Truncate pulse for unvoiced sections */
				if(params->two_pitch_period_diff_pulse == 0) {
					pulse = Truncate_pulse(pulse,fundf,sample_index,frame_index,params);
					N_orig = pulse->size;
					N = pulse->size;
				}

				/* Prevent going over signal length */
				if(params->two_pitch_period_diff_pulse == 0) {
					if(sample_index+N > params->signal_length) {
						gsl_vector_free(pulse);
						break;
					}
				}

				/* Normalize gain according to one pitch period */
				if(params->two_pitch_period_diff_pulse == 0) {
					gsl_vector *pulse_d = gsl_vector_alloc(pulse->size);
					gsl_vector_memcpy(pulse_d,pulse);
					LipRadiation(pulse_d);
					sum_d = 0;
					for(i=0;i<N;i++)
						sum_d = sum_d + gsl_vector_get(pulse_d,i)*gsl_vector_get(pulse_d,i);
					gsl_vector_free(pulse_d);
				} else {
					sum = 0;
					for(j=0;j<pulse->size;j++)
						sum = sum + gsl_vector_get(pulse,j)*gsl_vector_get(pulse,j);
				}

				/* Evaluate gain */
				if(params->two_pitch_period_diff_pulse == 0) {
					frame_index_center = GSL_MIN(floor(params->n_frames*((sample_index + 0.5*N_orig)/(double)params->signal_length)),params->n_frames-1);
					gain_v = sqrt(N*E_REF*powf(10.0,gsl_vector_get(gain,frame_index_center)/10.0)/sum_d);
				} else
					gain_v = sqrt(pulse->size*E_REF*powf(10.0,gsl_vector_get(gain,frame_index)/10.0)/(sum/0.375));

				/* Modulation */
				modg = 1 + sin(M_PI*mod_pulse_index/2.0)*mod_pulse_gain;
				gain_v = modg*gain_v;
				mod_pulse_index++;

				/* Set pulse to excitation */
				if(params->two_pitch_period_diff_pulse == 0) {
					for(i=0;i<N;i++)
						gsl_vector_set(excitation_voiced, sample_index+i, gain_v*(gsl_vector_get(pulse,i) - gsl_vector_get(pulse,0)));
				} else {
					for(j=0;j<pulse->size;j++) {
						gsl_vector_set(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(pulse->size/2),0),excitation_voiced->size-1),
						   gsl_vector_get(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(pulse->size/2),0),excitation_voiced->size-1)) +
						   gain_v*gsl_vector_get(pulse,j));
					}
				}

				/* Free memory */
				gsl_vector_free(pulse);

				/* Increment sample index */
				sample_index += N_orig;


			} else if(params->use_pulselib_pca == 0 && params->use_dnn_pulsegen == 0) {


				/*******************************************/
				/* Use voice source unit selection         */
				/*******************************************/

				/* Get the length of voiced section in frames and pulses */
				int cur_frame_index = frame_index;
				int cur_sample_index = sample_index;
				int n_pulses_in_section = 0;
				double mf0 = 0;
				while(1) {
					if(cur_frame_index > params->n_frames-1) {
						end_flag = 1;
						break;
					}
					n_pulses_in_section++;
					mf0 += gsl_vector_get(fundf,cur_frame_index);
					N = rint(params->FS/gsl_vector_get(fundf,cur_frame_index)/params->pitch);
					cur_frame_index = floor(params->n_frames*((cur_sample_index+N)/(double)params->signal_length));
					cur_sample_index += N;
					if(cur_frame_index > params->n_frames-1) {
						end_flag = 1;
						break;
					}
					if(gsl_vector_get(fundf,cur_frame_index) == 0 || gsl_vector_get(fundf,GSL_MIN(floor(params->n_frames*((cur_sample_index + rint(0.5*params->FS/gsl_vector_get(fundf,cur_frame_index)/params->pitch))/(double)params->signal_length)),params->n_frames-1)) == 0)
						break;
				}

				/* Evaluate mean f0 for the voiced section */
				mf0 = mf0/(double)n_pulses_in_section;

				/* Variables for Viterbi */
				gsl_matrix *v_trellis;
				gsl_matrix *v_indices;
				gsl_vector *final_pulses = gsl_vector_alloc(n_pulses_in_section);

				/* If pulse indices not known */
				if(params->resynth == 0) {

					/* Print number of pulses in voiced section */
					//printf("              Voiced section: %i pulse(s)\n",n_pulses_in_section);

					/* Allocate lattices for viterbi */
					v_trellis = gsl_matrix_alloc(n_pulses_in_section, params->n_pulsecandidates+1); // Onko oikein +1?
					v_indices = gsl_matrix_alloc(n_pulses_in_section, params->n_pulsecandidates+1);
					gsl_matrix_set_all(v_trellis, BIGGER_POS_NUMBER);
					gsl_matrix_set_all(v_indices, 1);

					/* Populate targets according to lowest target error */
					cur_frame_index = frame_index;
					cur_sample_index = sample_index;
					for(i=0;i<n_pulses_in_section;i++) {
						gsl_vector_view inds = gsl_matrix_row(v_indices,i);
						gsl_vector_view eparam = gsl_matrix_row(v_trellis,i);
						Evaluate_target_error(lsf,glflowsp,harmonics,hnr,waveform,h1h2,naq,oldgain,fundf,plsf,ptilt,pharm,phnr,pwaveform,ph1h2,pnaq,
								pgain,pca_w,pca_w_lib,pulse_lengths,cur_frame_index,&inds.vector,&eparam.vector,pulse_clus_id,pulse_clusters,params);
						N = rint(params->FS/gsl_vector_get(fundf,cur_frame_index)/params->pitch);
						cur_frame_index = floor(params->n_frames*((cur_sample_index+N)/(double)params->signal_length));
						cur_sample_index += N;
					}
					if(n_pulses_in_section > 1) {

						/* Viterbi, forward */
						/* Cumulative score */
						gsl_matrix *v_scores = gsl_matrix_alloc(n_pulses_in_section, params->n_pulsecandidates+1);
						gsl_matrix_set_all(v_scores, BIGGER_POS_NUMBER);

						/* Best node indices leading to each node */
						gsl_matrix *v_best = gsl_matrix_alloc(n_pulses_in_section, params->n_pulsecandidates+1);
						gsl_matrix_set_all(v_best, -1);

						/* Evaluate concatenation error */
						cur_frame_index = frame_index;
						cur_sample_index = sample_index;
						int dist_to_unvoiced;
						double additional_cost;
						for(i=0;i<n_pulses_in_section-1;i++) {
							gsl_vector_view inds = gsl_matrix_row(v_indices,i);
							gsl_vector_view inds_next = gsl_matrix_row(v_indices,i+1);
							gsl_vector_view eparam = gsl_matrix_row(v_trellis,i);
							gsl_vector_view eparam_next = gsl_matrix_row(v_trellis,i+1);

							/* Define distance to unvoiced */
							dist_to_unvoiced = i+1;

							/* Define additional cost (consisting of HNR and F0) */
							if(params->use_hnr == 1) {
								additional_cost = 0;
								for(j=0;j<hnr->size2;j++)
									additional_cost += 1.0-pow(10,gsl_matrix_get(hnr,cur_frame_index,j)/20.0);
								additional_cost = powf(additional_cost/hnr->size2,3.0);
							} else
								additional_cost = 1;
							additional_cost *= gsl_vector_get(fundf, frame_index)/120.0; // Add cost if high f0 (many concatenation points)

							/* Evaluate concatenation cost */
							Evaluate_concatenation_error_viterbi(&eparam.vector,&eparam_next.vector,&inds.vector,&inds_next.vector,
									pulses_rs, v_scores, v_best, i, dist_to_unvoiced, additional_cost, params);
							N = rint(params->FS/gsl_vector_get(fundf,cur_frame_index)/params->pitch);
							cur_frame_index = floor(params->n_frames*((cur_sample_index+N)/(double)params->signal_length));
							cur_sample_index += N;
						}

						/* Viterbi, get best path by backtracking */
						double best_sc = BIGGER_POS_NUMBER;
						int best_i = -1;
						for(i=0;i<params->n_pulsecandidates;i++) {
							if(gsl_matrix_get(v_scores,n_pulses_in_section-1, i) < best_sc) {
								best_sc = gsl_matrix_get(v_scores, n_pulses_in_section-1,i);
								best_i = i;
							}
						}
						gsl_vector_set(final_pulses, n_pulses_in_section-1, gsl_matrix_get(v_indices, n_pulses_in_section-1, best_i));
						for(i=n_pulses_in_section-2;i>= 0;i--) {
							gsl_vector_set(final_pulses, i, gsl_matrix_get(v_indices, i+1, best_i));
							best_i = (int)gsl_matrix_get(v_best, i+1,best_i);
						}

						/* Free viterbi structures */
						gsl_matrix_free(v_best);
						gsl_matrix_free(v_scores);

					} else {

						/* Only one pulse in section */
						gsl_vector_set(final_pulses,0,gsl_matrix_get(v_indices,0,i)); // TODO: Should be best_i(i) ?
					}

					/* Free viterbi structures */
					gsl_matrix_free(v_trellis);
					gsl_matrix_free(v_indices);


				} /* After this, same for resynthesis */

				/* Make pulse train */
				double prev_f0 = 0.0;
				for(i=0;i<=n_pulses_in_section-1;i++) {

					/* Define pulse index */
					if(params->resynth == 0) {
						pulseind = (int)gsl_vector_get(final_pulses, i);
						gsl_vector_set(resynthesis_pulse_index,frame_index,gsl_vector_get(final_pulses,i));
					} else {
						pulseind = gsl_vector_get(resynthesis_pulse_index, frame_index);
					}

					/* Print pulse index */
					//if(params->resynth == 0)
					//	printf("%i\n",pulseind);

					/* Get pulse length and f0 shift */
					double f0 = gsl_vector_get(fundf,frame_index);
					if(f0 == 0) f0 = prev_f0;
					if(f0 == 0) f0 = mf0;
					N = rint(params->FS/f0/params->pitch);
					NN = gsl_vector_get(pulse_lengths,pulseind);

					/* Get pulse to vector */
					gsl_vector *pulse;
					if(params->pulse_interpolation == 0) {
						if(sample_index+rint(NN/2) > params->signal_length)
							break;
						pulse = gsl_vector_alloc(NN);
						for(j=0;j<NN;j++)
							gsl_vector_set(pulse,j,gsl_matrix_get(pulses,pulseind,j));
					} else {
						if(sample_index+N > params->signal_length)
							break;
						gsl_vector *pulse_orig = gsl_vector_calloc(NN);
						pulse = gsl_vector_calloc(2*N);
						for(j=0;j<NN;j++)
							gsl_vector_set(pulse_orig,j,gsl_matrix_get(pulses,pulseind,j));
						Interpolate(pulse_orig,pulse);
						gsl_vector_free(pulse_orig);
						NN = 2*N;
					}

					/* Average pulses */
					Average_pulses(pulse, fundf, resynthesis_pulse_index, pulses, pulse_lengths, frame_index, params);

					/* Add jitter */
					if(params->jitter > 0) {
						NN = pulse->size + rint(RAND()*2*N*params->jitter);
						gsl_vector *pulse_jitter = gsl_vector_alloc(NN);
						Interpolate(pulse,pulse_jitter);
						gsl_vector_free(pulse);
						pulse = pulse_jitter;
					}

					/* Phase manipulation of the library pulse */
					if(params->add_noise_pulselib == 1 && params->resynth == 1) {
						Integrate(pulse, LEAK);
						Phase_manipulation(pulse,hnr,harmonics,frame_index,params);
						Differentiate(pulse,LEAK);
						gsl_vector_set(pulse,0,0); // Fix DC error
					}

					/* Evaluate gain */
					sum = 0;
					for(j=0;j<NN;j++)
						sum = sum + gsl_vector_get(pulse,j)*gsl_vector_get(pulse,j);
					gain_v = sqrt(NN*E_REF*powf(10.0,gsl_vector_get(gain,frame_index)/10.0)/(sum/0.375));

					/* Modulation */
					modg = 1 + sin(M_PI*mod_pulse_index/2.0)*mod_pulse_gain;
					gain_v = modg*gain_v;
					mod_pulse_index++;

					/* Set pulse to excitation (OLA) */
					for(j=0;j<NN;j++) {
						gsl_vector_set(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(NN/2),0),excitation_voiced->size-1),
						   gsl_vector_get(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(NN/2),0),excitation_voiced->size-1)) +
						   gain_v*gsl_vector_get(pulse,j));
					}
					gsl_vector_free(pulse);

					/* Go directly to the beginning of unvoiced segent if pulse exceeds the boundary */
					int frame_index_tmp = frame_index;
					int next_unvoiced_sample = sample_index;
					while(gsl_vector_get(fundf,frame_index_tmp) > 0) {
						frame_index_tmp = floor(params->n_frames*((next_unvoiced_sample)/(double)params->signal_length));
						next_unvoiced_sample++;
					}
					N = GSL_MIN(N, next_unvoiced_sample-sample_index);

					/* Previous F0, increment sample index */
					prev_f0 = gsl_vector_get(fundf, frame_index);
					sample_index += N;
					frame_index = floor(params->n_frames*((sample_index)/(double)params->signal_length));
				}

				/* Free memory */
				gsl_vector_free(final_pulses);

				/* Stop in the end */
				if(end_flag == 1)
					break;
	    		}


			/****************************************/
			/* Construct pulses from PCA parameters */
			/****************************************/

			else if(params->use_pulselib_pca == 1 && params->use_dnn_pulsegen == 0) {

	    			/* Allocate PCA pulse */
	    			gsl_vector *pca_pulse = gsl_vector_calloc(params->pca_pulse_length);

	    			/* Define pulse from PC and weights */
	    			if(params->pca_order_synthesis < 0 || params->pca_order_synthesis > params->pca_order) {
	    				printf("PCA_ORDER_SYNTHESIS must be between 0 and PCA_ORDER. PCA_ORDER_SYNTHESIS set to PCA_ORDER.\n");
	    				params->pca_order_synthesis = params->pca_order;
	    			}
	    			for(i=0;i<params->pca_order_synthesis;i++)
	    				for(j=0;j<params->pca_pulse_length;j++)
	    					gsl_vector_set(pca_pulse,j,gsl_vector_get(pca_pulse,j) + gsl_matrix_get(pca_w,frame_index,i)*gsl_matrix_get(pca_pc,j,i));

	    			/* Add mean */
	    			for(i=0;i<params->pca_pulse_length;i++)
	    				gsl_vector_set(pca_pulse,i,gsl_vector_get(pca_pulse,i) + gsl_vector_get(pca_mean,i));

	    			/* Define pulse length, add jitter */
	    			N = rint(params->FS/gsl_vector_get(fundf,frame_index)/params->pitch);
				NN = 2*N;
				if(params->jitter > 0)
					NN = rint(NN*(1.0 + RAND()*params->jitter));
				if(sample_index+NN/2.0 > params->signal_length) {
					gsl_vector_free(pca_pulse);
					break;
				}

				/* Allocate and interpolate pulse */
				gsl_vector *pulse = gsl_vector_alloc(NN);
				Interpolate(pca_pulse,pulse);

				/* Add noise */
				Integrate(pulse,LEAK);
				Phase_manipulation(pulse,hnr,harmonics,frame_index,params);
				Differentiate(pulse,LEAK);
				gsl_vector_set(pulse,0,0); // Fix DC error

				/* Create synthetic pulse train, analyse pulse train spectrum */
				pulse_train = Create_pulse_train_diff(pulse,pca_pulse,hnr,fundf,harmonics,frame_index,sample_index,params);
				Integrate(pulse_train,LEAK);
				Analyse_pulse_train_spectrum(pulse_train,glflowsp_new,N,sample_index,frame_index,params);
				frame_index_old = frame_index;
				gsl_vector_free(pulse_train);

				/* Evaluate gain */
				sum = 0;
				for(j=0;j<NN;j++)
					sum = sum + gsl_vector_get(pulse,j)*gsl_vector_get(pulse,j);
				gain_v = sqrt(NN*E_REF*powf(10.0,gsl_vector_get(gain,frame_index)/10.0)/(sum/0.375));

				/* Modulation (NOT USED) */
				modg = 1;
				gain_v = modg*gain_v;
				mod_pulse_index++;

				/* Set pulse to excitation (OLA) */
				for(j=0;j<NN;j++) {
					gsl_vector_set(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(NN/2),0),excitation_voiced->size-1),
					   gsl_vector_get(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(NN/2),0),excitation_voiced->size-1)) +
					   gain_v*gsl_vector_get(pulse,j));
				}

				/* Free memory */
				gsl_vector_free(pulse);
				gsl_vector_free(pca_pulse);

				/* Increment sample index */
				sample_index += N;
			}
			

			/****************************************/
			/* DNN pulse generation                 */
			/****************************************/

			// LOWER FOR DNN-PCA (switch commenting)
			else if(params->use_pulselib_pca == 0 && params->use_dnn_pulsegen == 1) {
			//else if(params->use_pulselib_pca == 1 && params->use_dnn_pulsegen == 1) {

				/* Generate pulse from DNNs */
				N = rint(params->FS/gsl_vector_get(fundf,frame_index)/params->pitch);
				gsl_vector *pulse = gsl_vector_calloc(2*N);

				/* If first synthesis, generate DNN pulse and search for closest library pulse (index) */
				if(params->resynth == 0) {

					/* Generate DNN pulse */
					/* FOR DNN-PCA, switch comments below and in Generate_excitation, set one switch as well */
					Generate_DNN_pulse(pulse,fundf,gain,naq,h1h2,hnr,glflowsp,lsf,DNN_W,input_minmax,frame_index,params);
					//Generate_DNN_pulse_PCA(pulse,fundf,gain,naq,h1h2,hnr,glflowsp,lsf,DNN_W,pca_pc,frame_index,params);

					/* Select between using DNN pulse itself, or selecting a natural pulse from library */
					if(params->use_dnn_pulselib_sel == 0) {

						/* Use DNN pulse and save it to memory */
						for(i=0;i<pulse->size;i++)
							gsl_vector_set(dnnpulses,dnnpulseind+i,gsl_vector_get(pulse,i));
						dnnpulseind += pulse->size;
						dnnpind++;

					} else {

						/* Select pulse from library, resample DNN pulse and normalize energy */
						gsl_vector *dnn_pulse_rs = gsl_vector_alloc(pulses_rs->size2);
						Interpolate(pulse,dnn_pulse_rs);
						double e = 0;
						for(i=0;i<dnn_pulse_rs->size;i++)
							e += gsl_vector_get(dnn_pulse_rs,i)*gsl_vector_get(dnn_pulse_rs,i);
						e = sqrt(e);
						for(i=0;i<dnn_pulse_rs->size;i++)
							gsl_vector_set(dnn_pulse_rs,i,gsl_vector_get(dnn_pulse_rs,i)/e);
			
						/* Search for a best matching pulse from pulse library */
						gsl_vector *error = gsl_vector_calloc(pulses_rs->size1);
						for(i=0;i<error->size;i++)
							for(j=0;j<dnn_pulse_rs->size;j++)
								gsl_vector_set(error,i,gsl_vector_get(error,i) + powf(gsl_matrix_get(pulses_rs,i,j)
									 - gsl_vector_get(dnn_pulse_rs,j),2));

 						/* Save pulse index */
						gsl_vector_set(dnnpulseindices,dnnpind,gsl_vector_min_index(error));

						/* Free memory */
						gsl_vector_free(error);
						gsl_vector_free(dnn_pulse_rs);
					}
				} else {

					/* Load DNN pulse from memory in resynthesis */
					if(params->use_dnn_pulselib_sel == 0) {

						/* Load DNN pulse from memory */
						for(i=0;i<pulse->size;i++)
							gsl_vector_set(pulse,i,gsl_vector_get(dnnpulses,dnnpulseind+i));
						dnnpulseind += pulse->size;
						dnnpind++;
					}
				}

				/* If pulse selection is used, get pulse from library according to index */
				if(params->use_dnn_pulselib_sel == 1) {	
					gsl_vector *libpulse_orig = gsl_vector_alloc(gsl_vector_get(pulse_lengths,gsl_vector_get(dnnpulseindices,dnnpind)));
					for(i=0;i<libpulse_orig->size;i++)
						gsl_vector_set(libpulse_orig,i,gsl_matrix_get(pulses,gsl_vector_get(dnnpulseindices,dnnpind),i));
					Interpolate(libpulse_orig,pulse);
					//if(params->resynth == 0)
					//	printf("%i\n",(int)gsl_vector_get(dnnpulseindices,dnnpind));
					gsl_vector_free(libpulse_orig);
					dnnpind++;
				}

				/* Evaluate voice source spectrum */
				gsl_vector *tmp_pulse = gsl_vector_alloc(pulse->size);
				gsl_vector_memcpy(tmp_pulse,pulse);
				pulse_train = Create_pulse_train_diff(pulse,tmp_pulse,hnr,fundf,harmonics,frame_index,sample_index,params);
				Integrate(pulse_train,LEAK);
				gsl_vector_free(tmp_pulse);

				/* Re-estimate HNR */
				if(params->hnr_reestimated == 0) {
					Upper_lower_envelope(pulse_train, hnr_new, gsl_vector_get(fundf,frame_index), frame_index, params);
					FillHNRValues(hnr_new,frame_index_old,frame_index);
					frame_index_old = frame_index;
					gsl_vector_free(pulse_train);
				} else { /* Add noise and analyse pulse train spectrum */
					Integrate(pulse,LEAK);
					Phase_manipulation(pulse,hnr,harmonics,frame_index,params);
					Differentiate(pulse,LEAK);
					gsl_vector_set(pulse,0,0); // Fix DC error
					Analyse_pulse_train_spectrum(pulse_train,glflowsp_new,N,sample_index,frame_index,params);
					frame_index_old = frame_index;
					gsl_vector_free(pulse_train);
				}

				/* Evaluate gain */
				sum = 0;
				for(j=0;j<pulse->size;j++)
					sum = sum + gsl_vector_get(pulse,j)*gsl_vector_get(pulse,j);
				gain_v = sqrt(pulse->size*E_REF*powf(10.0,gsl_vector_get(gain,frame_index)/10.0)/(sum/0.375));

				/* Set pulse to excitation (OLA) */
				for(j=0;j<pulse->size;j++) {
					gsl_vector_set(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(pulse->size/2),0),excitation_voiced->size-1),
					   gsl_vector_get(excitation_voiced,GSL_MIN(GSL_MAX(sample_index+j-rint(pulse->size/2),0),excitation_voiced->size-1)) +
					   gain_v*gsl_vector_get(pulse,j));
				}

				/* Free memory */
				gsl_vector_free(pulse);

				/* Increment sample index */
				sample_index += N;

			} // End single pulse / source unit selection / PCA pulse / DNN pulsegen

		} else {

			/****************************************************/
			/* 2. Segment is unvoiced                           */
			/****************************************************/


	    	/* Allocate noise vector */
			if(sample_index + params->shift/params->speed < params->signal_length)
				noise = gsl_vector_alloc(params->shift/params->speed);
			else
				noise = gsl_vector_alloc(params->signal_length-sample_index);

			/* Create noise excitation */
			double sum = 0;
			for(i=0;i<noise->size;i++) {
				gsl_vector_set(noise,i,RAND());
				sum += gsl_vector_get(noise,i)*gsl_vector_get(noise,i);
			}

			/* Add noise to voiced excitation */
			frame_index_center = GSL_MIN(floor(params->n_frames*((sample_index + 0.5*params->shift/params->speed)/(double)params->signal_length)),params->n_frames-1);
			ngain_uv = params->gain_unvoiced*sqrt(params->shift/params->speed*E_REF*powf(10.0,gsl_vector_get(gain,frame_index_center)/10.0)/sum);
			for(i=0;i<noise->size;i++)
				gsl_vector_set(excitation_unvoiced, sample_index+i, ngain_uv*gsl_vector_get(noise,i));
			gsl_vector_free(noise);

			/* Set synthetic source spectrum the same as original when unvoiced (FLAT FOR TESTING) */
			FillUnvoicedSyntheticSourceSpectrum(glflowsp,glflowsp_new,frame_index,frame_index_old);
			//FillUnvoicedSyntheticSourceSpectrum_FLAT(glflowsp,glflowsp_new,frame_index,frame_index_old);

			/* Update, increment sample index */
			frame_index_old = frame_index;
			sample_index += rint(params->shift/params->speed);
		}
	}

	/* Free memory */
	gsl_vector_free(prevpulse);
	gsl_vector_free(inds);
	gsl_vector_free(inds_next);
	gsl_vector_free(eparam);
	gsl_vector_free(eparam_next);
	gsl_vector_free(efinal);
	gsl_vector_free(epulse);
}







/**
 * Generate_DNN_pulse
 *
 * Generate glottal flow pulse from DNN through mapping of the synthesis parameters to the pulse waveform.
 *
 */
void Generate_DNN_pulse(gsl_vector *pulse, gsl_vector *fundf, gsl_vector *gain, gsl_vector *naq, gsl_vector *h1h2,
	gsl_matrix *hnr, gsl_matrix *glflowsp, gsl_matrix *lsf, gsl_matrix **DNN_W, gsl_vector **input_minmax, int index, PARAM *params) {

	/* Initialize */
	int i,j,L,numlayers = params->dnn_weight_dims->size/2;

	/* Allocate memory */
	gsl_vector *inputdata = gsl_vector_calloc(DNN_W[0]->size1);
	gsl_vector **WProbs = (gsl_vector**)malloc((numlayers-1)*sizeof(gsl_vector*));
	for(i=0;i<numlayers-1;i++)
		WProbs[i] = gsl_vector_calloc(DNN_W[i]->size2+1);
	gsl_vector *pulse_tmp = gsl_vector_calloc(DNN_W[numlayers-1]->size2);

	/* Assign input data to vector [F0 Gain HNR LSFsource LSF (Bias)] */
	/* Dimensions = [1 1 5 10 30 (1)] = 47 (+1) */
	int NPAR = 47;
	int f,cur_index,stack = 0;
	int lim1 = -floor(params->dnn_number_of_stacked_frames/2.0);
	int lim2 = ceil(params->dnn_number_of_stacked_frames/2.0);
	for(f=lim1;f<lim2;f++) {

		/* Evaluate current frame */
		cur_index = GSL_MIN(GSL_MAX(index+f,0),params->n_frames-1);

		/* Assign parameters */
		gsl_vector_set(inputdata,0+stack*NPAR,gsl_vector_get(fundf,cur_index));   			// F0
		gsl_vector_set(inputdata,1+stack*NPAR,gsl_vector_get(gain,cur_index));    	 		// Gain
		for(i=0;i<5;i++)
			gsl_vector_set(inputdata,2+stack*NPAR+i,gsl_matrix_get(hnr,cur_index,i));      	// HNR
		for(i=0;i<10;i++)
			gsl_vector_set(inputdata,7+stack*NPAR+i,gsl_matrix_get(glflowsp,cur_index,i)); 	// LSFsource
		for(i=0;i<30;i++)
			gsl_vector_set(inputdata,17+stack*NPAR+i,gsl_matrix_get(lsf,cur_index,i));     	// LSF

		/* Normalize input data to values between [0.1,0.9] (if done so in DNN training) */
		/* data_norm = 0.1 + 0.8*(data - min)/(max - min); */
		if(params->dnn_input_normalized == 1) {
			double min,max;
			for(i=0;i<NPAR;i++) {
				min = gsl_vector_get(input_minmax[0],i);
				max = gsl_vector_get(input_minmax[0],i+NPAR);
				gsl_vector_set(inputdata,i+stack*NPAR, 0.1 + 0.8*(gsl_vector_get(inputdata,i+stack*NPAR) - min)/(max-min));
			}
		}

		/* Increment stack */
		stack++;
	}

	/* Set the bias at the end of the vector */
	gsl_vector_set(inputdata,inputdata->size-1,1);

	/* Evaluate WProbs[0] = 1./(1 + exp(-data*DNN_W[0])) (sigmoid layer) */
	for(i=0;i<DNN_W[0]->size2;i++) {
		for(j=0;j<DNN_W[0]->size1;j++)
			gsl_vector_set(WProbs[0],i,gsl_vector_get(WProbs[0],i) + gsl_vector_get(inputdata,j)*gsl_matrix_get(DNN_W[0],j,i));
		gsl_vector_set(WProbs[0],i,1.0/(1.0 + exp(-gsl_vector_get(WProbs[0],i))));
	}
	gsl_vector_set(WProbs[0],WProbs[0]->size-1,1); // Set bias term

	/* Evaluate WPprobs[i] = 1./(1 + exp(-WProbs[i-1]*DNN_W[i])) (sigmoid layer) */
	for(L=1;L<numlayers-1;L++) {
		for(i=0;i<DNN_W[L]->size2;i++) {
			for(j=0;j<DNN_W[L]->size1;j++)
				gsl_vector_set(WProbs[L],i,gsl_vector_get(WProbs[L],i) + gsl_vector_get(WProbs[L-1],j)*gsl_matrix_get(DNN_W[L],j,i));
			gsl_vector_set(WProbs[L],i,1.0/(1.0 + exp(-gsl_vector_get(WProbs[L],i))));
		}
		gsl_vector_set(WProbs[L],WProbs[L]->size-1,1); // Set bias term
	}

	/* Evaluate pulses = WProbs[end]*DNN_W[end] (linear layer) */
	for(i=0;i<DNN_W[numlayers-1]->size2;i++)
		for(j=0;j<DNN_W[numlayers-1]->size1;j++)
			gsl_vector_set(pulse_tmp,i,gsl_vector_get(pulse_tmp,i) + gsl_vector_get(WProbs[numlayers-2],j)*gsl_matrix_get(DNN_W[numlayers-1],j,i));

	/* Interpolate pulse to desired length */
	Interpolate(pulse_tmp,pulse);

	/* Free memory */
	gsl_vector_free(inputdata);
	for(i=0;i<numlayers-1;i++)
		gsl_vector_free(WProbs[i]);
	free(WProbs);
	gsl_vector_free(pulse_tmp);
}







/**
 * Generate_DNN_pulse_PCA
 *
 * Generate glottal flow pulse by mapping the synthesis parameters to the PC weights and then reconstructing the pulse through PC weights and PCs.
 *
 */
void Generate_DNN_pulse_PCA(gsl_vector *pulse, gsl_vector *fundf, gsl_vector *gain, gsl_vector *naq, gsl_vector *h1h2,
	gsl_matrix *hnr, gsl_matrix *glflowsp, gsl_matrix *lsf, gsl_matrix **DNN_W, gsl_vector **input_minmax, gsl_matrix *pca_pc, int index, PARAM *params) {

	/* Initialize */
	int i,j,L,numlayers = params->dnn_weight_dims->size/2;

	/* Allocate memory */
	gsl_vector *inputdata = gsl_vector_calloc(DNN_W[0]->size1);
	gsl_vector **WProbs = (gsl_vector**)malloc((numlayers-1)*sizeof(gsl_vector*));
	for(i=0;i<numlayers-1;i++)
		WProbs[i] = gsl_vector_calloc(DNN_W[i]->size2+1);
	gsl_vector *pcw = gsl_vector_calloc(DNN_W[numlayers-1]->size2);

	/***** NORMAL *****/
	/* Assign input data to vector [F0 Gain HNR LSFsource LSF Bias] */
	/* Dimensions = [1 1 5 10 30 1] = 48 */
	gsl_vector_set(inputdata,0,gsl_vector_get(fundf,index));    // F0
	gsl_vector_set(inputdata,1,gsl_vector_get(gain,index));     // Gain
	for(i=0;i<5;i++)
		gsl_vector_set(inputdata,2+i,gsl_matrix_get(hnr,index,i));      // HNR
	for(i=0;i<10;i++)
		gsl_vector_set(inputdata,7+i,gsl_matrix_get(glflowsp,index,i)); // LSFsource
	for(i=0;i<30;i++)
		gsl_vector_set(inputdata,17+i,gsl_matrix_get(lsf,index,i));     // LSF
	gsl_vector_set(inputdata,47,1);                                         // Bias

	/***** INCLUDING NAQ AND H1H2 *****/
	/* Assign input data to vector [F0 Gain NAQ H1H2 HNR LSFsource LSF Bias] */
	/* Dimensions = [1 1 1 1 5 10 30 1] = 50 
	gsl_vector_set(inputdata,0,gsl_vector_get(fundf,index));    // F0
	gsl_vector_set(inputdata,1,gsl_vector_get(gain,index));     // Gain
	gsl_vector_set(inputdata,2,gsl_vector_get(naq,index));      // NAQ
	gsl_vector_set(inputdata,3,gsl_vector_get(h1h2,index));     // H1H2
	for(i=0;i<5;i++)
		gsl_vector_set(inputdata,4+i,gsl_matrix_get(hnr,index,i));      // HNR
	for(i=0;i<10;i++)
		gsl_vector_set(inputdata,9+i,gsl_matrix_get(glflowsp,index,i)); // LSFsource
	for(i=0;i<30;i++)
		gsl_vector_set(inputdata,19+i,gsl_matrix_get(lsf,index,i));     // LSF
	gsl_vector_set(inputdata,49,1);                                         // Bias
	*/

	/* Evaluate WProbs[0] = 1./(1 + exp(-data*DNN_W[0])) (sigmoid layer) */
	for(i=0;i<DNN_W[0]->size2;i++) {
		for(j=0;j<DNN_W[0]->size1;j++)
			gsl_vector_set(WProbs[0],i,gsl_vector_get(WProbs[0],i) + gsl_vector_get(inputdata,j)*gsl_matrix_get(DNN_W[0],j,i));
		gsl_vector_set(WProbs[0],i,1.0/(1.0 + exp(-gsl_vector_get(WProbs[0],i))));
	}
	gsl_vector_set(WProbs[0],WProbs[0]->size-1,1); // Set bias term

	/* Evaluate WPprobs[i] = 1./(1 + exp(-WProbs[i-1]*DNN_W[i])) (sigmoid layer) */
	for(L=1;L<numlayers-1;L++) {
		for(i=0;i<DNN_W[L]->size2;i++) {
			for(j=0;j<DNN_W[L]->size1;j++)
				gsl_vector_set(WProbs[L],i,gsl_vector_get(WProbs[L],i) + gsl_vector_get(WProbs[L-1],j)*gsl_matrix_get(DNN_W[L],j,i));
			gsl_vector_set(WProbs[L],i,1.0/(1.0 + exp(-gsl_vector_get(WProbs[L],i))));
		}
		gsl_vector_set(WProbs[L],WProbs[L]->size-1,1); // Set bias term
	}

	/* Evaluate pulses = WProbs[end]*DNN_W[end] (linear layer) */
	for(i=0;i<DNN_W[numlayers-1]->size2;i++)
		for(j=0;j<DNN_W[numlayers-1]->size1;j++)
			gsl_vector_set(pcw,i,gsl_vector_get(pcw,i) + gsl_vector_get(WProbs[numlayers-2],j)*gsl_matrix_get(DNN_W[numlayers-1],j,i));

	/* Allocate PCA pulse */
	gsl_vector *pca_pulse = gsl_vector_calloc(params->pca_pulse_length);

	/* Define pulse from PC and weights */
	if(params->pca_order_synthesis < 0 || params->pca_order_synthesis > params->pca_order) {
		printf("PCA_ORDER_SYNTHESIS must be between 0 and PCA_ORDER. PCA_ORDER_SYNTHESIS set to PCA_ORDER.\n");
		params->pca_order_synthesis = params->pca_order;
	}
	for(i=0;i<params->pca_order_synthesis;i++)
		for(j=0;j<params->pca_pulse_length;j++)
			gsl_vector_set(pca_pulse,j,gsl_vector_get(pca_pulse,j) + gsl_vector_get(pcw,i)*gsl_matrix_get(pca_pc,j,i));

	/* Interpolate pulse to desired length */
	Interpolate(pca_pulse,pulse);

	/* Free memory */
	gsl_vector_free(inputdata);
	for(i=0;i<numlayers-1;i++)
		gsl_vector_free(WProbs[i]);
	free(WProbs);
	gsl_vector_free(pcw);
	gsl_vector_free(pca_pulse);
}






/**
 * Shift_right_and_add
 *
 * Shift all the elements of the vector one index to the right, and add the new value in the beginning.
 *
 * @param vector
 * @param value
 */
void Shift_right_and_add(gsl_vector *vector, double value) {

	int i;
	for(i=vector->size-1;i>0;i--)
		gsl_vector_set(vector,i,gsl_vector_get(vector,i-1));
	gsl_vector_set(vector,0,value);
}




/**
 * Truncate_pulse
 *
 * Truncate pulse for unvoiced regions
 *
 * @param pulse
 * @param fundf
 * @param sample_index
 * @param frame_index
 * @param params
 */
gsl_vector *Truncate_pulse(gsl_vector *pulse, gsl_vector *fundf, int sample_index, int frame_index, PARAM *params) {

	/* Initialize */
	int sample_index_tmp,frame_index_tmp,i = 1,edit_flag = 0;
	int uvshift = rint(params->shift/params->speed);

	/* Find unvoiced */
	while(i*uvshift < pulse->size) {
		sample_index_tmp = sample_index + i*uvshift;
		frame_index_tmp = floor(params->n_frames*(sample_index_tmp/(double)params->signal_length));
		if(frame_index_tmp > fundf->size-1)
			break;
		if(gsl_vector_get(fundf,frame_index_tmp) == 0) {
			edit_flag = 1;
			break;
		}
		i = i + 1;
	}

	/* Truncate */
	if(edit_flag == 1) {
		gsl_vector *pulse_short = gsl_vector_alloc(i*uvshift);
		for(i=0;i<pulse_short->size;i++)
			gsl_vector_set(pulse_short,i,gsl_vector_get(pulse,i)*HANN(pulse_short->size+i,2*pulse_short->size));
		gsl_vector_free(pulse);
		pulse = pulse_short;
	}

	/* Return */
	return pulse;
}





/**
 * Average_pulses
 *
 * Average adjacent pulses in time and quality
 *
 * @param ...
 */
void Average_pulses(gsl_vector *pulse, gsl_vector *fundf, gsl_vector *resynthesis_pulse_index, gsl_matrix *pulses, gsl_vector *pulse_lengths, int frame_index, PARAM *params) {

	/* Perform only in resynthesis */
	if(params->resynth == 0)
		return;

	/* Sanity check, the number of pulses must be odd */
	if(params->average_n_adjacent_pulses < 2)
		params->average_n_adjacent_pulses = 1;
	if(params->average_n_adjacent_pulses%2 == 0)
		params->average_n_adjacent_pulses++;

	/* Initialize */
	int j,k,N = pulse->size;
	gsl_matrix *adjacent_pulses = gsl_matrix_alloc(params->average_n_adjacent_pulses,N);
	gsl_vector *pulse_N = gsl_vector_alloc(N);
	gsl_vector *pulse_temp;

	/* Fill pulse matrix with the original pulse */
	Interpolate(pulse,pulse_N);
	for(j=0;j<params->average_n_adjacent_pulses;j++)
		for(k=0;k<N;k++)
			gsl_matrix_set(adjacent_pulses,j,k,gsl_vector_get(pulse_N,k));

	/* Fill possible earlier pulses */
	j = 1;
	int added_pulses = 0;
	while(frame_index-j >= 0 && gsl_vector_get(fundf,frame_index-j) > 0 && added_pulses != floor(params->average_n_adjacent_pulses/2)) {
		if(gsl_vector_get(resynthesis_pulse_index,frame_index-j) != 0) {
			int pulseind_back = gsl_vector_get(resynthesis_pulse_index,frame_index-j);
			int pulselen = gsl_vector_get(pulse_lengths,pulseind_back);
			pulse_temp = gsl_vector_alloc(pulselen);
			for(k=0;k<pulselen;k++)
				gsl_vector_set(pulse_temp,k,gsl_matrix_get(pulses,pulseind_back,k));
			Interpolate(pulse_temp,pulse_N);
			gsl_vector_free(pulse_temp);
			for(k=0;k<N;k++)
				gsl_matrix_set(adjacent_pulses,floor(params->average_n_adjacent_pulses/2)-added_pulses-1,k,gsl_vector_get(pulse_N,k));
			added_pulses++;
		}
		j++;
	}

	/* Fill possible latter pulses */
	j = 1;
	added_pulses = 0;
	while(frame_index+j < params->n_frames && gsl_vector_get(fundf,frame_index+j) > 0 && added_pulses != floor(params->average_n_adjacent_pulses/2)) {
		if(gsl_vector_get(resynthesis_pulse_index,frame_index+j) != 0) {
			int pulseind_forward = gsl_vector_get(resynthesis_pulse_index,frame_index+j);
			int pulselen = gsl_vector_get(pulse_lengths,pulseind_forward);
			pulse_temp = gsl_vector_alloc(pulselen);
			for(k=0;k<pulselen;k++)
				gsl_vector_set(pulse_temp,k,gsl_matrix_get(pulses,pulseind_forward,k));
			Interpolate(pulse_temp,pulse_N);
			gsl_vector_free(pulse_temp);
			for(k=0;k<N;k++)
				gsl_matrix_set(adjacent_pulses,floor(params->average_n_adjacent_pulses/2)+added_pulses+1,k,gsl_vector_get(pulse_N,k));
			added_pulses++;
		}
		j++;
	}

	/* Average adjacent pulses to pulse */
	gsl_vector_set_zero(pulse);
	for(k=0;k<N;k++)
		for(j=0;j<params->average_n_adjacent_pulses;j++)
			gsl_vector_set(pulse,k,gsl_vector_get(pulse,k) + gsl_matrix_get(adjacent_pulses,j,k)/params->average_n_adjacent_pulses);

	/* Free memory */
	gsl_vector_free(pulse_N);
	gsl_matrix_free(adjacent_pulses);
}













/**
 * Function Fill_pulse_indices
 *
 * Fill pulse indices with 0 values with nearest pulsi indices
 *
 * @param ...
 */
void Fill_pulse_indices(gsl_vector *ind) {

	int i,k,k1,k2;

	/* Copy original vector */
	gsl_vector *ind_new = gsl_vector_alloc(ind->size);
	gsl_vector_memcpy(ind_new,ind);

	/* Fill */
	for(i=0;i<ind->size;i++) {
		if(gsl_vector_get(ind,i) == 0) {

			/* If gap is found, find the nearest non-zero value */
			if(gsl_vector_get(ind,i) == 0) {
				k1 = i+1;
				while(k1 < ind->size-1 && gsl_vector_get(ind,k1) == 0)
					k1++;
				k2 = i-1;
				while(k2 > 0 && gsl_vector_get(ind,k2) == 0)
					k2--;
				if(fabs(k1-i) < fabs(k2-i)) {
					if(k1 < ind->size-1)
						k = k1;
					else
						k = k2;
				} else {
					if(k2 > 0)
						k = k2;
					else
						k = k1;
				}

				/* Sanity check */
				if(k < 0 || k > ind->size-1)
					break;

				/* Fill the gap with the nearest found non-zero value */
				gsl_vector_set(ind_new,i,gsl_vector_get(ind,k));
			}
		}
	}

	/* Copy fixed vector to the original one and free memory */
	gsl_vector_memcpy(ind,ind_new);
	gsl_vector_free(ind_new);
}











/**
 * Function Evaluate_target_error
 *
 * Evaluate error between target parameters and pulse parameters.
 *
 * @param ...
 */
void Evaluate_target_error(gsl_matrix *lsf, gsl_matrix *glflowsp, gsl_matrix *harmonics, gsl_matrix *hnr_i, gsl_matrix *waveform, gsl_vector *h1h2, gsl_vector *naq, gsl_vector *gain, gsl_vector *fundf,
		gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *pharm, gsl_matrix *phnr, gsl_matrix *pwaveform, gsl_vector *ph1h2, gsl_vector *pnaq, gsl_vector *pgain, gsl_matrix *pca_w, gsl_matrix *pca_w_lib, gsl_vector *pulse_lengths,
		int frame_index, gsl_vector *inds, gsl_vector *eparam, gsl_vector *pulse_clus_id, gsl_matrix *pulse_clusters, PARAM *params) {

	int i,j,k,use_clusters = 0,npulses = params->number_of_pulses;
	int n_pulsecandidates = params->n_pulsecandidates;
	int nparams = params->paramweights->size;
	double m,std;
	gsl_vector_view cur_pulses;
	gsl_vector *lsfw = gsl_vector_calloc(lsf->size2);
	gsl_vector *lsfsourcew = gsl_vector_calloc(glflowsp->size2);
	gsl_vector_set_all(eparam, BIGGER_POS_NUMBER);

	/* Use pulse clusters */
	if(params->pulse_clustering == 1) {

		use_clusters = 1;
		cur_pulses = gsl_matrix_row(pulse_clusters, (int)(gsl_vector_get(pulse_clus_id, frame_index)));

		/* Get number of pulses in cluster (awkward) */
		npulses = 0;
		for(i=0;i<(&cur_pulses.vector)->size;i++) {
			if((int)gsl_vector_get(&cur_pulses.vector, i) == -1)
				break;
			npulses++;
		}
		if(npulses == 0) {
			use_clusters = 0;
			npulses = params->number_of_pulses;
		}
		if(params->n_pulsecandidates > npulses){
			n_pulsecandidates = npulses;
		}
	}

	/* Evaluate LSF weighting vector */
	for(i=1;i<lsf->size2-1;i++)
		gsl_vector_set(lsfw,i,1.0/(gsl_matrix_get(lsf,frame_index,i)-gsl_matrix_get(lsf,frame_index,i-1)) + 1.0/(gsl_matrix_get(lsf,frame_index,i+1)-gsl_matrix_get(lsf,frame_index,i)));
	double mean = Mean(lsfw);
	gsl_vector_set(lsfw,0,mean);
	gsl_vector_set(lsfw,lsfw->size-1,mean);

	/* Evaluate source LSF weighting vector */
	for(i=1;i<glflowsp->size2-1;i++)
		gsl_vector_set(lsfsourcew,i,1.0/(gsl_matrix_get(glflowsp,frame_index,i)-gsl_matrix_get(glflowsp,frame_index,i-1)) + 1.0/(gsl_matrix_get(glflowsp,frame_index,i+1)-gsl_matrix_get(glflowsp,frame_index,i)));
	mean = Mean(lsfsourcew);
	gsl_vector_set(lsfsourcew,0,gsl_vector_get(lsfsourcew,1));
	gsl_vector_set(lsfsourcew,lsfsourcew->size-1,mean);

	/* Initialize */
	gsl_matrix *error = gsl_matrix_calloc(npulses,nparams);
	gsl_vector *E = gsl_vector_calloc(npulses);
	gsl_permutation *p = gsl_permutation_alloc(E->size);

	/* FOR PULSE CLUSTERS! */
	if(use_clusters == 1) {

		/* Evaluate absolute error of every parameter: 0: LSF, 1:tilt, 2:harm, 3:hnr, 4:gain, 5:f0, 6:waveform, 7: h1h2, 8: naq*/
		for(k=0;k<npulses;k++){
			i = (int)gsl_vector_get((&cur_pulses.vector), k);
			if(i == -1)
				break;

			/* LSF */
			if(gsl_vector_get(params->paramweights,0) > 0) {
				for(j=0;j<plsf->size2;j++)
					gsl_matrix_set(error,k,0, gsl_matrix_get(error,k,0) + gsl_vector_get(lsfw,j)*powf(gsl_matrix_get(lsf,frame_index,j)-gsl_matrix_get(plsf,i,j),2));
				gsl_matrix_set(error,k,0,sqrt(gsl_matrix_get(error,k,0)));
			}

			/* Tilt */
			if(gsl_vector_get(params->paramweights,1) > 0 && params->use_tilt == 1) {
				for(j=0;j<ptilt->size2;j++)
				  gsl_matrix_set(error,k,1, gsl_matrix_get(error,k,1) + gsl_vector_get(lsfsourcew,j)*powf(gsl_matrix_get(glflowsp,frame_index,j)-gsl_matrix_get(ptilt,i,j),2));
				gsl_matrix_set(error,k,1,sqrt(gsl_matrix_get(error,k,1)));
			}

			/* Harmonics */
			if(gsl_vector_get(params->paramweights,2) > 0 && params->use_harmonics == 1) {
				for(j=0;j<pharm->size2;j++)
				  gsl_matrix_set(error,k,2, gsl_matrix_get(error,k,2) + powf(gsl_matrix_get(harmonics,frame_index,j)-gsl_matrix_get(pharm,i,j),2));
				gsl_matrix_set(error,k,2,sqrt(gsl_matrix_get(error,k,2)));
			}

			/* HNR */
			if(gsl_vector_get(params->paramweights,3) > 0 && params->use_hnr == 1) {
				for(j=0;j<phnr->size2;j++)
				  gsl_matrix_set(error,k,3, gsl_matrix_get(error,k,3) + powf(gsl_matrix_get(hnr_i,frame_index,j)-gsl_matrix_get(phnr,i,j),2));
				gsl_matrix_set(error,k,3,sqrt(gsl_matrix_get(error,k,3)));
			}

			/* Gain */
			if(gsl_vector_get(params->paramweights,4) > 0) {
				gsl_matrix_set(error,k,4,fabs(gsl_vector_get(gain,frame_index)-gsl_vector_get(pgain,i)));
			}

			/* F0 */
			if(gsl_vector_get(params->paramweights,5) > 0) {
				gsl_matrix_set(error,k,5,fabs(gsl_vector_get(fundf,frame_index)*params->pitch-2.0*params->FS/gsl_vector_get(pulse_lengths,i)));
			}

			/* Waveform */
			if(gsl_vector_get(params->paramweights,6) > 0 && params->use_waveform == 1) {
				for(j=0;j<pwaveform->size2;j++)
				  gsl_matrix_set(error,k,6, gsl_matrix_get(error,k,6) + powf(gsl_matrix_get(waveform,frame_index,j)-gsl_matrix_get(pwaveform,i,j),2));
				gsl_matrix_set(error,k,6,sqrt(gsl_matrix_get(error,k,6)));
			}

			/* H1H2 */
			if(gsl_vector_get(params->paramweights,7) > 0 && params->use_h1h2 == 1) {
				gsl_matrix_set(error,k,7,fabs(gsl_vector_get(h1h2,frame_index)-gsl_vector_get(ph1h2,i)));
			}

			/* NAQ */
			if(gsl_vector_get(params->paramweights,8) > 0 && params->use_naq == 1) {
				gsl_matrix_set(error,k,8,fabs(gsl_vector_get(naq,frame_index)-gsl_vector_get(pnaq,i)));
			}

			/* PCA */
			if(gsl_vector_get(params->paramweights,9) > 0 && params->use_pulselib_pca == 1) {
				for(j=0;j<pca_w->size2;j++)
				  gsl_matrix_set(error,k,9, gsl_matrix_get(error,k,9) + powf(gsl_matrix_get(pca_w,frame_index,j)-gsl_matrix_get(pca_w_lib,i,j),2));
				gsl_matrix_set(error,k,9,sqrt(gsl_matrix_get(error,k,9)));
			}
		}

		/* Normalize error of every parameter (subtract mean and divide by standard deviation) and weight */
		for(i=0;i<nparams;i++) {
			if(gsl_vector_get(params->paramweights,i) > 0) {

				/* Evaluate mean of each error */
				m = 0;
				for(j=0;j<npulses;j++)
					m += gsl_matrix_get(error,j,i);
				m = m/npulses;

				/* Evaluate standard deviation or each error */
				std = 0;
				for(j=0;j<npulses;j++)
					std += (gsl_matrix_get(error,j,i)-m)*(gsl_matrix_get(error,j,i)-m);
				std = sqrt(std/(npulses-1));
				if(std == 0)
					std = 1;

				/* Weight */
				for(j=0;j<npulses;j++)
					gsl_matrix_set(error,j,i,(gsl_matrix_get(error,j,i)-m)/std*gsl_vector_get(params->paramweights, i));
			}
		}

		/* Evaluate total error */
		double min_err = BIG_POS_NUMBER;
		int min_i = 0;
		gsl_vector *avg_err = gsl_vector_alloc(nparams);
		gsl_vector_set_zero(avg_err);
		for(i=0;i<npulses;i++){
			for(j=0;j<nparams;j++){
				gsl_vector_set(E,i,gsl_vector_get(E,i) + gsl_matrix_get(error,i,j));
				gsl_vector_set(avg_err, j, gsl_vector_get(avg_err, j) + gsl_matrix_get(error,i,j));
			}
			if (gsl_vector_get(E,i) < min_err){
				min_err = gsl_vector_get(E,i);
				min_i = i;
			}
		}
		gsl_vector_free(avg_err);
		gsl_sort_vector_index(p,E);
		for(i=0;i<n_pulsecandidates;i++) {
			gsl_vector_set(inds,i,gsl_permutation_get(p,i));
			gsl_vector_set(eparam,i,gsl_vector_get(E,gsl_vector_get(inds,i)));
			gsl_vector_set(inds,i,(int)gsl_vector_get(&cur_pulses.vector,gsl_permutation_get(p,i)));
		}

		/* Add best candidates to cluster of next frame */
		if(gsl_vector_get(pulse_clus_id, frame_index) != gsl_vector_get(pulse_clus_id, frame_index+1)) {
			if(npulses > 1000)
				npulses = 1000;
			for(i=0;i<n_pulsecandidates;i++) {
				gsl_matrix_set(pulse_clusters, (int)(gsl_vector_get(pulse_clus_id, frame_index+1)), npulses+i, gsl_vector_get(inds, i));
			}
		}

		/* Sanity check */
		int from_cluster = 0;
		for(i=0;i<(&cur_pulses.vector)->size;i++) {
			if(gsl_vector_get(&cur_pulses.vector, i) == gsl_vector_get(inds, 0))
				from_cluster = 1;
		}




	} else { /* FOR NORMAL ERROR CALCULATION */


		// TODO: SET LIMITS FOR EACH PARAMETER THAT PREVENT SELECTING A PULSE!
		// TODO: MAE OR RMSE

	    /* Evaluate RMSE of every parameter: 0: LSF, 1:tilt, 2:harm, 3:hnr, 4:gain, 5:f0, 6:waveform, 7:h1h2, 8: naq */
		for(i=0;i<npulses;i++) {

			/* LSF, with weighting */
			if(gsl_vector_get(params->paramweights,0) > 0) {
				for(j=0;j<plsf->size2;j++)
					gsl_matrix_set(error,i,0, gsl_matrix_get(error,i,0) + gsl_vector_get(lsfw,j)*powf(gsl_matrix_get(lsf,frame_index,j)-gsl_matrix_get(plsf,i,j),2));
				gsl_matrix_set(error,i,0,sqrt(gsl_matrix_get(error,i,0)));
			}

			/* Tilt, with weighting */
			if(gsl_vector_get(params->paramweights,1) > 0 && params->use_tilt == 1) {
				for(j=0;j<ptilt->size2;j++)
					gsl_matrix_set(error,i,1, gsl_matrix_get(error,i,1) + gsl_vector_get(lsfsourcew,j)*powf(gsl_matrix_get(glflowsp,frame_index,j)-gsl_matrix_get(ptilt,i,j),2));
				gsl_matrix_set(error,i,1,sqrt(gsl_matrix_get(error,i,1)));
			}

			/* Harmonics */
			if(gsl_vector_get(params->paramweights,2) > 0 && params->use_harmonics == 1) {
				for(j=0;j<pharm->size2;j++)
					gsl_matrix_set(error,i,2, gsl_matrix_get(error,i,2) + powf(gsl_matrix_get(harmonics,frame_index,j)-gsl_matrix_get(pharm,i,j),2));
				gsl_matrix_set(error,i,2,sqrt(gsl_matrix_get(error,i,2)));
			}

			/* HNR */
			if(gsl_vector_get(params->paramweights,3) > 0 && params->use_hnr == 1) {
				for(j=0;j<phnr->size2;j++)
					gsl_matrix_set(error,i,3, gsl_matrix_get(error,i,3) + powf(gsl_matrix_get(hnr_i,frame_index,j)-gsl_matrix_get(phnr,i,j),2));
				gsl_matrix_set(error,i,3,sqrt(gsl_matrix_get(error,i,3)));
			}

			/* Gain */
			if(gsl_vector_get(params->paramweights,4) > 0) {
				gsl_matrix_set(error,i,4,fabs(gsl_vector_get(gain,frame_index)-gsl_vector_get(pgain,i)));
			}

			/* F0 */
			if(gsl_vector_get(params->paramweights,5) > 0) {
				gsl_matrix_set(error,i,5,fabs(gsl_vector_get(fundf,frame_index)*params->pitch-2.0*params->FS/gsl_vector_get(pulse_lengths,i)));
			}

			/* Waveform */
			if(gsl_vector_get(params->paramweights,6) > 0 && params->use_waveform == 1) {
				for(j=0;j<pwaveform->size2;j++)
				  gsl_matrix_set(error,i,6, gsl_matrix_get(error,i,6) + powf(gsl_matrix_get(waveform,frame_index,j)-gsl_matrix_get(pwaveform,i,j),2));
				gsl_matrix_set(error,i,6,sqrt(gsl_matrix_get(error,i,6)));
			}

			/* H1H2 */
			if(gsl_vector_get(params->paramweights,7) > 0 && params->use_h1h2 == 1) {
				gsl_matrix_set(error,i,7,fabs(gsl_vector_get(h1h2,frame_index)-gsl_vector_get(ph1h2,i)));
			}

			/* Gain */
			if(gsl_vector_get(params->paramweights,8) > 0 && params->use_naq == 1) {
				gsl_matrix_set(error,i,8,fabs(gsl_vector_get(naq,frame_index)-gsl_vector_get(pnaq,i)));
			}

			/* PCA */
			if(gsl_vector_get(params->paramweights,9) > 0 && params->use_pulse_pca == 1) {
				for(j=0;j<pca_w->size2;j++)
				  gsl_matrix_set(error,i,9, gsl_matrix_get(error,i,9) + powf(gsl_matrix_get(pca_w,frame_index,j)-gsl_matrix_get(pca_w_lib,i,j),2));
				gsl_matrix_set(error,i,9,sqrt(gsl_matrix_get(error,i,9)));
			}
		}

	    /* Normalize error of every parameter (subtract mean and divide by standard deviation) and weight */
	    for(i=0;i<nparams;i++) {
	    	if(gsl_vector_get(params->paramweights,i) > 0) {

				/* Evaluate mean error of each parameter */
				m = 0;
				for(j=0;j<npulses;j++)
					m += gsl_matrix_get(error,j,i);
				m = m/npulses;

				/* Evaluate standard deviation of error of each parameter */
				std = 0;
				for(j=0;j<npulses;j++)
					std += (gsl_matrix_get(error,j,i)-m)*(gsl_matrix_get(error,j,i)-m);
				std = sqrt(std/(npulses-1));
				if(std == 0)
					std = 1;

				/* Normalize and weight */
				for(j=0;j<npulses;j++)
					gsl_matrix_set(error,j,i,(gsl_matrix_get(error,j,i)-m)/std*gsl_vector_get(params->paramweights,i));
	    	}
	    }

	    /* Evaluate total error, apply penalty for gross errors */
	    // TODO: DOES THIS YIELD BETTER QUALITY?
	    double gross_error_penalty_limit = 0; // For normally distributed data (mean = 0, sigma = 1)
	    double gross_error_penalty = 1.0;
		for(i=0;i<npulses;i++) {
			for(j=0;j<nparams;j++) {
				gsl_vector_set(E,i,gsl_vector_get(E,i) + gsl_matrix_get(error,i,j));
				if(gsl_matrix_get(error,i,j) > gross_error_penalty_limit)
					gsl_vector_set(E,i,gsl_vector_get(E,i) + gross_error_penalty);
			}
		}

		/* Select best candidates according to lowest target error and compare
		 * them to the previous pulse */
		gsl_sort_vector_index(p,E);
		for(i=0;i<n_pulsecandidates;i++) {
			gsl_vector_set(inds,i,gsl_permutation_get(p,i));
			gsl_vector_set(eparam,i,gsl_vector_get(E,gsl_vector_get(inds,i)));
		}
	}

	/* Free memory */
	gsl_matrix_free(error);
	gsl_vector_free(E);
	gsl_vector_free(lsfw);
	gsl_vector_free(lsfsourcew);
	gsl_permutation_free(p);
}




















/**
 * Function Evaluate_concatenation_error_viterbi
 *
 * Eval concatenation error for viterbi pulse search
 *
 * @param ...
 */
void Evaluate_concatenation_error_viterbi(gsl_vector *target, gsl_vector* target_next, gsl_vector *inds, gsl_vector *inds_next,gsl_matrix *pulses_rs,
		gsl_matrix *v_scores, gsl_matrix *v_best, int index, int dist_to_unvoiced, double additional_cost, PARAM *params) {

	int i, j, k = 0, num_skipped = 0;
	double e = 0,ccost;
	//double em = 0;

	/* Define concatenation cost according to distance to unvoiced and HNR */
	ccost = params->concatenation_cost*additional_cost*(1.0-(1.0/sqrt((double)(dist_to_unvoiced))));

	/* Current pulses */
	for(i=0;i<params->n_pulsecandidates;i++) {

		/* Next pulses */
		for(j=0;j<params->n_pulsecandidates;j++){

			/* Possible optimization: skip comparison if highly unlikely to beat the best score */
			if(index > 0) {
				double current_best = gsl_matrix_get(v_scores, index+1, j);
				if(ccost * 0.1 + params->target_cost * gsl_vector_get(target,i) + gsl_matrix_get(v_scores,index,i) >= current_best) {
					num_skipped++;
					continue;
				}
			}

		    	/* Compare downsampled waveforms */
			e = 0;
			for(k = 0;k<pulses_rs->size2;k++){
				e += powf(gsl_matrix_get(pulses_rs,gsl_vector_get(inds,i),k) - gsl_matrix_get(pulses_rs,gsl_vector_get(inds_next,j),k),2);
			}
			e = sqrt(e);

			/* Bias against same pulse */
			if(e == 0) e = params->pulse_error_bias;

			/* Combine the two error measures */
			//e = (e + em)/2.0;

			/* Evaluate combined score */
			double combined_score;
			if(index == 0) {
				combined_score  = params->target_cost * gsl_vector_get(target,i);
				gsl_matrix_set(v_scores, index,i, combined_score);
			} else
				combined_score = ccost * e + params->target_cost * gsl_vector_get(target,i) + gsl_matrix_get(v_scores, index, i);

			/* Set scores and indices */
			if(combined_score < gsl_matrix_get(v_scores, index+1, j)) {
				gsl_matrix_set(v_scores, index+1, j, combined_score);
				gsl_matrix_set(v_best, index+1, j, i);
			}
		}
	}
}















/**
 * Function Reestimate_hnr
 *
 * Re-estimate HNR for pulse library
 *
 * @param excitation_voiced
 * @param excitation_unvoiced
 * @param hnr HNR matrix
 * @param fundf F0 vector
 * @param params
 */
void Reestimate_hnr(gsl_vector *excitation_voiced, gsl_vector *excitation_unvoiced, gsl_matrix *hnr, gsl_vector *fundf, PARAM *params) {

	/* Only if pulse library and harmonics are used */
	if(params->use_pulselib == 0 && params->use_hnr == 1)
		return;

	/* Integrate voiced excitation signal */
	int i;
	gsl_vector *excitation_voiced_flow = gsl_vector_alloc(excitation_voiced->size);
	gsl_vector_memcpy(excitation_voiced_flow,excitation_voiced);
	for(i=1;i<excitation_voiced_flow->size;i++)
		gsl_vector_set(excitation_voiced_flow,i,gsl_vector_get(excitation_voiced_flow,i-1)*LEAK + gsl_vector_get(excitation_voiced_flow,i));

	/* Evaluate HNR */
	HNR_eval_vu(excitation_voiced_flow, excitation_unvoiced, hnr, fundf, params);
	Smooth_matrix(hnr,params->hnr_smooth_len);
	gsl_vector_free(excitation_voiced_flow);
}







/**
 * Function HNR_compensation
 *
 * Compensate HNR values based on the new values
 *
 * @param hnr old HNR values
 * @param hnr_new new HNR values
 */
void HNR_compensation(gsl_matrix *hnr, gsl_matrix *hnr_new, PARAM *params) {

	/* Set params */
	params->hnr_reestimated = 1;

	/* Return if HNR compensation is not used, or HNR is not used, or pulse library is used */
	if(params->hnr_compensation == 0 || hnr == NULL || params->use_pulselib == 1)
		return;

	/* Make compensation based on the estimated HNR of synthesized excitation */
	Smooth_matrix(hnr_new,params->hnr_smooth_len);
	int i,j;
	for(j=0;j<hnr->size1;j++) {
		for(i=0;i<hnr->size2;i++)
			gsl_matrix_set(hnr,j,i,gsl_matrix_get(hnr,j,i) + gsl_matrix_get(hnr,j,i) - gsl_matrix_get(hnr_new,j,i));
	}
}











/**
 * Function HNR_eval_vu
 *
 * Evaluate HNR of signal (voiced/unvoiced)
 *
 * @param signal_v input signal voiced
 * @param signal_uv input signal unvoiced
 * @param hnr HNR matrix
 * @param fundf fundamental frequency
 * @param frame_length frame length in samples
 * @paran shift shift length in samples
 * @param speed synthesis speed
 *
 */
void HNR_eval_vu(gsl_vector *signal_v,gsl_vector *signal_uv,gsl_matrix *hnr,gsl_vector *fundf, PARAM *params) {

	int i,j,add = rint((params->f0_frame_length/(double)params->shift/params->speed-1.0)*params->shift/params->speed/2.0);
	gsl_vector *frame = gsl_vector_alloc(rint(params->f0_frame_length/params->speed));
	gsl_vector *signal = gsl_vector_calloc(signal_v->size);
	gsl_vector_add(signal,signal_v);
	gsl_vector_add(signal,signal_uv);

	/* Zeropad signal */
	gsl_vector *signal_zp = gsl_vector_calloc(signal->size + 2*add);
	for(i=0;i<signal->size;i++)
		gsl_vector_set(signal_zp,i+add,gsl_vector_get(signal,i));

	/* Window and calculate glottal flow spectrum */
	for(i=0;i<hnr->size1;i++) {
		for(j=rint(i*params->shift/params->speed);j<GSL_MIN(rint(i*params->shift/params->speed+params->f0_frame_length/params->speed),signal_zp->size-1);j++) {
			gsl_vector_set(frame,GSL_MIN(j-rint(i*params->shift/params->speed),frame->size-1),gsl_vector_get(signal_zp,j));
		}
		Upper_lower_envelope(frame, hnr, gsl_vector_get(fundf,i), i, params);
	}
	gsl_vector_free(frame);
	gsl_vector_free(signal_zp);
	gsl_vector_free(signal);
}


































/**
 * Function Phase_manipulation
 *
 * Modify the phase and magnitude of the pulse to create noise.
 * The modified pulse is saved inplace to pulse.
 *
 * @param pulse original pulse
 * @param hnr HNR matrix
 * @param harmonics harmonic magnitudes
 * @param index current frame index
 * @param params parameter structure
 */
void Phase_manipulation(gsl_vector *pulse, gsl_matrix *hnr, gsl_matrix *harmonics, int index, PARAM *params) {

	/* Initialize */
	int i,n = pulse->size;
	double data[n];

	/* Set pulse data to array "data" */
	for(i=0;i<n;i++)
		data[i] = gsl_vector_get(pulse,i);

	/* FFT */
	gsl_fft_real_wavetable *wreal = gsl_fft_real_wavetable_alloc(n);
	gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(n);
	gsl_fft_real_transform(data,1,n,wreal,work);
	gsl_fft_real_wavetable_free(wreal);

	/* Extract real and imaginary parts to vectors */
	double x[2*n];
	gsl_complex_packed_array complex_coefficients = x;
	gsl_vector *real = gsl_vector_alloc(n);
	gsl_vector *imag = gsl_vector_alloc(n);
	gsl_fft_halfcomplex_unpack(data,complex_coefficients,1,n);
	for(i=0;i<n;i++) {
		gsl_vector_set(real,i,REAL(complex_coefficients,i));
		gsl_vector_set(imag,i,IMAG(complex_coefficients,i));
	}

	/* Get HNR values in ERB scale to vector, convert to Hz scale */
	gsl_vector *hnr_values_erb = gsl_vector_alloc(params->hnr_channels);
	gsl_vector *hnr_values = gsl_vector_alloc(ceil((imag->size-1)/2.0));
	for(i=0;i<params->hnr_channels;i++)
		gsl_vector_set(hnr_values_erb,i,gsl_matrix_get(hnr,index,i));
	Convert_ERB2Hz(hnr_values_erb,hnr_values,params);

	/* Noise power equals to pulse power minus the difference between harmonics and noise (indicated by HNR values) */
	gsl_vector *pulse_power = gsl_vector_alloc(hnr_values->size);
	for(i=0;i<hnr_values->size;i++)
		gsl_vector_set(pulse_power, i, 20*log10(sqrt(pow(gsl_vector_get(real,i+1), 2) + pow(gsl_vector_get(imag,i+1), 2))));
	MA(pulse_power,11);
	for(i=0;i<hnr_values->size;i++)
		gsl_vector_set(hnr_values, i, gsl_vector_get(hnr_values,i) + gsl_vector_get(pulse_power,i));
	gsl_vector_free(hnr_values_erb);
	gsl_vector_free(pulse_power);

	/* Convert HNR values from logarithmic scale to linear scale (actual noise amplitudes) */
	for(i=0;i<hnr_values->size;i++)
		gsl_vector_set(hnr_values,i,pow(10,(gsl_vector_get(hnr_values,i)/20.0)));

	/* Modification of the amplitudes of the harmonics */
	if(params->use_harmonic_modification != 0) {

		/* Extract radius and phase */
		gsl_vector *radius = gsl_vector_alloc(imag->size);
		gsl_vector *phase = gsl_vector_alloc(imag->size);
		for(i=0;i<radius->size;i++)
			gsl_vector_set(radius,i,sqrt(pow(gsl_vector_get(real,i),2) + pow(gsl_vector_get(imag,i),2)));
		for(i=0;i<phase->size;i++)
			gsl_vector_set(phase,i,atan2(gsl_vector_get(imag,i),gsl_vector_get(real,i)));

		/* Get pulse magnitude */
		gsl_vector *magnitude = gsl_vector_alloc(radius->size);
		for(i=0;i<magnitude->size;i++)
			gsl_vector_set(magnitude,i,20*log10(gsl_vector_get(radius,i)));

		/* Modify harmonics, use harmonics or modify by rule */
		if(params->use_harmonics == 1) { /* Use harmonics */

			/* Get harmonic log-magnitudes from matrix "harmonics" */
			gsl_vector *rnew = gsl_vector_alloc(params->number_of_harmonics+1);
			gsl_vector_set(rnew,0,20*log10(gsl_vector_get(radius,1)));
			for(i=0;i<params->number_of_harmonics;i++)
				gsl_vector_set(rnew,i+1,gsl_matrix_get(harmonics,index,i) + gsl_vector_get(rnew,0));

			/* Modify pulse magnitude */
			double diff_mag = gsl_vector_get(rnew,params->number_of_harmonics)-20.0*log10(gsl_vector_get(radius,params->number_of_harmonics+1));
			for(i=0;i<params->number_of_harmonics+1;i++)
				gsl_vector_set(magnitude,i+1,gsl_vector_get(rnew,i));
			for(i=params->number_of_harmonics+2;i<magnitude->size;i++)
				gsl_vector_set(magnitude,i,gsl_vector_get(magnitude,i)+diff_mag);

			/* Free memory */
			gsl_vector_free(rnew);

		} else if(params->pulse_tilt_decrease_coeff < 1) { /* Decrease spectral tilt by rule */

			/* Modify magnitude */
			double m0 = gsl_vector_get(magnitude,0);
			for(i=0;i<magnitude->size;i++)
				gsl_vector_set(magnitude,i,m0 + params->pulse_tilt_decrease_coeff*(gsl_vector_get(magnitude,i)-m0));
		}

		/* Set log-magnitude back to magnitude */
		for(i=0;i<radius->size;i++)
			gsl_vector_set(radius,i,pow(10,gsl_vector_get(magnitude,i)/20.0));

		/* Reconstruct imag and real */
		for(i=1;i<imag->size;i++) {
			gsl_vector_set(real,i,gsl_vector_get(radius,i)*cos(gsl_vector_get(phase,i)));
			gsl_vector_set(imag,i,gsl_vector_get(radius,i)*sin(gsl_vector_get(phase,i)));
		}

		/* Free memory */
		gsl_vector_free(radius);
		gsl_vector_free(phase);
		gsl_vector_free(magnitude);
	}

	/* Change noise low frequency limit from Hz to relative to FS */
	double noise_low_freq_limit_rel = params->noise_low_freq_limit/params->FS*2;

	/* Add noise by modifying both the magnitude and the phase of the spectrum of the pulse */
	for(i=rint(noise_low_freq_limit_rel*hnr_values->size);i<hnr_values->size;i++) {
		gsl_vector_set(imag,i+1,gsl_vector_get(imag,i+1) + params->noise_gain_voiced*RAND()*gsl_vector_get(hnr_values,i));
		gsl_vector_set(real,i+1,gsl_vector_get(real,i+1) + params->noise_gain_voiced*RAND()*gsl_vector_get(hnr_values,i));
	}

	/* Copy the noise to the folded spectrum as well */
	if(imag->size%2 == 0) {
		for(i=0;i<hnr_values->size;i++) {
			gsl_vector_set(imag,hnr_values->size+i,-gsl_vector_get(imag,hnr_values->size-i));
			gsl_vector_set(real,hnr_values->size+i,gsl_vector_get(real,hnr_values->size-i));
		}
	} else {
		for(i=0;i<hnr_values->size;i++) {
			gsl_vector_set(imag,hnr_values->size+1+i,-gsl_vector_get(imag,hnr_values->size-i));
			gsl_vector_set(real,hnr_values->size+1+i,gsl_vector_get(real,hnr_values->size-i));
		}
	}

	/* Create halfcomplex data */
	data[0] = gsl_vector_get(real,0);
	for(i=1;i<n;i++) {
		if(i%2 == 1)
			data[i] = gsl_vector_get(real,(i+1)/2);
		else
			data[i] = gsl_vector_get(imag,i/2);
	}

	/* Inverse FFT */
	gsl_fft_halfcomplex_wavetable *whc = gsl_fft_halfcomplex_wavetable_alloc(n);
	gsl_fft_halfcomplex_inverse(data,1,n,whc,work);
	gsl_fft_halfcomplex_wavetable_free(whc);
	gsl_fft_real_workspace_free(work);

	/* Set data from array "data" to vector "pulse" */
	for(i=0;i<n;i++)
		gsl_vector_set(pulse,i,data[i]);

	/* Free memory */
	gsl_vector_free(real);
	gsl_vector_free(imag);
	gsl_vector_free(hnr_values);
}









/**
 * Function Highpassfilter_fft
 *
 */
void Highpassfilter_fft(gsl_vector *signal) {

	/* Initialize */
	int i,n = signal->size;
	double data[n];

	/* Set signal to array "data" */
	for(i=0;i<n;i++)
		data[i] = gsl_vector_get(signal,i);

	/* FFT */
	gsl_fft_real_wavetable *wreal = gsl_fft_real_wavetable_alloc(n);
	gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(n);
	gsl_fft_real_transform(data,1,n,wreal,work);
	gsl_fft_real_wavetable_free(wreal);

	/* Extract real and imaginary parts to vectors */
	double x[2*n];
	gsl_complex_packed_array complex_coefficients = x;
	gsl_vector *real = gsl_vector_alloc(n);
	gsl_vector *imag = gsl_vector_alloc(n);
	gsl_fft_halfcomplex_unpack(data,complex_coefficients,1,n);
	for(i=0;i<n;i++) {
		gsl_vector_set(real,i,REAL(complex_coefficients,i));
		gsl_vector_set(imag,i,IMAG(complex_coefficients,i));
	}

	/* Extract radius and phase */
	gsl_vector *radius = gsl_vector_alloc(imag->size);
	gsl_vector *phase = gsl_vector_alloc(imag->size);
	for(i=0;i<radius->size;i++)
		gsl_vector_set(radius,i,sqrt(pow(gsl_vector_get(real,i),2) + pow(gsl_vector_get(imag,i),2)));
	for(i=0;i<phase->size;i++)
		gsl_vector_set(phase,i,atan2(gsl_vector_get(imag,i),gsl_vector_get(real,i)));

	/* Modify radius for high-pass filtering: Remove first half of the spectrum */
	int Nrem = rint(radius->size/4);
	for(i=0;i<Nrem;i++) {
		gsl_vector_set(radius,i+1,0);
		gsl_vector_set(radius,radius->size-1-i,0);
	}

	/* Reconstruct imag and real */
	for(i=1;i<imag->size;i++) {
		gsl_vector_set(real,i,gsl_vector_get(radius,i)*cos(gsl_vector_get(phase,i)));
		gsl_vector_set(imag,i,gsl_vector_get(radius,i)*sin(gsl_vector_get(phase,i)));
	}
	gsl_vector_set(imag,0,0);
	gsl_vector_set(real,0,0);

	/* Create halfcomplex data */
	data[0] = gsl_vector_get(real,0);
	for(i=1;i<n;i++) {
		if(i%2 == 1)
			data[i] = gsl_vector_get(real,(i+1)/2);
		else
			data[i] = gsl_vector_get(imag,i/2);
	}

	/* Inverse FFT */
	gsl_fft_halfcomplex_wavetable *whc = gsl_fft_halfcomplex_wavetable_alloc(n);
	gsl_fft_halfcomplex_inverse(data,1,n,whc,work);
	gsl_fft_halfcomplex_wavetable_free(whc);
	gsl_fft_real_workspace_free(work);

	/* Set data from array "data" to vector "signal" */
	for(i=0;i<n;i++)
		gsl_vector_set(signal,i,data[i]);

	/* Free memory */
	gsl_vector_free(radius);
	gsl_vector_free(phase);
	gsl_vector_free(real);
	gsl_vector_free(imag);
}





/**
 * Function Highpassfilter_fir
 *
 */
void Highpassfilter_fir(gsl_vector *signal, gsl_vector *coeffs) {

	/* Initialize */
	int i,j,n = coeffs->size;
	double sum;
	gsl_vector *temp = gsl_vector_calloc(signal->size+round(coeffs->size/2.0)-1);
	gsl_vector *signal_temp = gsl_vector_calloc(signal->size+round(coeffs->size/2.0)-1);
	for(i=0;i<signal->size;i++)
		gsl_vector_set(signal_temp,i,gsl_vector_get(signal,i));

	/* Filter signal */
	for(i=0; i<signal_temp->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, n-1); j++)
			sum += gsl_vector_get(signal_temp,i-j)*gsl_vector_get(coeffs,j);
		gsl_vector_set(temp, i, sum);
	}

	/* Copy "temp" samples to "signal" */
	for(i=0; i<signal->size; i++)
		gsl_vector_set(signal, i, gsl_vector_get(temp, i+coeffs->size/2));

	/* Free memory */
	gsl_vector_free(temp);
	gsl_vector_free(signal_temp);
}






/**
 * Function Convert_ERB2Hz
 *
 * Convert vector scale from ERB to Hz
 *
 * @param vector_erb pointer to vector of ERB-scale HNR values
 * @param vector pointer to reconstructed HNR vector
 *
 */
void Convert_ERB2Hz(gsl_vector *vector_erb, gsl_vector *vector, PARAM *params) {

	int i,j,hnr_channels = vector_erb->size;
	gsl_vector *erb = gsl_vector_alloc(vector->size);

	/* Evaluate ERB scale indices for vector */
	for(i=0;i<vector->size;i++)
		gsl_vector_set(erb,i,log10(0.00437*(i/(vector->size-1.0)*(params->FS/2.0))+1.0)/log10(0.00437*(params->FS/2.0)+1.0)*(hnr_channels-SMALL_NUMBER));

	/* Evaluate values according to ERB rate, smooth */
	for(i=0;i<vector->size;i++) {
		j = floor(gsl_vector_get(erb,i));
		gsl_vector_set(vector,i,gsl_vector_get(vector_erb,j));
	}
	MA(vector,3);

	/* Free memory */
	gsl_vector_free(erb);
}





















/**
 * Function Upper_lower_envelope
 *
 * Extract the amount of aperiodicity according to smoothed upper and lower spectral envelopes.
 *
 * @param frame pointer to frame
 * @param hnr pointer to matrix where the results are saved
 * @param f0 fundamental frequency
 * @param index time index
 *
 */
void Upper_lower_envelope(gsl_vector *frame, gsl_matrix *hnr, double f0, int index, PARAM *params) {

	/* Define FFT length */
	int FFT_LENGTH = MIN_FFT_LENGTH;
	while(FFT_LENGTH < frame->size)
		FFT_LENGTH = FFT_LENGTH*2;

	/* Variables */
	int i,j,guess_index = 0,h_values_size,harmonic_search_range;
	int hnr_channels = hnr->size2;
	double ind[MAX_HARMONICS] = {0};
	double guess_index_double = 0;
	double data[FFT_LENGTH];
	double h_values[MAX_HARMONICS] = {0};
	double n_values[MAX_HARMONICS] = {0};
	gsl_vector *fft = gsl_vector_alloc(FFT_LENGTH/2);
	gsl_vector *find;

	/* Initialize data */
	for(i=0;i<FFT_LENGTH;i++)
		data[i] = 0;

	/* FFT (with windowing) */
	for (i=0; i<frame->size; i++)
		data[i] = gsl_vector_get(frame, i)*HANN(i,frame->size);
	gsl_fft_real_radix2_transform(data, 1, FFT_LENGTH);
	for(i=1; i<FFT_LENGTH/2; i++)
		gsl_vector_set(fft, i, 20*log10(sqrt(pow(data[i], 2) + pow(data[FFT_LENGTH-i], 2))));
	gsl_vector_set(fft, 0, 20*log10(abs(data[0])));

	/* Set (possible) infinity values to zero */
	for(i=0;i<fft->size;i++) {
		if(!gsl_finite(gsl_vector_get(fft,i)))
			gsl_vector_set(fft,i,MIN_LOG_POWER);
	}

	/* Find the indices and magnitudes of the harmonic peaks */
	i = 0;
	while(1) {

		/* Define harmonics search range, decreasing to the higher frequencies */
		harmonic_search_range = GSL_MAX(HARMONIC_SEARCH_COEFF*f0/(params->FS/(double)FFT_LENGTH)*((fft->size-1-guess_index)/(double)(fft->size-1)),1.0);
		find = gsl_vector_alloc(harmonic_search_range);

		/* Estimate the index of the i_th harmonic
		 * Use an iterative estimation based on earlier values */
		if(i > 0) {
			guess_index_double = 0;
			for(j=0;j<i;j++)
				guess_index_double += ind[j]/(j+1.0)*(i+1.0);
			guess_index = (int)GSL_MAX(guess_index_double/j - (harmonic_search_range-1)/2.0,0);
		} else
			guess_index = (int)GSL_MAX(f0/(params->FS/(double)FFT_LENGTH) - (harmonic_search_range-1)/2.0,0);

		/* Stop search if the end of the fft vector or the maximum number of harmonics is reached */
		if(guess_index + rint(HNR_UNCERTAINTY_COEFF*FFT_LENGTH) > fft->size-1 || i > MAX_HARMONICS-1) {
			gsl_vector_free(find);
			break;
		}

		/* Find the maximum of the i_th harmonic */
		for(j=0; j<harmonic_search_range; j++) {
			if(guess_index+j < fft->size)
				gsl_vector_set(find, j, gsl_vector_get(fft, guess_index+j));
			else
				gsl_vector_set(find, j, BIG_NEG_NUMBER);
		}
		ind[i] = guess_index + gsl_vector_max_index(find);
		h_values[i] = gsl_vector_get(fft, ind[i]);
		gsl_vector_free(find);
		i++;
	}
	h_values_size = i;

	/* Estimate the level of interharmonic noise */
	double ind_n[MAX_HARMONICS] = {0};
	for(i=0;i<h_values_size-1;i++) {

		/* Evaluate value exactly between the harmonics */
		ind_n[i] = rint((ind[i]+ind[i+1])/2.0);
		n_values[i] = gsl_vector_get(fft,ind_n[i]);
	}

	/* Postfilter and iterpolate vectors */
	gsl_vector *hnr_est = gsl_vector_alloc(h_values_size-1);
	gsl_vector *hnr_est_erb = gsl_vector_calloc(hnr_channels);
	for(i=0;i<hnr_est->size;i++)
		gsl_vector_set(hnr_est,i,n_values[i]-h_values[i]);
	MedFilt3(hnr_est);
	Convert_Hz2ERB(hnr_est, hnr_est_erb, params);

	/* Set values to hnr matrix */
	for(i=0;i<hnr_channels;i++)
		gsl_matrix_set(hnr,index,i,gsl_vector_get(hnr_est_erb,i));

	/* Free memory */
	gsl_vector_free(fft);
	gsl_vector_free(hnr_est);
	gsl_vector_free(hnr_est_erb);
}










/**
 * Function Convert_Hz2ERB
 *
 * Convert vector scale from Hz to ERB
 *
 * @param vector pointer to original HNR vector
 * @param vector_erb pointer to vector of ERB-scale HNR values
 *
 */
void Convert_Hz2ERB(gsl_vector *vector, gsl_vector *vector_erb, PARAM *params) {

	int i,j,hnr_channels = vector_erb->size;
	gsl_vector *erb = gsl_vector_alloc(vector->size);
	gsl_vector *erb_sum = gsl_vector_calloc(hnr_channels);

	/* Evaluate ERB scale indices for vector */
	for(i=0;i<vector->size;i++)
		gsl_vector_set(erb,i,log10(0.00437*(i/(vector->size-1.0)*(params->FS/2.0))+1.0)/log10(0.00437*(params->FS/2.0)+1.0)*(hnr_channels-SMALL_NUMBER));

	/* Evaluate values according to ERB rate */
	for(i=0;i<vector->size;i++) {
		j = floor(gsl_vector_get(erb,i));
		gsl_vector_set(vector_erb,j,gsl_vector_get(vector_erb,j)+gsl_vector_get(vector,i));
		gsl_vector_set(erb_sum,j,gsl_vector_get(erb_sum,j)+1);
	}

	/* Average values */
	for(i=0;i<hnr_channels;i++)
			gsl_vector_set(vector_erb,i,gsl_vector_get(vector_erb,i)/gsl_vector_get(erb_sum,i));

	/* Prevent NaN-values (due to division by zero) */
	for(i=0;i<hnr_channels;i++) {
		if(gsl_vector_get(erb_sum,i) == 0) {
			j = 1;
			while(gsl_vector_get(erb_sum,i+j) == 0)
				j++;
			gsl_vector_set(vector_erb,i,0.5*gsl_vector_get(vector_erb,i-1)+0.5*gsl_vector_get(vector_erb,i+j));
		}
	}

	/* Free memory */
	gsl_vector_free(erb);
	gsl_vector_free(erb_sum);
}









/**
 * Function Create_pulse_train
 *
 * Reconstruct voice source in a frame
 *
 * @param ...
 */
gsl_vector *Create_pulse_train(gsl_vector *pulse, gsl_vector *original_pulse, gsl_matrix *hnr, gsl_vector *fundf, gsl_matrix *harmonics, int frame_index, int sample_index, PARAM *params) {

	/* Initialize */
	int i,frame_index_new,sample_index_new,N,pulse_start,pulse_end;
	gsl_vector *pulsetrain = gsl_vector_calloc(params->f0_frame_length+pulse->size);
	gsl_vector *pulse_copy = gsl_vector_alloc(pulse->size);
	int ptlen = pulsetrain->size;

	/* Set first pulse to pulsetrain in the middle of the frame */
	gsl_vector_memcpy(pulse_copy,pulse);
	Phase_manipulation(pulse_copy,hnr,harmonics,frame_index,params);
	for(i=0;i<pulse->size;i++)
		gsl_vector_set(pulsetrain,i+ptlen/2.0-pulse->size/2.0,gsl_vector_get(pulse_copy,i));
	gsl_vector_free(pulse_copy);
	pulse_start = ptlen/2.0-pulse->size/2.0+1;
	pulse_end = ptlen/2.0-pulse->size/2.0+i-1;

	/* Set earlier pulses to pulsetrain */
	frame_index_new = frame_index;
	sample_index_new = sample_index;
	N = pulse->size;
	while(1) {
		sample_index_new = GSL_MAX(sample_index_new - N,0);
		frame_index_new = GSL_MIN(rint(params->n_frames*((sample_index_new+0.5*N)/(double)params->signal_length)),params->n_frames-1);
		N = rint(params->FS/gsl_vector_get(fundf,frame_index_new)/params->pitch);
		if(gsl_vector_get(fundf,frame_index_new)/params->pitch == 0) {
			N = params->shift/params->speed;
			pulse_start = pulse_start-N;
			if(pulse_start-N >= 0)
				continue;
			else
				break;
		}
		gsl_vector *pulse_new = gsl_vector_alloc(N);
		Interpolate(original_pulse,pulse_new);
		Phase_manipulation(pulse_new,hnr,harmonics,frame_index_new,params);
		for(i=0;i<N;i++)
			gsl_vector_set(pulsetrain,GSL_MAX(pulse_start-N+i,0),gsl_vector_get(pulse_new,i));
		pulse_start = pulse_start-N+1;
		gsl_vector_free(pulse_new);
		if(pulse_start < 0)
			break;
	}

	/* Set later pulses to pulsetrain */
	frame_index_new = frame_index;
	sample_index_new = sample_index;
	N = pulse->size;
	while(1) {
		sample_index_new = GSL_MIN(sample_index_new + N,params->signal_length-1);
		frame_index_new = GSL_MIN(rint(params->n_frames*((sample_index_new-0.5*N)/(double)params->signal_length)),params->n_frames-1);
		N = rint(params->FS/gsl_vector_get(fundf,frame_index_new)/params->pitch);
		if(gsl_vector_get(fundf,frame_index_new)/params->pitch == 0) {
			N = params->shift/params->speed;
			pulse_end = pulse_end+N;
			if(pulse_end+N < pulsetrain->size)
				continue;
			else
				break;
		}
		gsl_vector *pulse_new = gsl_vector_alloc(N);
		Interpolate(original_pulse,pulse_new);
		Phase_manipulation(pulse_new,hnr,harmonics,frame_index_new,params);
		for(i=0;i<N;i++) {
			if(pulse_end+i > pulsetrain->size-1)
				break;
			gsl_vector_set(pulsetrain,pulse_end+i,gsl_vector_get(pulse_new,i));
		}
		pulse_end = pulse_end+N-1;
		gsl_vector_free(pulse_new);
		if(pulse_end > pulsetrain->size-1)
			break;
	}
	/* Return result */
	return pulsetrain;
}



/**
 * Function Create_pulse_train_diff
 *
 * Reconstruct voice source in a frame (for differentiated two-pitch-period pulses)
 *
 * @param ...
 */
gsl_vector *Create_pulse_train_diff(gsl_vector *pulse, gsl_vector *original_pulse, gsl_matrix *hnr, gsl_vector *fundf, gsl_matrix *harmonics, int frame_index, int sample_index, PARAM *params) {

	/* Initialize */
	int i,frame_index_new,sample_index_new,N,pulse_start,pulse_end;
	gsl_vector *pulsetrain = gsl_vector_calloc(params->f0_frame_length+round(pulse->size/2));
	gsl_vector *pulse_copy = gsl_vector_alloc(pulse->size);
	int ptlen = pulsetrain->size;

	/* Add noise to pulse copy */
	gsl_vector_memcpy(pulse_copy,pulse);
	Integrate(pulse_copy,LEAK);
	Phase_manipulation(pulse_copy,hnr,harmonics,frame_index,params);
	Differentiate(pulse_copy,LEAK);
	gsl_vector_set(pulse_copy,0,0);

	/* Set first pulse to pulsetrain in the middle of the frame */
	for(i=0;i<GSL_MIN(pulse->size,ptlen);i++)
		gsl_vector_set(pulsetrain,GSL_MAX(i+ptlen/2.0-pulse->size/2.0,0),gsl_vector_get(pulse_copy,i));
	gsl_vector_free(pulse_copy);
	pulse_start = ptlen/2.0-pulse->size/2.0+1;
	pulse_end = ptlen/2.0-pulse->size/2.0+i-1;

	/* Set earlier pulses to pulsetrain */
	frame_index_new = frame_index;
	sample_index_new = sample_index;
	N = pulse->size;
	while(1) {
		sample_index_new = GSL_MAX(sample_index_new - round(N/2),0);
		frame_index_new = GSL_MIN(rint(params->n_frames*((sample_index_new+0.25*N)/(double)params->signal_length)),params->n_frames-1);
		N = rint(2.0*params->FS/gsl_vector_get(fundf,frame_index_new)/params->pitch);
		if(gsl_vector_get(fundf,frame_index_new)/params->pitch == 0) {
			N = params->shift/params->speed;
			pulse_start = pulse_start-N;
			if(pulse_start-N >= 0)
				continue;
			else
				break;
		}
		gsl_vector *pulse_new = gsl_vector_alloc(N);
		Interpolate(original_pulse,pulse_new);

		/* Add noise */
		Integrate(pulse_new, LEAK);
		Phase_manipulation(pulse_new,hnr,harmonics,frame_index_new,params);
		Differentiate(pulse_new,LEAK);
		gsl_vector_set(pulse_new,0,0); // Fix DC error

		/* Set to excitation */
		for(i=0;i<N;i++)
			gsl_vector_set(pulsetrain,GSL_MAX(pulse_start-round(N/2)+i,0),gsl_vector_get(pulsetrain,GSL_MAX(pulse_start-round(N/2)+i,0)) + gsl_vector_get(pulse_new,i));
		pulse_start = pulse_start-round(N/2+1);
		gsl_vector_free(pulse_new);
		if(pulse_start < 0)
			break;
	}

	/* Set later pulses to pulsetrain */
	frame_index_new = frame_index;
	sample_index_new = sample_index;
	N = pulse->size;
	while(1) {
		sample_index_new = GSL_MIN(sample_index_new + round(N/2),params->signal_length-1);
		frame_index_new = GSL_MIN(rint(params->n_frames*((sample_index_new-0.25*N)/(double)params->signal_length)),params->n_frames-1);
		N = rint(2*params->FS/gsl_vector_get(fundf,frame_index_new)/params->pitch);
		if(gsl_vector_get(fundf,frame_index_new)/params->pitch == 0) {
			N = params->shift/params->speed;
			pulse_end = pulse_end+N;
			if(pulse_end+N < pulsetrain->size)
				continue;
			else
				break;
		}
		gsl_vector *pulse_new = gsl_vector_alloc(N);
		Interpolate(original_pulse,pulse_new);

		/* Add noise */
		Integrate(pulse_new, LEAK);
		Phase_manipulation(pulse_new,hnr,harmonics,frame_index_new,params);
		Differentiate(pulse_new,LEAK);
		gsl_vector_set(pulse_new,0,0); // Fix DC error

		/* Set to excitation */
		for(i=0;i<N;i++) {
			if(pulse_end-round(N/2)+i > pulsetrain->size-1 || pulse_end-round(N/2)+i < 0)
				break;
			gsl_vector_set(pulsetrain,pulse_end-round(N/2)+i,gsl_vector_get(pulsetrain,pulse_end-round(N/2)+i) + gsl_vector_get(pulse_new,i));
		}
		pulse_end = pulse_end+round(N/2)-1;
		gsl_vector_free(pulse_new);
		if(pulse_end > pulsetrain->size-1)
			break;
	}

	/* Return result */
	return pulsetrain;
}








/**
 * Function Analyse_pulse_train_spectrum
 *
 * Analyse LP spectrum of the pulse train
 *
 * @param ...
 */
void Analyse_pulse_train_spectrum(gsl_vector *pulse_train, gsl_matrix *glflowsp_new, int N, int sample_index, int frame_index, PARAM *params) {

	int i;
	gsl_vector *a;

	/* Apply pre-emphasis */
	if(USE_PRE_EMPH == 1) {
		gsl_vector *e_pulse_train = gsl_vector_alloc(pulse_train->size);
		gsl_vector_set(e_pulse_train,0,gsl_vector_get(pulse_train,0));
		for(i=1;i<pulse_train->size;i++)
			gsl_vector_set(e_pulse_train,i,gsl_vector_get(pulse_train,i) - 0.99*gsl_vector_get(pulse_train,i-1));
		a = WLPC(e_pulse_train, params->lpc_order_gl, params->lambda_gl);
		gsl_vector_free(e_pulse_train);
	} else {
		a = WLPC(pulse_train, params->lpc_order_gl, params->lambda_gl);
	}

	/* Set spectrum to matrix, fill in empty values */
	for(i=frame_index;i<GSL_MIN(floor(params->n_frames*((sample_index+N)/(double)params->signal_length)),params->n_frames);i++)
		Convert_to_LSF(glflowsp_new,a,i);

	/* Free memory */
	gsl_vector_free(a);
}
















/**
 * Function Mean
 *
 * Evaluate mean of a vector
 *
 * @param vector pointer to vector
 * @return mean
 *
 */
double Mean(gsl_vector *vector) {

	int i;
	double mean = 0;
	for(i=0; i<vector->size; i++) {
		mean += gsl_vector_get(vector, i);
	}
	return (mean/(double)vector->size);
}

/**
 * Function NonZeroMean
 *
 * Evaluate mean of a vector from elements that are not zeros
 *
 * @param vector pointer to vector
 * @return mean
 *
 */
double NonZeroMean(gsl_vector *vector) {

	int i,ind = 0;
	double mean = 0;
	for(i=0; i<vector->size; i++) {
		if(gsl_vector_get(vector,i) != 0) {
			mean += gsl_vector_get(vector, i);
			ind++;
		}
	}
	return (mean/(double)ind);
}








/**
 * Function FillHNRValues
 *
 * Fill missing HNR values
 *
 * @param hnr HNR values
 * @param frame_index_old old index
 * @param frame_index current index
 */
void FillHNRValues(gsl_matrix *hnr, int frame_index_old, int frame_index) {

	int i,j;
	for(i=frame_index_old+1;i<frame_index;i++) {
		for(j=0;j<hnr->size2;j++)
			gsl_matrix_set(hnr,i,j,gsl_matrix_get(hnr,frame_index,j));
	}
}






/**
 * Function FillUnvoicedSyntheticSourceSpectrum
 *
 * Fill current and previous source spectrum values
 *
 * @param glflowsp original spectrum
 * @param glflowsp_new new spectrum
 * @param frame_index_old old index
 * @param frame_index current index
 */
void FillUnvoicedSyntheticSourceSpectrum(gsl_matrix *glflowsp, gsl_matrix *glflowsp_new, int frame_index, int frame_index_old) {

	if(glflowsp == NULL)
		return;

	/* Set */
	int i,j;
	for(i=0;i<glflowsp->size2;i++)
		gsl_matrix_set(glflowsp_new, frame_index, i, gsl_matrix_get(glflowsp, frame_index, i));
	for(j=frame_index_old+1;j<frame_index;j++) {
		if(gsl_matrix_get(glflowsp_new,j,0) == 0) {
			for(i=0;i<glflowsp->size2;i++)
				gsl_matrix_set(glflowsp_new, j, i, gsl_matrix_get(glflowsp, j, i));
		}
	}

	/* Double check */
	if(frame_index+1 < glflowsp_new->size1)
		for(i=0;i<glflowsp_new->size2;i++)
			gsl_matrix_set(glflowsp_new, frame_index+1, i, gsl_matrix_get(glflowsp, frame_index+1, i));
}



/**
 * Function FillUnvoicedSyntheticSourceSpectrum_FLAT
 *
 * Fill current and previous source spectrum values
 *
 * @param glflowsp original spectrum
 * @param glflowsp_new new spectrum
 * @param frame_index_old old index
 * @param frame_index current index
 */
void FillUnvoicedSyntheticSourceSpectrum_FLAT(gsl_matrix *glflowsp, gsl_matrix *glflowsp_new, int frame_index, int frame_index_old) {

	if(glflowsp == NULL)
		return;

	/* Set */
	int i,j;
	for(i=0;i<glflowsp->size2;i++)
		gsl_matrix_set(glflowsp_new, frame_index, i, (i+1.0)*(M_PI/glflowsp->size2));
	for(j=frame_index_old+1;j<frame_index;j++) {
		if(gsl_matrix_get(glflowsp_new,j,0) == 0) {
			for(i=0;i<glflowsp->size2;i++)
				gsl_matrix_set(glflowsp_new, j, i, (i+1.0)*(M_PI/glflowsp->size2));
		}
	}

	/* Double check */
	if(frame_index+1 < glflowsp_new->size1)
		for(i=0;i<glflowsp_new->size2;i++)
			gsl_matrix_set(glflowsp_new, frame_index+1, i, gsl_matrix_get(glflowsp, frame_index+1, i));
}













/**
 * Function MedFilt5_matrix
 *
 * 5-point median filtering for matrices.
 * Filter along the first dimension.
 *
 * @param frame pointer to matrix to be filtered
 *
 */
void MedFilt5_matrix(gsl_matrix *matrix) {

	int i,j;
	gsl_vector *temp = gsl_vector_alloc(matrix->size1);
	for(i=0;i<matrix->size2;i++) {
		for(j=0;j<matrix->size1;j++) {
			gsl_vector_set(temp,j,gsl_matrix_get(matrix,j,i));
		}
		MedFilt5(temp);
		for(j=0;j<matrix->size1;j++) {
			gsl_matrix_set(matrix,j,i,gsl_vector_get(temp,j));
		}
	}
	gsl_vector_free(temp);
}

















/**
 * Function Modify_harmonics
 *
 * Modification of the harmonic information for speech in the presence of noise
 *
 * @param harmonics pointer to harmonics matrix
 * @param fundf pointer to F0 vector
 *
 */
void Modify_harmonics(gsl_matrix *harmonics, double scale) {

	/* Decrease the tilt of the spectrum */
	gsl_matrix_scale(harmonics, scale);

}





/**
 * Function Compression
 *
 * Dynamic range compression
 * Nonlinear power of k, 0<k<1, and linear response below threshold
 *
 * @param signal pointer to signal
 *
 */
void Compression(gsl_vector *signal, double k) {

	/* Check validity */
	if(k >= 1 || k <= 0)
		return;

	/* Dynamic range compression */
	double th = powf((1.0/k),(1.0/(k-1.0))); // Point where the derivative of the curve is one -> start nonlinearity after this threshold
	double a = powf(th,k) - th;              // The difference at the the turning point
	int i;
	for(i=0;i<signal->size;i++) {
		if(gsl_vector_get(signal,i) < -th)
			gsl_vector_set(signal,i,-powf(-gsl_vector_get(signal,i),k) + a);
		else if(gsl_vector_get(signal,i) > th)
			gsl_vector_set(signal,i,powf(gsl_vector_get(signal,i),k) - a);
	}

	/* Scale maximum to one */
	double absmax = GSL_MAX(gsl_vector_max(signal),-gsl_vector_min(signal));
	for(i=0;i<signal->size;i++)
		gsl_vector_set(signal,i,WAV_SCALE*gsl_vector_get(signal,i)/absmax);
}
















/**
 * Function MedFilt5
 *
 * 5-point median filtering
 *
 * @param frame pointer to vector to be filtered
 *
 */
void MedFilt5(gsl_vector *frame) {

	int i,j,n = 5;
	double temp1 = gsl_vector_get(frame, 0);
	double temp2 = gsl_vector_get(frame,1);
	double temp3;
	gsl_vector *samples = gsl_vector_alloc(n);

	for(i=0; i<frame->size-(n-1); i++) {
		for(j=0; j<n; j++) {
			gsl_vector_set(samples, j, gsl_vector_get(frame, i+j));
		}
		gsl_sort_vector(samples);
		temp3 = gsl_vector_get(samples, 2);
		gsl_vector_set(frame, i, temp1);
		temp1 = temp2;
		temp2 = temp3;
	}
	gsl_vector_free(samples);
}




/**
 * Function MedFilt3
 *
 * 3-point median filtering
 *
 * @param frame pointer to vector to be filtered
 *
 */
void MedFilt3(gsl_vector *frame) {

	int i,j,n;
	double temp1 = 0;
	double temp2 = gsl_vector_get(frame, 0);
	n = 3;
	gsl_vector *samples = gsl_vector_alloc(n);

	for(i=0; i<frame->size-2; i++) {
		for(j=0; j<n; j++) {
			gsl_vector_set(samples, j, gsl_vector_get(frame, i+j));
		}
		gsl_sort_vector(samples);
		temp1 = gsl_vector_get(samples, 1);
		gsl_vector_set(frame, i, temp2);
		temp2 = temp1;
	}
	gsl_vector_free(samples);
}

























/**
 * Function Smooth_interp_lsf
 *
 * Smooth and interpolate LSF-matrix
 *

 * @param LSF_interp pointer to interpolated LSF-matrix
 * @param LSF original LSF-matrix
 * @param signal_len signal length
 * @param use_hmm HMM-switch
 * @param lsf_smooth_len smoothing length in samples
 *
 */
void Smooth_interp_lsf(gsl_matrix *LSF_i, gsl_matrix *LSF, int signal_len, int use_hmm, int lsf_smooth_len) {

	int i,j;
	gsl_vector *temp = gsl_vector_alloc(LSF->size1);
	gsl_vector *temp_i = gsl_vector_alloc(signal_len);

	for(i=0;i<LSF->size2;i++) {
		for(j=0;j<LSF->size1;j++) {
			gsl_vector_set(temp,j,gsl_matrix_get(LSF,j,i));
		}

		/* Smooth vectors in time */
		if(use_hmm == 0)
			MA(temp,lsf_smooth_len);

		/* Interpolate vector */
		Interpolate(temp,temp_i);

		/* Set smoothed and interpolated LSFs to matrices */
		for(j=0;j<signal_len;j++) {
			gsl_matrix_set(LSF_i,j,i,gsl_vector_get(temp_i,j));
		}
	}

	/* Free memory */
	gsl_vector_free(temp);
	gsl_vector_free(temp_i);
}





/**
 * Function Filter_excitation
 *
 * Filter excitation (normal/warped)
 *
 * @param excitation excitation
 * @param LSF lsf matrix
 * @param params
 *
 */
void Filter_excitation(gsl_vector *excitation, gsl_matrix *LSF, PARAM *params) {

	/* Normal filtering */
	if(params->lambda_vt == 0) {

		int i,j;
		double sum;

		gsl_vector *A = gsl_vector_alloc(params->lpc_order_vt+1);
		for(i=0;i<params->signal_length;i++) {

			/* Update filter coeffs: convert LSF to poly */
			if(i%params->filter_update_interval_vt == 0) {
		        lsf2poly(LSF,A,i,params->use_hmm);
				for(j=0;j<params->lpc_order_vt+1;j++) {
					gsl_vector_set(A,j,gsl_vector_get(A,j)*(-1));
				}
				gsl_vector_set(A,0,1);
			}

			/* Filter */
	        sum = 0;
	        for(j=0;j<GSL_MIN(params->lpc_order_vt+1,i);j++) {
	        	sum += gsl_vector_get(excitation,i-j)*gsl_vector_get(A,j);
	        }
	        gsl_vector_set(excitation,i,sum);
		}
		gsl_vector_free(A);

	/* Warped filtering */
	} else {

		int i,q,mlen;
	    long int o;
	    double xr,x,ffr,tmpr,Bb;
	    double *sigma;
	    long int len = params->signal_length;
	    int adim = 1;
	    int bdim = LSF->size2 + 1;
	    double *Ar = (double *)calloc(adim,sizeof(double));
	    double *Br = (double *)calloc(bdim,sizeof(double));
	    double *ynr = (double *)calloc(excitation->size,sizeof(double));
	    double *rsignal = (double *)calloc(excitation->size,sizeof(double));
	    double *rmem = (double *)calloc(GSL_MAX(adim,bdim)+2,sizeof(double));
	    gsl_vector *B = gsl_vector_alloc(bdim);

	    /* Set excitation to array */
	    for(i=0;i<len;i++) {
			rsignal[i] = gsl_vector_get(excitation,i);
		}

	    /* Initialize */
	    Ar[0] = 1;
	    Bb = 0;
	    sigma = NFArray(bdim+2);
	    if(adim >= bdim)
	    	mlen = adim;
	    else
	    	mlen = bdim + 1;

	    /* Warped filtering */
	    for(o=0;o<len;o++) {

	    	/* Update filter coefficients */
	    	if(o%params->filter_update_interval_vt == 0) {
	    		lsf2poly(LSF,B,o,params->use_hmm);
	    		for(i=0;i<bdim;i++) {
	    			Br[i] = gsl_vector_get(B,i);
	    		}
	    		alphas2sigmas(Br,sigma,params->lambda_vt,bdim-1);
	    		Bb = 1/Br[0];
	    	}

	    	xr = rsignal[o]*Bb;

	    	/* Update feedbackward sum */
	    	for(q=0;q<bdim;q++) {
	    		xr -= sigma[q]*rmem[q];
	    	}
	    	xr = xr/sigma[bdim];
	    	x = xr*Ar[0];

	    	/* Update inner states */
	    	for(q=0;q<mlen;q++) {
	    		tmpr = rmem[q] + params->lambda_vt*(rmem[q+1] - xr);
	    		rmem[q] = xr;
	    		xr = tmpr;
	    	}

	    	/* Update feedforward sum */
	    	for(q=0,ffr=0.0;q<adim-1;q++) {
	    		ffr += Ar[q+1]*rmem[q+1];
	    	}

	       /* Update output */
	       ynr[o] = x + ffr;
		}

	    /* Set output to vector */
	    for(i=0;i<len;i++) {
	    	gsl_vector_set(excitation,i,ynr[i]);
	    }

	    /* Free memory */
		free(ynr);
		free(rsignal);
		free(Ar);
		free(Br);
		free(rmem);
		free(sigma);
		gsl_vector_free(B);
	}
}




















/**
 * Function LipRadiation
 *
 * Lip radiation
 *
 * @param excitation_voiced excitation
 *
 */
void LipRadiation(gsl_vector *excitation_voiced) {

	int i,j;
	double coeffs[2] = {1,-LIP_RADIATION};
	double sum;
	gsl_vector *temp = gsl_vector_alloc(excitation_voiced->size);
	gsl_vector_memcpy(temp, excitation_voiced);

	for(i=0; i<excitation_voiced->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, 1); j++) {
			sum += gsl_vector_get(excitation_voiced, i-j)*coeffs[j];
		}
		gsl_vector_set(temp, i, sum);
	}
	for(i=0; i<excitation_voiced->size; i++) {
		gsl_vector_set(excitation_voiced, i, gsl_vector_get(temp, i));
	}
	gsl_vector_free(temp);
}




/**
 * Function Differentiate
 *
 * Differentiate signal with leak
 *
 * @param signal signal to be differentiated
 * @param leak leaky parameter (e.g. 0.99)
 *
 */
void Differentiate(gsl_vector *signal, double leak) {

	int i,j;
	double coeffs[2] = {1,-leak};
	double sum;
	gsl_vector *temp = gsl_vector_alloc(signal->size);
	gsl_vector_memcpy(temp, signal);

	for(i=0; i<signal->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, 1); j++) {
			sum += gsl_vector_get(signal, i-j)*coeffs[j];
		}
		gsl_vector_set(temp, i, sum);
	}
	for(i=0; i<signal->size; i++) {
		gsl_vector_set(signal, i, gsl_vector_get(temp, i));
	}
	gsl_vector_free(temp);
}








/**
 * Function Differentiate_noleak
 *
 * Differentiate signal (no leakage)
 *
 * @param signal signal to be differentiated
 *
 */
void Differentiate_noleak(gsl_vector *signal) {

	int i,j;
	double coeffs[2] = {1,-1};
	double sum;
	gsl_vector *temp = gsl_vector_alloc(signal->size);
	gsl_vector_memcpy(temp, signal);

	for(i=0; i<signal->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, 1); j++) {
			sum += gsl_vector_get(signal, i-j)*coeffs[j];
		}
		gsl_vector_set(temp, i, sum);
	}
	for(i=0; i<signal->size; i++) {
		gsl_vector_set(signal, i, gsl_vector_get(temp, i));
	}
	gsl_vector_free(temp);
}



















/**
 * Function Evaluate_new_gain
 *
 * Evaluate new gain for synthesis by comparing synthetic and original gains
 *
 * @param gain original gain
 * @param gain_new new gain
 * @param params
 */
void Evaluate_new_gain(gsl_vector *signal, gsl_vector *gain_new, gsl_vector *gain, gsl_vector *fundf, PARAM *params) {

	/* Estimate the gain of the synthetic signal */
	Gain_eval(signal,gain_new,fundf,params);
	MA_voiced(gain_new,fundf,params->norm_gain_smooth_v_len);
	MA_unvoiced(gain_new,fundf,params->norm_gain_smooth_uv_len);

	/* Evaluate new gain for synthesis by comparing synthetic and original gains */
	Modify_gain(gain_new,gain);
}




/**
 * Function Gain_eval
 *
 * Evaluate gain of a signal
 *
 * @param signal input signal
 * @param gain gain vector
 * @param fundf F0 vector
 * @param params
 *
 */
void Gain_eval(gsl_vector *signal, gsl_vector *gain, gsl_vector *fundf, PARAM *params) {

	int i,j;
	gsl_vector *frame = gsl_vector_alloc(rint(params->frame_length/params->speed));
	int add = rint((params->frame_length/(double)params->shift/params->speed-1.0)*params->shift/params->speed/2.0);

	/* Zeropad signal */
	gsl_vector *signal_zp = gsl_vector_calloc(signal->size + 2*add);
	for(i=0;i<signal->size;i++)
		gsl_vector_set(signal_zp,i+add,gsl_vector_get(signal,i));

	/* Calculate gain vector */
	for(i=0;i<gain->size;i++) {
		for(j=rint(i*params->shift/params->speed);j<GSL_MIN(rint(i*params->shift/params->speed+params->frame_length/params->speed),signal_zp->size-1);j++) {
			gsl_vector_set(frame,GSL_MIN(j-rint(i*params->shift/params->speed),frame->size-1),gsl_vector_get(signal_zp,j));
		}
		if(gsl_vector_get(fundf,i) > 0)
			uvGain(frame,gain,i,NO_WINDOWING,rint(params->gain_voiced_frame_length/params->speed));
		else
			uvGain(frame,gain,i,NO_WINDOWING,rint(params->gain_unvoiced_frame_length/params->speed));
	}
	gsl_vector_free(frame);
	gsl_vector_free(signal_zp);
}







/**
 * Function Modify_gain
 *
 * Compare gain (gain) and synthetic gain (gain_new), and normalize the result inplace to gain_new
 *
 * @param gain_new synthetic gain
 * @param gain original gain
 *
 */
void Modify_gain(gsl_vector *gain_new, gsl_vector *gain) {

	int i;
	for(i=0;i<gain->size;i++) {
		gsl_vector_set(gain_new,i,gsl_vector_get(gain,i) + gsl_vector_get(gain,i) - gsl_vector_get(gain_new,i));
	}
}




/**
 * Function Noise_reduction
 *
 * Reduce noise by reducing the gain of low level parts of the signal
 *
 * @param gain vector
 * @param params
 */
void Noise_reduction(gsl_vector *gain, PARAM *params) {

	/* Apply noise reduction */
	if(params->noise_reduction_synthesis == 1) {
		int i;
		for(i=0;i<gain->size;i++)
			if(gsl_vector_get(gain,i) < params->noise_reduction_limit_db)
				gsl_vector_set(gain,i,gsl_vector_get(gain,i)-params->noise_reduction_db);
	}
}


/**
 * Function uvGain
 *
 * Calculate energy from potentially unvoiced frames (shorter window)
 *
 * @param frame pointer to the samples
 * @param gain vector for results
 * @param index time index
 * @param windowing switch for using windowing
 * @param unvoiced_frame_length
 */
void uvGain(gsl_vector *frame, gsl_vector *gain, int index, int windowing, int unvoiced_frame_length) {

	int i,cntr;
	double sum;
	gsl_vector *uvframe = gsl_vector_alloc(unvoiced_frame_length);

	/* Take shorter frame for unvoiced analysis */
	cntr = rint((frame->size-unvoiced_frame_length)/2.0);
	for(i=0;i<unvoiced_frame_length;i++)
		gsl_vector_set(uvframe,i,gsl_vector_get(frame,i+cntr));

	/* Windowing switch */
	if(windowing == 1) {

		/* Windowing */
		for(i=0;i<uvframe->size;i++)
			gsl_vector_set(uvframe,i,gsl_vector_get(uvframe,i)*HANN(i,uvframe->size));

		/* Evaluate gain of uvframe, normalize energy per sample basis */
		sum = 0;
		for(i=0;i<uvframe->size;i++) {
			sum = sum + gsl_vector_get(uvframe,i)*gsl_vector_get(uvframe,i);
		}
		gsl_vector_set(gain, index, 10.0*log10((8.0/3.0)*sum/E_REF/((double)(uvframe->size))));

	} else {

		/* Evaluate gain of frame, normalize energy per sample basis */
		sum = 0;
		for(i=0;i<uvframe->size;i++) {
			sum = sum + gsl_vector_get(uvframe,i)*gsl_vector_get(uvframe,i);
		}
		gsl_vector_set(gain, index, 10.0*log10(sum/E_REF/((double)(uvframe->size))));
	}

	/* Ensure non-infinity values */
	if(isinf(gsl_vector_get(gain,index)) != 0)
		gsl_vector_set(gain,index,MIN_LOG_POWER);

	/* Free memory */
	gsl_vector_free(uvframe);
}






























/**
 * Function Smooth_matrix
 *
 * Moving average smoothing for matrix parameters
 *
 * @param matrix matrix to smooth
 * @param len smoothing length in samples
 *
 */
void Smooth_matrix(gsl_matrix *matrix, int len) {

	if(matrix == NULL)
		return;

	int i,j;
	gsl_vector *temp = gsl_vector_alloc(matrix->size1);
	for(i=0;i<matrix->size2;i++) {
		for(j=0;j<matrix->size1;j++) {
			gsl_vector_set(temp,j,gsl_matrix_get(matrix,j,i));
		}
		MA(temp,len);
		for(j=0;j<matrix->size1;j++) {
			gsl_matrix_set(matrix,j,i,gsl_vector_get(temp,j));
		}
	}
	gsl_vector_free(temp);
}







/**
 * Function Spectral_match
 *
 * Match the synthetic excitation spectrum to real one (normal/warped)
 *
 * @param signal excitation vector
 * @param flow original glottal flow spectrum
 * @param flow_new synthetic glottal flow spectrum
 * @param
 *
 */
void Spectral_match(gsl_vector *signal, gsl_matrix *flow, gsl_matrix *flow_new, PARAM *params) {

	/* Do not perform if noise robust speech is used */
	if(params->noise_robust_speech == 1)
		return;

	/* Normal filtering */
	if(params->lambda_gl == 0) {

		int i,j;
		double sum;

		/* Check compatibility */
		if(params->lpc_order_gl < MIN_P_TILT) {
			printf("\n	Warning: The degree of spectral model is too low - spectral matching is not performed!\n\n");
			return;
		}

		/* Smooth and interpolate glottal flow spectra */
		gsl_matrix *flow_i = gsl_matrix_alloc(params->signal_length,params->lpc_order_gl);
		gsl_matrix *flow_new_i = gsl_matrix_alloc(params->signal_length,params->lpc_order_gl);
		Smooth_interp_lsf(flow_i,flow,params->signal_length,params->use_hmm,params->glflowsp_smooth_len);
		Smooth_interp_lsf(flow_new_i,flow_new,params->signal_length,0,params->glflowsp_smooth_len); // Smooth always flow_new

		/* Initialize spectral correction */
		gsl_vector *A = gsl_vector_alloc(params->lpc_order_gl+1);
		gsl_vector *B = gsl_vector_alloc(params->lpc_order_gl+1);
		gsl_vector *signal_orig = gsl_vector_alloc(params->signal_length);
		gsl_vector_memcpy(signal_orig, signal);

		/* Spectral correction */
		for(i=0;i<params->signal_length;i++) {

			/* Update filter coeffs: convert LSF to poly */
			if(i%params->filter_update_interval_gl == 0) {
		        lsf2poly(flow_new_i,A,i,params->use_hmm);
		        lsf2poly(flow_i,B,i,params->use_hmm);
				gsl_vector_set(B,0,0);
			}

			/* Filter */
	    	sum = 0;
	    	for(j=0;j<GSL_MIN(params->lpc_order_gl+1,i);j++) {
	    		sum += gsl_vector_get(signal_orig,i-j)*gsl_vector_get(A,j) - gsl_vector_get(signal,i-j)*gsl_vector_get(B,j);
	    	}
	    	gsl_vector_set(signal,i,sum);
		}

		/* Free memory */
		gsl_matrix_free(flow_new_i);
		gsl_matrix_free(flow_i);
		gsl_vector_free(signal_orig);
		gsl_vector_free(A);
		gsl_vector_free(B);


	/* Warped filtering */
	} else {

		/* Check compatibility */
		if(params->lpc_order_gl+1 < MIN_P_TILT-1) {
			printf("\n  Warning: The degree of spectral model is too low - spectral matching is not performed!\n\n");
			return;
		}

		int i,q,mlen;
	    long int o;
	    double xr,x,ffr,tmpr,Bb;
	    double *sigma;
	    long int len = params->signal_length;
	    int adim = params->lpc_order_gl + 1;
	    int bdim = params->lpc_order_gl + 1;
	    double *Ar = (double *)calloc(adim,sizeof(double));
	    double *Br = (double *)calloc(bdim,sizeof(double));
	    double *ynr = (double *)calloc(len,sizeof(double));
	    double *rsignal = (double *)calloc(len,sizeof(double));
	    double *rmem = (double *)calloc((bdim+2),sizeof(double));
	    gsl_vector *A = gsl_vector_alloc(adim);
	    gsl_vector *B = gsl_vector_alloc(bdim);

		/* Smooth and interpolate glottal flow spectrums */
		gsl_matrix *flow_i = gsl_matrix_alloc(params->signal_length,params->lpc_order_gl);
		gsl_matrix *flow_new_i = gsl_matrix_alloc(params->signal_length,params->lpc_order_gl);
		Smooth_interp_lsf(flow_i,flow,len,params->use_hmm,params->glflowsp_smooth_len);
		Smooth_interp_lsf(flow_new_i,flow_new,len,0,params->glflowsp_smooth_len); // Smooth always flow_new

	    /* Set signal to array */
	    for(i=0;i<len;i++) {
			rsignal[i] = gsl_vector_get(signal,i);
		}

	    /* Initialize */
	    sigma = NFArray(bdim+2);
	    Bb = 0;
	    if(adim >= bdim)
	    	mlen = adim;
	    else
	    	mlen = bdim + 1;

	    /* Warped filtering */
	    for(o=0;o<len;o++) {

	    	/* Update filter coefficients */
	    	if(o%params->filter_update_interval_gl == 0) {
	    		lsf2poly(flow_new_i,A,o,params->use_hmm);
	    		lsf2poly(flow_i,B,o,params->use_hmm);
	    		for(i=0;i<adim;i++) {
	    			Ar[i] = gsl_vector_get(A,i);
	    		}
	    		for(i=0;i<bdim;i++) {
	    			Br[i] = gsl_vector_get(B,i);
	    		}
	    		alphas2sigmas(Br,sigma,params->lambda_gl,bdim-1);
	    		Bb = 1/Br[0];
	    	}

	    	xr = rsignal[o]*Bb;

	    	/* Update feedbackward sum */
	    	for(q=0;q<bdim;q++) {
	    		xr -= sigma[q]*rmem[q];
	    	}
	    	xr = xr/sigma[bdim];
	    	x = xr*Ar[0];

	    	/* Update inner states */
	    	for(q=0;q<mlen;q++) {
	    		tmpr = rmem[q] + params->lambda_gl*(rmem[q+1] - xr);
	    		rmem[q] = xr;
	    		xr = tmpr;
	    	}

	    	/* Update feedforward sum */
	    	for(q=0,ffr=0.0;q<adim-1;q++) {
	    		ffr += Ar[q+1]*rmem[q+1];
	    	}

	       /* Update output */
	       ynr[o] = x + ffr;
		}

	    /* Set output to vector */
	    for(i=0;i<len;i++) {
	    	gsl_vector_set(signal,i,ynr[i]);
	    }

	    /* Free memory */
		free(ynr);
		free(rsignal);
		free(Ar);
		free(Br);
		free(rmem);
		free(sigma);
		gsl_matrix_free(flow_i);
		gsl_matrix_free(flow_new_i);
		gsl_vector_free(A);
		gsl_vector_free(B);
	}
}










/**
 * Function alphas2sigmas
 *
 * Convert alhas to sigmas.
 *
 * @param alp alphas
 * @param sigm sigmas
 * @param lambda warping coefficient
 * @param dim dimension
 *
 */
void alphas2sigmas(double *alp, double *sigm, double lambda, int dim) {

	int q;
	double S=0,Sp;

	sigm[dim] = lambda*alp[dim]/alp[0];
	Sp = alp[dim]/alp[0];
	for(q=dim;q>1;q--) {
		S = alp[q-1]/alp[0] - lambda*Sp;
		sigm[q-1] = lambda*S + Sp;
		Sp = S;
	}
	sigm[0] = S;
	sigm[dim+1] = 1 - lambda*S;
}


/**
 * Function NFArray
 *
 * Create array.
 *
 * @param size
 * @return array
 *
 */
double *NFArray(int size) {

	double *p;
	p = (double *)calloc(sizeof(*p),size);
	return p;
}



















/**
 * Function Scale_signal
 *
 * Scale signal if maximum value is greater than 1.0
 *
 * @param signal
 */
void Scale_signal(gsl_vector *signal, int mode) {

	/* Evaluate the absolute maximum of the signal */
	int i;
	double absmax = GSL_MAX(gsl_vector_max(signal),-gsl_vector_min(signal));

	/* Scale maximum to one if absmax is greater than one */
	if(mode == SCALE_IF_GREATER_THAN_ONE && absmax > 1.0) {
		printf("		Maximum value of the signal (%1.3lf) is greater than 1.0. Signal values rescaled!\n",absmax);
		absmax = WAV_SCALE/absmax;
		for(i=0;i<signal->size;i++)
			gsl_vector_set(signal,i,gsl_vector_get(signal,i)*absmax);
	}

	/* Scale absmax of signal to one */
	if(mode == FORCE_MAX_TO_ONE) {
		absmax = WAV_SCALE/absmax;
		for(i=0;i<signal->size;i++)
			gsl_vector_set(signal,i,gsl_vector_get(signal,i)*absmax);
	}
}





/**
 * Function Save_signal_to_file
 *
 * Save signal to file
 *
 * @param signal
 * @param params
 */
int Save_signal_to_file(gsl_vector *signal, PARAM *params, char *alternative_filename_ending) {

	/* Copy values to array */
	double *samples = (double *)calloc(params->signal_length, sizeof(double));
	int i;
	for(i=0;i<params->signal_length;i++)
		samples[i] = gsl_vector_get(signal,i);

	/* Save synthesized speech to wav-file */
	SNDFILE *soundfile;
	SF_INFO sfinfo;
	sfinfo.samplerate = params->FS;
	sfinfo.channels = 1;
	sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	char temp[DEF_STRING_LEN];
	strcpy(temp, params->synlist[params->synfilenumber]);

	/* Open file with default ending or givend ending */
	if(alternative_filename_ending == NULL)
		soundfile = sf_open(strcat(temp,FILENAME_ENDING_SYNTHESIS), SFM_WRITE, &sfinfo);
	else
		soundfile = sf_open(strcat(temp,alternative_filename_ending), SFM_WRITE, &sfinfo);

	/* Check if success */
	if(soundfile==NULL) {
		printf("\n\nError creating file \"%s\": %s\n\n",temp,strerror(errno));
		return EXIT_FAILURE;
	}

	/* Write to file */
	sf_write_double(soundfile, samples, params->signal_length);

	/* Free memory */
	free(samples);
	sf_close(soundfile);

	return EXIT_SUCCESS;
}








/**
 * Function WLPC
 *
 * Calculate Warped Linear Prediction (WLP) coefficients using
 * autocorrelation method.
 *
 * @param frame pointer to the samples
 * @param a pointer to coefficiets
 * @param p LPC degree
 * @param lambda warping coefficient
 * @return pointer to the WLP-coefficients
 */
gsl_vector *WLPC(gsl_vector *frame, int p, double lambda) {

	int i,j,s;
	double win = 0;
	gsl_vector *a = gsl_vector_alloc(p+1);
	gsl_vector *wframe = gsl_vector_alloc(frame->size);
	gsl_vector *a_temp = gsl_vector_calloc(p);
	gsl_vector *r = gsl_vector_calloc(p+1);
	gsl_vector *b = gsl_vector_alloc(p);
	gsl_matrix *R = gsl_matrix_alloc (p, p);
	gsl_permutation *perm = gsl_permutation_alloc(p);

	/* Windowing (choose window) */
	for(i=0; i<frame->size; i++) {
		if (WIN_TYPE == HANN_WIN) win = HANN(i,frame->size);
		else if (WIN_TYPE == BLACKMAN_WIN) win = BLACKMAN(i,frame->size);
		else if (WIN_TYPE == HAMMING_WIN) win = HAMMING(i,frame->size);
		else win = 1.0;
		gsl_vector_set(wframe, i, gsl_vector_get(frame, i)*win);
	}

	/* Copy warped frame */
	gsl_vector *wframe_w = gsl_vector_alloc(wframe->size);
	gsl_vector_memcpy(wframe_w,wframe);

	/* Set r(0) */
	for(i=0;i<wframe->size;i++) {
		gsl_vector_set(r, 0, gsl_vector_get(r,0) + gsl_vector_get(wframe,i)*gsl_vector_get(wframe,i));
	}

	/* Evaluate r */
	for(i=1;i<p+1;i++) {
		AllPassDelay(wframe_w,lambda);
		for(j=0;j<wframe->size;j++) {
			gsl_vector_set(r, i, gsl_vector_get(r,i) + gsl_vector_get(wframe,j)*gsl_vector_get(wframe_w,j));
		}
	}

	/* Autocorrelation matrix (Toeplitz) */
	for(i=0; i<p;i++) {
	    for(j=0; j<p; j++) {
	        gsl_matrix_set(R, i, j, gsl_vector_get(r, abs(i-j)));
	    }
	}

	/* Vector b */
	for(i=1; i<p+1; i++) {
		gsl_vector_set(b, i-1, gsl_vector_get(r, i));
	}

	/* Ra=r solver (LU-decomposition) */
	gsl_linalg_LU_decomp(R, perm, &s);
	gsl_linalg_LU_solve(R, perm, b, a_temp);

	/* Construct vector a and return */
	for(i=1;i<p+1;i++) {
		gsl_vector_set(a, i, -1.0*gsl_vector_get(a_temp,i-1));
	}
	gsl_vector_set(a,0,1);

	/* Replace NaN-values with zeros in case of all-zero frames */
	for(i=0;i<a->size;i++) {
		if(gsl_isnan(gsl_vector_get(a,i)))
			gsl_vector_set(a,i,0);
	}

	/* Free memory */
	gsl_vector_free(wframe);
	gsl_vector_free(wframe_w);
	gsl_vector_free(a_temp);
	gsl_vector_free(r);
	gsl_vector_free(b);
	gsl_matrix_free(R);
	gsl_permutation_free(perm);

	return a;
}





/**
 * Function AllPassDelay
 *
 * All pass delay filter for WLPC.
 *
 * @param signal pointer to the samples
 * @param lambda all-pass filter coefficient
 */
void AllPassDelay(gsl_vector *signal, double lambda) {

	int i,j,n=2;
	double sum;

	/* Create coefficient arrays: B = [-lambda 1], A = [1 -lambda] */
	double B[2] = {-lambda,1.0};
	double A[2] = {0,lambda};

	/* Zeropad the beginning of the signal */
	gsl_vector *signal_zp = gsl_vector_calloc(signal->size+n);
	for(i=n;i<signal_zp->size;i++) {
		gsl_vector_set(signal_zp,i,gsl_vector_get(signal,i-n));
	}

	/* Copy vector signal */
	gsl_vector *signal_zp_orig = gsl_vector_alloc(signal_zp->size);
	gsl_vector_memcpy(signal_zp_orig,signal_zp);

	/* Filter */
	for(i=n;i<signal_zp->size;i++) {
    	sum = 0;
    	for(j=0;j<n;j++) {
    		sum += gsl_vector_get(signal_zp_orig,i-j)*B[j] + gsl_vector_get(signal_zp,i-j)*A[j];
    	}
    	gsl_vector_set(signal_zp,i,sum);
	}

	/* Remove zeropadding from the signal */
	for(i=n;i<signal_zp->size;i++) {
		gsl_vector_set(signal,i-n,gsl_vector_get(signal_zp,i));
	}

	/* Free memory */
	gsl_vector_free(signal_zp);
	gsl_vector_free(signal_zp_orig);
}





/**
 * Function SWLP
 *
 * Calculate Stabilized Weighted Linear Prediction (SWLP) coefficients
 * (unstabilized can be evaluated by switching "stabilized" to zero)
 *
 * @param frame pointer to the samples
 * @param a pointer to coefficiets
 * @param p LPC degree
 * @param M weighting window length
 * @param lag lag of the weighting window
 */
void SWLP(gsl_vector *frame, gsl_vector *a, int M, int lag, int weighting, gsl_vector *fundf, int FS, int index, gsl_vector *glottsig, int stabilized) {

	int i,j,k,s,p = a->size-1;
	double win = 0,sum = 0;
	gsl_vector *wframe = gsl_vector_alloc(frame->size);
	gsl_vector *weight = gsl_vector_calloc(frame->size+p);

	/* Windowing (choose window) */
	for(i=0; i<frame->size; i++) {
		if (WIN_TYPE == HANN_WIN) win = HANN(i,frame->size);
		else if (WIN_TYPE == BLACKMAN_WIN) win = BLACKMAN(i,frame->size);
		else if (WIN_TYPE == HAMMING_WIN) win = HAMMING(i,frame->size);
		else win = 1.0;
		gsl_vector_set(wframe, i, gsl_vector_get(frame, i)*win);
	}

	/* Evaluate weighting function */
	if(weighting == 0)
		Eval_STE_weight(wframe,weight,M,lag);
	else
		Eval_GCI_weight(glottsig,weight,fundf,FS,index);

	/* Create partial weights */
	gsl_matrix *Z = gsl_matrix_calloc(wframe->size+p,p+1); // Partial weights
	gsl_matrix *Y = gsl_matrix_calloc(wframe->size+p,p+1); // Delayed and weighted versions of the signal
	for(i=0;i<weight->size;i++)
		gsl_matrix_set(Z,i,0,sqrt(gsl_vector_get(weight,i)));
	for(i=0;i<wframe->size;i++)
		gsl_matrix_set(Y,i,0,gsl_vector_get(wframe,i)*sqrt(gsl_vector_get(weight,i)));
	for(i=0;i<p;i++) {
		if(stabilized == 1) {
			for(j=i+1;j<Z->size1;j++) {
				gsl_matrix_set(Z,j,i+1,GSL_MAX(sqrt(gsl_vector_get(weight,j)/gsl_vector_get(weight,j-1)),1)*gsl_matrix_get(Z,j-1,i));
			}
		} else {
			for(j=i+1;j<Z->size1;j++) {
				gsl_matrix_set(Z,j,i+1,sqrt(gsl_vector_get(weight,j)));
			}
		}
		for(j=0;j<wframe->size;j++) {
			gsl_matrix_set(Y,j+i+1,i+1,gsl_matrix_get(Z,j+i+1,i+1)*gsl_vector_get(wframe,j));
		}
	}

	/* Autocorrelation matrix R (R = (YT*Y)/N, size p*p) and vector b (size p) */
	gsl_matrix *R = gsl_matrix_calloc(p,p);
	gsl_vector *b = gsl_vector_calloc(p);
	for(i=0;i<Y->size2;i++) {
		for(j=0;j<Y->size2;j++) {
			if(i > 0 && j > 0) {
				for(k=0;k<Y->size1;k++) {
					gsl_matrix_set(R,i-1,j-1,gsl_matrix_get(R,i-1,j-1) + gsl_matrix_get(Y,k,i)*gsl_matrix_get(Y,k,j));
				}
				gsl_matrix_set(R,i-1,j-1,gsl_matrix_get(R,i-1,j-1)/wframe->size);
			}
			if(i > 0 && j == 0) {
				for(k=0;k<Y->size1;k++)
					gsl_vector_set(b,i-1,gsl_vector_get(b,i-1) + gsl_matrix_get(Y,k,i)*gsl_matrix_get(Y,k,j));
				gsl_vector_set(b,i-1,gsl_vector_get(b,i-1)/wframe->size);
				sum += gsl_vector_get(b,i-1);
			}
		}
	}

	/* Ra=r solver (LU-decomposition) (Do not evaluate LU if sum = 0) */
	gsl_vector *a_temp = gsl_vector_calloc(p);
	gsl_permutation *perm = gsl_permutation_alloc(p);
	if(sum != 0) {
		gsl_linalg_LU_decomp(R, perm, &s);
		gsl_linalg_LU_solve(R, perm, b, a_temp);
	}

	/* Set LP-coefficients to vector "a" */
	for(i=1; i<a->size; i++) {
		gsl_vector_set(a, i, (-1)*gsl_vector_get(a_temp, i-1));
	}
	gsl_vector_set(a, 0, 1);

	/* Stabilize through LSFs if the method itself is not guaranteed to produce a stable filter */
	if(stabilized == 0)
		LSF_stabilize(a);

	/* Free memory */
	gsl_vector_free(weight);
	gsl_vector_free(wframe);
	gsl_matrix_free(Z);
	gsl_matrix_free(Y);
	gsl_matrix_free(R);
	gsl_vector_free(b);
	gsl_vector_free(a_temp);
	gsl_permutation_free(perm);
}







/**
 * Function LSF_stabilize
 *
 * Check the validity of polynomial through LSFs and fix found errors
 *
 * @param lsf vector
 */
void LSF_stabilize(gsl_vector *a) {

	gsl_vector *lsf = gsl_vector_calloc(a->size-1);
	Convert_vector_to_LSF(a, lsf);
	LSF_fix_vector(lsf);
	lsf_vector2poly(lsf,a);
}




/**
 * Function Convert_vector_to_LSF
 *
 * Convert LPC-coefficients to Line Spectrum Frequencies (LSF)
 * Maximum LPC-polynomial degree is 36!
 *
 *
 * @param a pointer to LPC-vector
 * @param LSF pointer to the LSF matrix
 *
 */
void Convert_vector_to_LSF(gsl_vector *a, gsl_vector *LSF) {

	int i,n;

	/* Count the number of nonzero elements in "a" */
	n = 0;
	for(i=0; i<a->size; i++) {
		if(gsl_vector_get(a, i) != 0) {
			n++;
		}
	}

	/* In case of only one non-zero element */
	if(n == 1) {
		for(i=0; i<LSF->size; i++)
			gsl_vector_set(LSF, i, (i+1)*M_PI/LSF->size);
		return;
	}

	gsl_vector *aa = gsl_vector_calloc(n+1);
	gsl_vector *flip_aa = gsl_vector_calloc(n+1);
	gsl_vector *p = gsl_vector_calloc(n+1);
	gsl_vector *q = gsl_vector_calloc(n+1);

	/* Construct vectors aa=[a 0] and flip_aa=[0 flip(a)] */
	for(i=0; i<n; i++) {
		gsl_vector_set(aa, i, gsl_vector_get(a, i));
		gsl_vector_set(flip_aa, flip_aa->size-1-i, gsl_vector_get(a, i));
	}

	/* Construct vectors p and q */
	for(i=0; i<n+1; i++) {
		gsl_vector_set(p, i, gsl_vector_get(aa, i) + gsl_vector_get(flip_aa, i));
		gsl_vector_set(q, i, gsl_vector_get(aa, i) - gsl_vector_get(flip_aa, i));
	}

	/* Remove trivial zeros */
	if((n+1)%2 == 0) {
		double y;
		/* Deconvolve p with [1 1] */
		y = 0;
		for(i=0; i<p->size; i++) {
    		gsl_vector_set(p, i, gsl_vector_get(p, i)-y);
    		y = gsl_vector_get(p, i);
		}
		gsl_vector_set(p, p->size-1, 0);
		/* Deconvolve q with [1 -1] */
		y = 0;
		for(i=0; i<q->size; i++) {
    		gsl_vector_set(q, i, gsl_vector_get(q, i)+y);
    		y = gsl_vector_get(q, i);
		}
		gsl_vector_set(q, q->size-1, 0);
	} else {
		double y;
		/* Deconvolve q with [1 1] */
		y = 0;
		for(i=0; i<q->size; i++) {
    		gsl_vector_set(q, i, gsl_vector_get(q, i)-y);
    		y = gsl_vector_get(q, i);
		}
		/* Deconvolve q with [1 -1] */
		y = 0;
		for(i=0; i<q->size; i++) {
    		gsl_vector_set(q, i, gsl_vector_get(q, i)+y);
    		y = gsl_vector_get(q, i);
		}
		gsl_vector_set(q, q->size-1, 0);
		gsl_vector_set(q, q->size-2, 0);
	}

	/* Count the number of nonzero elements in "p" and "q" */
	int n_p = 0;
	int n_q = 0;
	for(i=0; i<p->size; i++) {
		if(gsl_vector_get(p, i) != 0)
			n_p++;
		if(gsl_vector_get(q, i) != 0)
			n_q++;
	}

	/* Take the last half of "p" and "q" to vectors "TP" and "TQ" */
	gsl_vector *TP = gsl_vector_alloc((n_p+1)/2);
	gsl_vector *TQ = gsl_vector_alloc((n_q+1)/2);
	for(i=0; i<TP->size; i++) {
		gsl_vector_set(TP, i, gsl_vector_get(p, i+(n_p-1)/2));
	}
	for(i=0; i<TQ->size; i++) {
		gsl_vector_set(TQ, i, gsl_vector_get(q, i+(n_q-1)/2));
	}

	/* Chebyshev transform */
	Chebyshev(TP);
	Chebyshev(TQ);

	/* Initialize root arrays */
	int nroots_TP = 2*(TP->size-1);
	int nroots_TQ = 2*(TQ->size-1);
	double TP_coeffs[TP->size];
	double TQ_coeffs[TQ->size];
	double TP_roots[nroots_TP];
	double TQ_roots[nroots_TQ];

	/* Copy coefficients to arrays */
    for(i=0; i<TP->size; i++) {
    	TP_coeffs[i] = gsl_vector_get(TP, i);
    }
    for(i=0; i<TQ->size; i++) {
    	TQ_coeffs[i] = gsl_vector_get(TQ, i);
    }

	/* Solve roots */
    gsl_poly_complex_workspace *w_TP = gsl_poly_complex_workspace_alloc(TP->size);
    gsl_poly_complex_workspace *w_TQ = gsl_poly_complex_workspace_alloc(TQ->size);
    gsl_poly_complex_solve(TP_coeffs, TP->size, w_TP, TP_roots);
    gsl_poly_complex_solve(TQ_coeffs, TQ->size, w_TQ, TQ_roots);

	/* Convert to LSF and sort */
	double sorted_LSF[(nroots_TP+nroots_TQ)/2];
	for(i=0; i<nroots_TP; i=i+2) {
		sorted_LSF[i/2] = acos(TP_roots[i]/2);
	}
	for(i=0; i<nroots_TQ; i=i+2) {
		sorted_LSF[(i+nroots_TP)/2] = acos(TQ_roots[i]/2);
	}
	gsl_sort(sorted_LSF, 1, (nroots_TP+nroots_TQ)/2);

	/* Copy LSFs to vector LSF */
	for(i=0; i<(nroots_TP+nroots_TQ)/2; i++)
		gsl_vector_set(LSF, i, sorted_LSF[i]);

	/* Free memory */
	gsl_poly_complex_workspace_free(w_TP);
    gsl_poly_complex_workspace_free(w_TQ);
	gsl_vector_free(aa);
	gsl_vector_free(flip_aa);
	gsl_vector_free(p);
	gsl_vector_free(q);
	gsl_vector_free(TP);
	gsl_vector_free(TQ);
}






/**
 * Function lsf_vector2poly
 *
 * Convert LSF vector to polynomial
 *
 * @param lsf_vector
 * @param poly
 * @param index
 * @param HMM switch for HMM
 */
void lsf_vector2poly(gsl_vector *lsf_vector, gsl_vector *poly) {

	int i,l = lsf_vector->size;
	gsl_vector *fi_p = NULL, *fi_q = NULL;

	/* Create fi_p and fi_q */
	if(l%2 == 0) {
		fi_p = gsl_vector_alloc(l/2);
		fi_q = gsl_vector_alloc(l/2);
		for(i=0;i<l;i=i+2) {
			gsl_vector_set(fi_p,i/2,gsl_vector_get(lsf_vector,i));
		}
		for(i=1;i<l;i=i+2) {
			gsl_vector_set(fi_q,(i-1)/2,gsl_vector_get(lsf_vector,i));
		}
	} else {
		fi_p = gsl_vector_alloc((l+1)/2);
		for(i=0;i<l;i=i+2) {
			gsl_vector_set(fi_p,i/2,gsl_vector_get(lsf_vector,i));
		}
		if((l-1)/2 > 0) {
			fi_q = gsl_vector_alloc((l-1)/2);
			for(i=1;i<l-1;i=i+2) {
				gsl_vector_set(fi_q,(i-1)/2,gsl_vector_get(lsf_vector,i));
			}
		}
	}

	/* Construct vectors P and Q */
	gsl_vector *cp = gsl_vector_calloc(3);
	gsl_vector *cq = gsl_vector_calloc(3);
	gsl_vector_add_constant(cp,1);
	gsl_vector_add_constant(cq,1);
	gsl_vector *P = gsl_vector_alloc(1);
	gsl_vector *Q = gsl_vector_alloc(1);
	gsl_vector_set(P,0,1);
	gsl_vector_set(Q,0,1);
	for(i=0;i<fi_p->size;i++) {
		gsl_vector_set(cp,1,-2*cos(gsl_vector_get(fi_p,i)));
		P = Conv(P,cp);
	}
	if((l-1)/2 > 0) {
		for(i=0;i<fi_q->size;i++) {
			gsl_vector_set(cq,1,-2*cos(gsl_vector_get(fi_q,i)));
			Q = Conv(Q,cq);
		}
	}

	/* Add trivial zeros */
	if(l%2 == 0) {
		gsl_vector *conv = gsl_vector_calloc(2);
		gsl_vector_add_constant(conv,1);
		P = Conv(P,conv);
		gsl_vector_set(conv,0,-1);
    	Q = Conv(Q,conv);
    	gsl_vector_free(conv);
	} else {
		gsl_vector *conv = gsl_vector_calloc(3);
		gsl_vector_set(conv,0,-1);
		gsl_vector_set(conv,2,1);
		Q = Conv(Q,conv);
		gsl_vector_free(conv);
	}

	/* Construct polynomial */
	for(i=1;i<P->size;i++) {
		gsl_vector_set(poly,P->size-i-1,0.5*(gsl_vector_get(P,i)+gsl_vector_get(Q,i)));
	}

	/* Free memory */
	gsl_vector_free(fi_p);
	if((l-1)/2 > 0)
		gsl_vector_free(fi_q);
	gsl_vector_free(cp);
	gsl_vector_free(cq);
	gsl_vector_free(P);
	gsl_vector_free(Q);
}









/**
 * Function Eval_STE_weight
 *
 * Evaluate short time energy (STE) function for SWLP weighting
 *
 * @param frame pointer to the frame
 * @param ste pointer to the STE function
 * @param M weighting window length
 * @param lag time lag of estimating the STE
 */
void Eval_STE_weight(gsl_vector *frame, gsl_vector *ste, int M, int lag) {

	int i,j;
	for(i=0;i<ste->size;i++) {
		for(j=GSL_MAX(i-lag-M+1,0);j<GSL_MIN(i-lag+1,(int)frame->size);j++) {
			gsl_vector_set(ste,i,gsl_vector_get(ste,i) + gsl_vector_get(frame,j)*gsl_vector_get(frame,j));
		}
		gsl_vector_set(ste,i,gsl_vector_get(ste,i) + DBL_EPSILON);
	}
}




/**
 * Function Eval_GCI_weight
 *
 * Find GCIs (Glottal Closure Instants) and construct weighting for de-emphasizing GCIs in (S)WLP
 *
 * @param ...
 *
 */
void Eval_GCI_weight(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index) {

	/* Estimate glottal closure instants */
	gsl_vector *inds = Find_GCI(glottsig, fundf, FS, index);

	/* If unvoiced or GCIs not found, simply set weight to 1 */
	if(inds == NULL || gsl_vector_get(fundf,index) == 0) {
		if(inds != NULL)
			gsl_vector_free(inds);
		gsl_vector_set_all(weight,1);
		return;
	}

	/* Algorithm parameters */
	double closed_time = 0.4;
	double gci_pos = 0.8;
	double deemph = 0.01;

	/* Initialize */
	int i,j;
	int csamples = rint(closed_time*FS/gsl_vector_get(fundf,index));

	/* Set weight according to GCIs */
	gsl_vector_set_all(weight,1);
	for(i=0;i<inds->size;i++) {
		for(j=0;j<csamples;j++) {
			gsl_vector_set(weight,GSL_MIN(GSL_MAX(gsl_vector_get(inds,i)-rint(csamples*gci_pos)+j,0),weight->size-1),deemph);
		}
	}

	/* Free memory */
	gsl_vector_free(inds);
}







/**
 * Function Find_GCI
 *
 * Find GCIs (Glottal Closure Instants)
 *
 * @param ...
 *
 */
gsl_vector *Find_GCI(gsl_vector *frame_orig, gsl_vector *fundf, int FS, int index) {

	int i,j,min_ind,ind_ind,t0_tmp;
	double t0,min_val,temp;
	gsl_vector *indices = gsl_vector_calloc(100);
	gsl_vector *final_inds;
	gsl_vector *frame = gsl_vector_calloc(frame_orig->size);
	gsl_vector_memcpy(frame,frame_orig);

	/* Differentiate frame and find minimum */
	Differentiate(frame,LEAK);
	min_ind = gsl_vector_min_index(frame);

	/* Find t0 minima of the glottal waveform in order to find GCIs,
	 * start evaluating from the mimina, first backward, then forward */
	t0 = FS/gsl_vector_get(fundf,index);
	gsl_vector_set(indices,50,min_ind);
	t0_tmp = min_ind;
	ind_ind = 51;

	/* Backward */
	temp = 0;
	ind_ind = 51;
	while(1) {
		t0_tmp = rint(gsl_vector_get(indices,ind_ind-1) - t0 - temp);
		if(t0_tmp < 0)
			break;
		min_val = BIG_POS_NUMBER;
		for(i=-20;i<21;i++) {
			if(gsl_vector_get(frame,GSL_MIN(GSL_MAX(t0_tmp + i,0),frame->size-1)) < min_val) {
				min_val = gsl_vector_get(frame,GSL_MIN(GSL_MAX(t0_tmp + i,0),frame->size-1));
				gsl_vector_set(indices,ind_ind,GSL_MIN(GSL_MAX(t0_tmp + i,0),frame->size-1));
			}
		}
		if(gsl_vector_get(indices,ind_ind-1)-gsl_vector_get(indices,ind_ind) < t0/2.0) {
			temp = temp + t0;
		} else {
			ind_ind++;
			temp = 0;
		}
	}

	/* Forward */
	temp = 0;
	ind_ind = 49;
	while(1) {
		t0_tmp = rint(gsl_vector_get(indices,ind_ind+1) + t0 + temp);
		if(t0_tmp > frame->size-1)
			break;
		min_val = BIG_POS_NUMBER;
		for(i=-20;i<21;i++) {
			if(gsl_vector_get(frame,GSL_MIN(GSL_MAX(t0_tmp + i,0),frame->size-1)) < min_val) {
				min_val = gsl_vector_get(frame,GSL_MIN(GSL_MAX(t0_tmp + i,0),frame->size-1));
				gsl_vector_set(indices,ind_ind,GSL_MIN(GSL_MAX(t0_tmp + i,0),frame->size-1));
			}
		}
		if(gsl_vector_get(indices,ind_ind)-gsl_vector_get(indices,ind_ind+1) < t0/2.0) {
			temp = temp + t0;
		} else {
			ind_ind--;
			temp = 0;
		}
	}

	/* Sort indices */
	gsl_sort_vector(indices);

	/* Allocate vector for non-zero indices and return */
	i = 0;
	while(gsl_vector_get(indices,indices->size-1-i) > 0)
		i += 1;
	if(i < 2) {
		gsl_vector_free(indices);
		gsl_vector_free(frame);
		return NULL;
	} else {
		final_inds = gsl_vector_alloc(i);
		for(j=0;j<i;j++)
			gsl_vector_set(final_inds,j,gsl_vector_get(indices,indices->size-i+j));
		gsl_vector_free(indices);
		gsl_vector_free(frame);
		return final_inds;
	}
}






























/**
 * Function Interpolate
 *
 * Interpolates given vector to new vector of given length
 *
 * @param vector original vector
 * @param i_vector interpolated vector
 */
void Interpolate(gsl_vector *vector, gsl_vector *i_vector) {

	int i,len = vector->size,length = i_vector->size;

	/* Read values to array */
	double x[len];
	double y[len];
	for(i=0; i<len; i++) {
		x[i] = i;
		y[i] = gsl_vector_get(vector,i);
	}
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,len);
    gsl_spline_init(spline, x, y, len);
    double xi;
    i = 0;

    /* New implementation (27.3.2009, bug fix 8.2.2010) */
    /* Bug fix to GSL v.1.15, 26.1.2012 */
    xi = x[0];
    while(i<length) {
    	gsl_vector_set(i_vector,i,gsl_spline_eval(spline, xi, acc));
    	xi += (len-1)/(double)(length-1);
    	if(xi > len-1)
    		xi = len-1;
    	i++;
    }

    /* Free memory */
    gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}

/**
 * Function Interpolate_matrix
 *
 * Interpolates given matrix to new matrix of given length
 *
 * @param matrix original matrix
 * @param imatrix interpolated matrix
 */
void Interpolate_matrix(gsl_matrix *matrix, gsl_matrix *imatrix) {

	int i,j;
	gsl_vector *ivec = gsl_vector_alloc(imatrix->size1);
	for(i=0;i<matrix->size2;i++) {
		gsl_vector_view column = gsl_matrix_column(matrix,i);
		Interpolate((gsl_vector *)(&column),ivec);
		for(j=0;j<ivec->size;j++)
			gsl_matrix_set(imatrix,j,i,gsl_vector_get(ivec,j));
	}
	gsl_vector_free(ivec);
}




/**
 * Function Interpolate_lin
 *
 * Interpolates linearly given vector to new vector of given length
 *
 * @param vector original vector
 * @param i_vector interpolated vector
 */
void Interpolate_lin(gsl_vector *vector, gsl_vector *i_vector) {

	int i,len = vector->size,length = i_vector->size;

	/* Read values to array */
	double x[len];
	double y[len];
	for(i=0; i<len; i++) {
		x[i] = i;
		y[i] = gsl_vector_get(vector,i);
	}
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, len);
	gsl_spline_init(spline, x, y, len);
	double xi;
    i = 0;

    /* New implementation (27.3.2009, bug fix 8.2.2010) */
    /* Bug fix to GSL v.1.15, 26.1.2012 */
    xi = x[0];
    while(i<length) {
    	gsl_vector_set(i_vector,i,gsl_spline_eval(spline, xi, acc));
    	xi += (len-1)/(double)(length-1);
    	if(xi > len-1)
    		xi = len-1;
    	i++;
    }

    /* Free memory */
    gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}




/**
 * Function Interpolate_fract
 *
 * Interpolates given vector to new vector of given fractional length
 *
 * @param vector original vector
 * @param i_vector interpolated vector
 * @param flen fractional length
 */
void Interpolate_fract(gsl_vector *vector, gsl_vector *i_vector, double flen) {

	int i,len = vector->size,length = i_vector->size;

	/* Read values to array */
	double x[len];
	double y[len];
	for(i=0; i<len; i++) {
		x[i] = i;
		y[i] = gsl_vector_get(vector,i);
	}
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,len);
    gsl_spline_init(spline, x, y, len);
    double xi = x[0];
    i = 0;
    while(i<length) {
    	gsl_vector_set(i_vector,i,gsl_spline_eval(spline, xi, acc));
    	xi += (len-1)/(flen-1.0);
    	if(xi > len-1)
    		xi = len-1;
    	i++;
    }

    /* Free memory */
    gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}











/**
 * Function Convert_to_LSF
 *
 * Converts the LPC-coefficients to Line Spectrum Frequencies (LSF)
 * Maximum LPC-polynomial degree is 36!
 *
 * @param LSF pointer to the LSF matrix
 * @param a pointer to LPC-vector
 * @param index time index
 *
 */
void Convert_to_LSF(gsl_matrix *LSF, gsl_vector *a, int index) {

	int i,n;

	/* Count the number of nonzero elements in "a" */
	n = 0;
	for(i=0; i<a->size; i++) {
		if(gsl_vector_get(a, i) != 0) {
			n++;
		}
	}

	/* In case of only one non-zero element */
	if(n == 1) {
		for(i=0; i<LSF->size2; i++)
			gsl_matrix_set(LSF, index, i, (i+1)*M_PI/LSF->size2);
		return;
	}

	gsl_vector *aa = gsl_vector_calloc(n+1);
	gsl_vector *flip_aa = gsl_vector_calloc(n+1);
	gsl_vector *p = gsl_vector_calloc(n+1);
	gsl_vector *q = gsl_vector_calloc(n+1);

	/* Construct vectors aa=[a 0] and flip_aa=[0 flip(a)] */
	for(i=0; i<n; i++) {
		gsl_vector_set(aa, i, gsl_vector_get(a, i));
		gsl_vector_set(flip_aa, flip_aa->size-1-i, gsl_vector_get(a, i));
	}

	/* Construct vectors p and q */
	for(i=0; i<n+1; i++) {
		gsl_vector_set(p, i, gsl_vector_get(aa, i) + gsl_vector_get(flip_aa, i));
		gsl_vector_set(q, i, gsl_vector_get(aa, i) - gsl_vector_get(flip_aa, i));
	}

	/* Remove trivial zeros */
	if((n+1)%2 == 0) {
		double y;
		/* Deconvolve p with [1 1] */
		y = 0;
		for(i=0; i<p->size; i++) {
    		gsl_vector_set(p, i, gsl_vector_get(p, i)-y);
    		y = gsl_vector_get(p, i);
		}
		gsl_vector_set(p, p->size-1, 0);
		/* Deconvolve q with [1 -1] */
		y = 0;
		for(i=0; i<q->size; i++) {
    		gsl_vector_set(q, i, gsl_vector_get(q, i)+y);
    		y = gsl_vector_get(q, i);
		}
		gsl_vector_set(q, q->size-1, 0);
	} else {
		double y;
		/* Deconvolve q with [1 1] */
		y = 0;
		for(i=0; i<q->size; i++) {
    		gsl_vector_set(q, i, gsl_vector_get(q, i)-y);
    		y = gsl_vector_get(q, i);
		}
		/* Deconvolve q with [1 -1] */
		y = 0;
		for(i=0; i<q->size; i++) {
    		gsl_vector_set(q, i, gsl_vector_get(q, i)+y);
    		y = gsl_vector_get(q, i);
		}
		gsl_vector_set(q, q->size-1, 0);
		gsl_vector_set(q, q->size-2, 0);
	}

	/* Count the number of nonzero elements in "p" and "q" */
	int n_p = 0;
	int n_q = 0;
	for(i=0; i<p->size; i++) {
		if(gsl_vector_get(p, i) != 0)
			n_p++;
		if(gsl_vector_get(q, i) != 0)
			n_q++;
	}

	/* Take the last half of "p" and "q" to vectors "TP" and "TQ" */
	gsl_vector *TP = gsl_vector_alloc((n_p+1)/2);
	gsl_vector *TQ = gsl_vector_alloc((n_q+1)/2);
	for(i=0; i<TP->size; i++) {
		gsl_vector_set(TP, i, gsl_vector_get(p, i+(n_p-1)/2));
	}
	for(i=0; i<TQ->size; i++) {
		gsl_vector_set(TQ, i, gsl_vector_get(q, i+(n_q-1)/2));
	}

	/* Chebyshev transform */
	Chebyshev(TP);
	Chebyshev(TQ);

	/* Initialize root arrays */
	int nroots_TP = 2*(TP->size-1);
	int nroots_TQ = 2*(TQ->size-1);
	double TP_coeffs[TP->size];
	double TQ_coeffs[TQ->size];
	double TP_roots[nroots_TP];
	double TQ_roots[nroots_TQ];

	/* Copy coefficients to arrays */
    for(i=0; i<TP->size; i++) {
    	TP_coeffs[i] = gsl_vector_get(TP, i);
    }
    for(i=0; i<TQ->size; i++) {
    	TQ_coeffs[i] = gsl_vector_get(TQ, i);
    }

	/* Solve roots */
    gsl_poly_complex_workspace *w_TP = gsl_poly_complex_workspace_alloc(TP->size);
    gsl_poly_complex_workspace *w_TQ = gsl_poly_complex_workspace_alloc(TQ->size);
    gsl_poly_complex_solve(TP_coeffs, TP->size, w_TP, TP_roots);
    if(nroots_TQ > 0) {
    	gsl_poly_complex_solve(TQ_coeffs, TQ->size, w_TQ, TQ_roots);
    }

	/* Convert to LSF and sort */
	double sorted_LSF[(nroots_TP+nroots_TQ)/2];
	for(i=0; i<nroots_TP; i=i+2) {
		sorted_LSF[i/2] = acos(TP_roots[i]/2);
	}
	for(i=0; i<nroots_TQ; i=i+2) {
		sorted_LSF[(i+nroots_TP)/2] = acos(TQ_roots[i]/2);
	}
	gsl_sort(sorted_LSF, 1, (nroots_TP+nroots_TQ)/2);

	/* Copy LSFs to vector LSF */
	for(i=0; i<(nroots_TP+nroots_TQ)/2; i++)
		gsl_matrix_set(LSF, index, i, sorted_LSF[i]);

	/* Free memory */
	gsl_poly_complex_workspace_free(w_TP);
    gsl_poly_complex_workspace_free(w_TQ);
	gsl_vector_free(aa);
	gsl_vector_free(flip_aa);
	gsl_vector_free(p);
	gsl_vector_free(q);
	gsl_vector_free(TP);
	gsl_vector_free(TQ);
}









/**
 * Function Chebyshev
 *
 * Chebyshev transformation
 *
 * @param T pointer to polynomial coefficients (result will be computed in-place)
 *
 */
void Chebyshev(gsl_vector *T) {

	int i,j,s;
	gsl_matrix *C;
	gsl_vector *cheb = gsl_vector_alloc(T->size);
	gsl_matrix *inv_C = gsl_matrix_alloc(T->size, T->size);
	gsl_permutation *perm = gsl_permutation_alloc(T->size);

	/* Create C, matrix inversion through LU-decomposition */
	C = Construct_C(T->size);
	gsl_linalg_LU_decomp(C, perm, &s);
	gsl_linalg_LU_invert(C, perm, inv_C);

	/* Evaluate r = inv(C)'*T */
	double sum;
	for(i=0; i<T->size; i++) {
		sum = 0;
		for(j=0; j<T->size; j++)
			sum += gsl_matrix_get(inv_C, j, i)*gsl_vector_get(T, j);
		gsl_vector_set(cheb, i, sum);
	}

	/* Copy coefficients to T */
	for(i=0; i<T->size; i++) {
		gsl_vector_set(T, i, gsl_vector_get(cheb, i));
	}

	/* Free memory */
	gsl_vector_free(cheb);
	gsl_matrix_free(C);
	gsl_matrix_free(inv_C);
	gsl_permutation_free(perm);
}





/**
 * Function Construct_C
 *
 * Construct data matrix C for Chebyshev transform
 *
 * @param size
 * @return C matrix
 *
 */
gsl_matrix *Construct_C(int size) {

	int i,j;
	gsl_matrix *C = gsl_matrix_calloc(size,size);
	gsl_vector *tmp = gsl_vector_alloc(3);
	gsl_vector *f = gsl_vector_alloc(3);

	/* Set tmp and f */
	gsl_vector_set(tmp,0,1);
	gsl_vector_set(tmp,1,0);
	gsl_vector_set(tmp,2,1);
	gsl_vector_memcpy(f,tmp);

	/* Set diagonal to 1 */
	for(i=0;i<size;i++)
		gsl_matrix_set(C,i,i,1);

	/* Construct C */
	for(i=2;i<size;i++) {
		f = Conv(f,tmp);
		for(j=0;j<i;j++)
			gsl_matrix_set(C,i,j,gsl_vector_get(f,(f->size-1)/2+j));
	}

	/* Free memory */
	gsl_vector_free(tmp);
	gsl_vector_free(f);

	return C;
}






/**
 * Function MA
 *
 * Moving average smoothing (inplace)
 *
 * @param vector original vector
 * @param length smoothing length in samples
 */
void MA(gsl_vector *vector, int length) {

	int i,j;
	double sum;

	if(length == 0)
		return;
	if(length%2 == 0) {
		printf("Warning: Span of the moving average filter must be odd.\n");
		printf("         Span is changed to N-1.\n");
		length = length - 1;
	}
	if(length < 2) {
		printf("Warning: Span of the moving average filter must be at least 3.\n");
		return;
	}
	if(vector->size < length)
		return;

	/* Copy vector */
	gsl_vector *vector_orig = gsl_vector_alloc(vector->size);
	gsl_vector_memcpy(vector_orig,vector);

	/* Filter */
	for(i=length;i<vector->size;i++) {
    	sum = 0;
    	for(j=0;j<length;j++) {
    		sum += gsl_vector_get(vector_orig,i-j);
    	}
    	gsl_vector_set(vector,i-length/2,sum/length);
	}

	/* Fix the beginning */
	for(i=1;i<length+1;i=i+2) {
		sum = 0;
		for(j=0;j<i;j++) {
			sum += gsl_vector_get(vector_orig,j);
		}
		gsl_vector_set(vector,(i-1)/2,sum/i);
	}

	/* Fix the end */
	for(i=1;i<length+1;i=i+2) {
		sum = 0;
		for(j=0;j<i;j++) {
			sum += gsl_vector_get(vector_orig,vector->size-1-j);
		}
		gsl_vector_set(vector,vector->size-1-(i-1)/2,sum/i);
	}

	/* Free memory */
	gsl_vector_free(vector_orig);
}








/**
 * Function MA_voiced
 *
 * Moving average smoothing (inplace) for voiced regions (fundf > 0)
 *
 * @param vector original vector
 * @param fundf f0 vector
 * @param length smoothing length in samples for each sample
 */
void MA_voiced(gsl_vector *vector, gsl_vector *fundf, int len) {

	int i,j,k;
	gsl_vector *temp;
	for(i=0;i<vector->size;i++) {
		if(gsl_vector_get(fundf,i) > 0) {
			j = 0;
			while(i+j < vector->size && gsl_vector_get(fundf,i+j) > 0)
				j++;
			if(j>2) {
				temp = gsl_vector_alloc(j);
				for(k=0;k<j;k++)
					gsl_vector_set(temp,k,gsl_vector_get(vector,i+k));
				MA(temp,len);
				for(k=0;k<j;k++)
					gsl_vector_set(vector,i+k,gsl_vector_get(temp,k));
				gsl_vector_free(temp);
			}
			i = i + j;
		}
	}
}




/**
 * Function MA_unvoiced
 *
 * Moving average smoothing (inplace) for unvoiced regions (fundf = 0)
 *
 * @param vector original vector
 * @param fundf f0 vector
 * @param length smoothing length in samples for each sample
 */
void MA_unvoiced(gsl_vector *vector, gsl_vector *fundf, int len) {

	int i,j,k;
	gsl_vector *temp;
	for(i=0;i<vector->size;i++) {
		if(gsl_vector_get(fundf,i) == 0) {
			j = 0;
			while(i+j < vector->size && gsl_vector_get(fundf,i+j) == 0)
				j++;
			if(j>2) {
				temp = gsl_vector_alloc(j);
				for(k=0;k<j;k++)
					gsl_vector_set(temp,k,gsl_vector_get(vector,i+k));
				MA(temp,len);
				for(k=0;k<j;k++)
					gsl_vector_set(vector,i+k,gsl_vector_get(temp,k));
				gsl_vector_free(temp);
			}
			i = i + j;
		}
	}
}












/**
 * Function lsf2poly
 *
 * Convert LSF to polynomial
 *
 * @param lsf_matrix
 * @param poly
 * @param index
 * @param HMM switch for HMM
 */
void lsf2poly(gsl_matrix *lsf_matrix, gsl_vector *poly, int index, int HMM) {

	int i,l = lsf_matrix->size2;
	gsl_vector *lsf_vector = gsl_vector_calloc(l);
	gsl_vector *fi_p = NULL, *fi_q = NULL;

	/* Copy values to vector */
	for(i=0;i<l;i++)
		gsl_vector_set(lsf_vector,i,gsl_matrix_get(lsf_matrix,index,i));

	/* Check the validity of LSF and fix found errors (HMM parameters) */
	if(HMM == 1) {
		LSF_fix_vector(lsf_vector);
	}

	/* Create fi_p and fi_q */
	if(l%2 == 0) {
		fi_p = gsl_vector_alloc(l/2);
		fi_q = gsl_vector_alloc(l/2);
		for(i=0;i<l;i=i+2) {
			gsl_vector_set(fi_p,i/2,gsl_vector_get(lsf_vector,i));
		}
		for(i=1;i<l;i=i+2) {
			gsl_vector_set(fi_q,(i-1)/2,gsl_vector_get(lsf_vector,i));
		}
	} else {
		fi_p = gsl_vector_alloc((l+1)/2);
		for(i=0;i<l;i=i+2) {
			gsl_vector_set(fi_p,i/2,gsl_vector_get(lsf_vector,i));
		}
		if((l-1)/2 > 0) {
			fi_q = gsl_vector_alloc((l-1)/2);
			for(i=1;i<l-1;i=i+2) {
				gsl_vector_set(fi_q,(i-1)/2,gsl_vector_get(lsf_vector,i));
			}
		}
	}

	/* Construct vectors P and Q */
	gsl_vector *cp = gsl_vector_calloc(3);
	gsl_vector *cq = gsl_vector_calloc(3);
	gsl_vector_add_constant(cp,1);
	gsl_vector_add_constant(cq,1);
	gsl_vector *P = gsl_vector_alloc(1);
	gsl_vector *Q = gsl_vector_alloc(1);
	gsl_vector_set(P,0,1);
	gsl_vector_set(Q,0,1);
	for(i=0;i<fi_p->size;i++) {
		gsl_vector_set(cp,1,-2*cos(gsl_vector_get(fi_p,i)));
		P = Conv(P,cp);
	}
	if((l-1)/2 > 0) {
		for(i=0;i<fi_q->size;i++) {
			gsl_vector_set(cq,1,-2*cos(gsl_vector_get(fi_q,i)));
			Q = Conv(Q,cq);
		}
	}

	/* Add trivial zeros */
	if(l%2 == 0) {
		gsl_vector *conv = gsl_vector_calloc(2);
		gsl_vector_add_constant(conv,1);
		P = Conv(P,conv);
		gsl_vector_set(conv,0,-1);
    	Q = Conv(Q,conv);
    	gsl_vector_free(conv);
	} else {
		gsl_vector *conv = gsl_vector_calloc(3);
		gsl_vector_set(conv,0,-1);
		gsl_vector_set(conv,2,1);
		Q = Conv(Q,conv);
		gsl_vector_free(conv);
	}

	/* Construct polynomial */
	for(i=1;i<P->size;i++) {
		gsl_vector_set(poly,P->size-i-1,0.5*(gsl_vector_get(P,i)+gsl_vector_get(Q,i)));
	}

	/* Free memory */
	gsl_vector_free(lsf_vector);
	gsl_vector_free(fi_p);
	if((l-1)/2 > 0)
		gsl_vector_free(fi_q);
	gsl_vector_free(cp);
	gsl_vector_free(cq);
	gsl_vector_free(P);
	gsl_vector_free(Q);
}















/**
 * Function Conv
 *
 * Convolve two vectors
 *
 * @param conv1
 * @param conv2
 */
gsl_vector *Conv(gsl_vector *conv1, gsl_vector *conv2) {

	int i,j,n = conv2->size;
	double sum;
	gsl_vector *result = gsl_vector_alloc(conv1->size+conv2->size-1);
	gsl_vector *temp = gsl_vector_calloc(conv1->size+conv2->size-1);

	/* Set coefficients to temp */
	for(i=0;i<conv1->size;i++) {
		gsl_vector_set(temp,i,gsl_vector_get(conv1,i));
	}

	/* FIR-filter (Convolution) */
	for(i=0;i<temp->size;i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, n-1); j++) {
			sum += gsl_vector_get(temp, i-j)*gsl_vector_get(conv2,j);
		}
		gsl_vector_set(result, i, sum);
	}

	/* Free memory */
	gsl_vector_free(temp);
	gsl_vector_free(conv1);

	return result;
}










/**
 * Function Postfilter
 *
 * Postfilter LSFs
 *
 * @param LSF
 * @param params
 */
void Postfilter(gsl_matrix *LSF, PARAM *params) {

	if(params->postfilter_method == POSTFILTER_ID_LSF)
		LSF_Postfilter(LSF,params->postfilter_alpha);
	else if(params->postfilter_method == POSTFILTER_ID_LPC) {
		LPC_Postfilter(LSF, params->postfilter_alpha, params->frame_length);
		Smooth_matrix(LSF,3);
	}
}









/**
 * Function LSF_Postfilter
 *
 * Apply formant enhancement to LSFs
 *
 * @param lsf LSF-matrix
 * @param alpha postfilter coefficient alpha
 */
void LSF_Postfilter(gsl_matrix *lsf, double alpha) {

	if(alpha > 0) {
		int i,j;
		double d[lsf->size2-1];
		for(i=0;i<lsf->size1;i++) {
			for(j=0;j<lsf->size2-1;j++) {
				d[j] = alpha*(gsl_matrix_get(lsf,i,j+1) - gsl_matrix_get(lsf,i,j));
				if(j>0) {
					gsl_matrix_set(lsf,i,j, gsl_matrix_get(lsf,i,j-1) + d[j-1] + (pow(d[j-1],2)/(pow(d[j-1],2) + pow(d[j],2))) * ( gsl_matrix_get(lsf,i,j+1) - gsl_matrix_get(lsf,i,j-1) - d[j] - d[j-1] ) );
				}
			}
		}
	}
}





/**
 * Function LPC_Postfilter
 *
 * Enhance Formants by modifying the re-evaluated the LPC power spectrum,
 * and evaluating the LPC-coefficients again
 *
 * @param LSF pointer to the LSF matrix
 * @param gamma enhancement coefficient
 * @param frame_length
 */
void LPC_Postfilter(gsl_matrix *LSF, double gamma, int frame_length) {

	int i,j,fi,s,nf,p = LSF->size2,n = POWER_SPECTRUM_FRAME_LEN;
	double A[n];
	double B[n];
	double data[n];
	gsl_vector *a = gsl_vector_alloc(p+1);
	gsl_vector *a_temp = gsl_vector_calloc(p);
	gsl_vector *r = gsl_vector_alloc(p+1);
	gsl_vector *rr = gsl_vector_alloc(n);
	gsl_vector *S = gsl_vector_alloc(n);
	gsl_vector *b = gsl_vector_alloc(p);
	gsl_matrix *R = gsl_matrix_alloc (p, p);
	gsl_vector *formants = gsl_vector_calloc(100);
	gsl_permutation *perm = gsl_permutation_alloc(p);
	gsl_complex ca;
	gsl_complex cb;
	gsl_fft_real_wavetable *wreal = gsl_fft_real_wavetable_alloc(n);
	gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(n);
	gsl_fft_complex_wavetable *cwt = gsl_fft_complex_wavetable_alloc(n);
	gsl_fft_complex_workspace *cwork = gsl_fft_complex_workspace_alloc(n);
	double xa[2*n];
	double xb[2*n];

	/* Initialize B */
	B[0] = 1;
	for(i=1;i<n;i++)
		B[i] = 0;
	gsl_fft_real_transform(B,1,n,wreal,work);
	gsl_complex_packed_array complex_coefficients_b = xb;
	gsl_fft_halfcomplex_unpack(B,complex_coefficients_b,1,n);

	/* Loop for every index of the LSF matrix */
	for(fi=0;fi<LSF->size1;fi++) {

		/* Convert LSF to LPC */
		lsf2poly(LSF,a,fi,1);

		/* Evaluate power spectrum S, this assumes an all-pole model (B = 1) */
		for(i=0;i<p+1;i++)
			A[i] = gsl_vector_get(a,i);
		for(i=p+1;i<n;i++)
			A[i] = 0;
		gsl_fft_real_transform(A,1,n,wreal,work);
		gsl_complex_packed_array complex_coefficients_a = xa;
		gsl_fft_halfcomplex_unpack(A,complex_coefficients_a,1,n);
		for(i=0;i<n;i++) {
			GSL_SET_COMPLEX(&ca, REAL(complex_coefficients_a,i), IMAG(complex_coefficients_a,i));
			GSL_SET_COMPLEX(&cb, REAL(complex_coefficients_b,i), IMAG(complex_coefficients_b,i));
			ca = gsl_complex_div(cb, ca);
			gsl_vector_set(S,i,gsl_complex_abs2(ca));
		}

		/* Modification of the power spectrum S */
		GetFormants(S,formants,&nf);
		ModPowerSpectrum(S,formants,nf,gamma);

		/* Construct autocorrelation r */
		for(i=0;i<n;i++)
			data[i] = gsl_vector_get(S,i);
		gsl_fft_real_unpack(data, complex_coefficients_a, 1, n);

		gsl_fft_complex_inverse(complex_coefficients_a, 1, n, cwt, cwork);
		for(i=0;i<2*n;i = i + 2)
			gsl_vector_set(rr,i/2,xa[i]);
		for(i=0;i<p+1;i++)
			gsl_vector_set(r,i,gsl_vector_get(rr,i));

		/* Construct LPC */
		for(i=0; i<p;i++)
			for(j=0; j<p; j++)
				gsl_matrix_set(R, i, j, gsl_vector_get(r, abs(i-j)));
		for(i=1; i<p+1; i++)
			gsl_vector_set(b, i-1, gsl_vector_get(r, i));
		gsl_linalg_LU_decomp(R, perm, &s);
		gsl_linalg_LU_solve(R, perm, b, a_temp);
		for(i=1; i<a->size; i++)
			gsl_vector_set(a, i, (-1.0)*gsl_vector_get(a_temp, i-1));
		gsl_vector_set(a, 0, 1);
		for(i=0;i<a->size;i++)
			if(gsl_isnan(gsl_vector_get(a,i)))
				gsl_vector_set(a,i,0);

		/* Convert LPC back to LSF */
		Convert_to_LSF(LSF, a, fi);
	}

	/* Free memory */
	gsl_vector_free(formants);
	gsl_vector_free(a);
	gsl_vector_free(a_temp);
	gsl_vector_free(r);
	gsl_vector_free(rr);
	gsl_vector_free(S);
	gsl_vector_free(b);
	gsl_matrix_free(R);
	gsl_permutation_free(perm);
	gsl_fft_real_wavetable_free(wreal);
	gsl_fft_real_workspace_free(work);
	gsl_fft_complex_wavetable_free(cwt);
	gsl_fft_complex_workspace_free(cwork);
}






/**
 * Function ModPowerSpectrum
 *
 * Modify power spectrum in order to enhance formants
 *
 * @param s pointer to power spectrum vector
 * @param gamma enhancement coefficient
 */
void ModPowerSpectrum(gsl_vector *s, gsl_vector *formants, int n, double gamma) {

	int i,j;

	/* Nonlinearity in power reduction depending on the width of the valley */
	double l = 150.0;
	double d = 40.0;
	double add = 0.5 + gamma;
	double c = 0.5;
	int dist;
	double mod;

	/* Modify spectrum between zero and the first formant */
	dist = gsl_vector_get(formants,0);
	mod = c*(-1.0/(1.0 + exp((-dist+l)/d))) + add;
	for(i=0;i<gsl_vector_get(formants,0) - POWER_SPECTRUM_WIN;i++)
		gsl_vector_set(s,i,gsl_vector_get(s,i)*gamma);

	/* Modify spectrum between the last formant and FS/2 */
	dist = floor(s->size/2)-gsl_vector_get(formants,n-1);
	mod = c*(-1.0/(1.0 + exp((-dist+l)/d))) + add;
	for(i=gsl_vector_get(formants,n-1) + POWER_SPECTRUM_WIN+1;i<s->size/2+1;i++)
		gsl_vector_set(s,i,gsl_vector_get(s,i)*gamma);

	/* Modify spectrum within a constant number of bins from the formant peaks */
	for(i=0;i<n-1;i++) {
		dist = gsl_vector_get(formants,i+1) - gsl_vector_get(formants,i);
		mod = c*(-1.0/(1.0 + exp((-dist+l)/d))) + add;
		for(j=gsl_vector_get(formants,i) + POWER_SPECTRUM_WIN+1;j<gsl_vector_get(formants,i+1) - POWER_SPECTRUM_WIN;j++)
			gsl_vector_set(s,j,gsl_vector_get(s,j)*gamma);
	}

	/* Reconstruct image spectrum */
	if((s->size)%2 == 0)
		for(i=1;i<s->size/2;i++)
			gsl_vector_set(s,s->size-i,gsl_vector_get(s,i));
	else
		for(i=1;i<ceil(s->size/2)+1;i++)
			gsl_vector_set(s,s->size-i,gsl_vector_get(s,i));
}







/**
 * Function GetFormants
 *
 * Get formant positions from smooth power spectrum
 *
 * @param s pointer to power spectrum vector
 * @param formants pointer to forman position vector
 */
void GetFormants(gsl_vector *s, gsl_vector *formants, int *n) {

	int i,ind;
	gsl_vector *sd = gsl_vector_alloc(s->size);

	/* Differentiate s */
	gsl_vector_memcpy(sd,s);
	Differentiate_noleak(sd);

	/* Find formant peaks */
	ind = 0;
	for(i=0;i<round(s->size/2)+1;i++) {
		if(gsl_vector_get(sd,i) >= 0 && gsl_vector_get(sd,i+1) <= 0) {
			gsl_vector_set(formants,ind,i);
			ind++;
		}
	}

	/* Set the number of formants */
	(*n) = ind;

	/* Free memory */
	gsl_vector_free(sd);
}


























/**
 * Function Hp_filt_f0
 *
 * High-pass filter speech below F0
 *
 * @param signal speech
 * @param fundf F0
 *
 */
void Hp_filt_below_f0(gsl_vector *signal, gsl_vector *fundf, PARAM *params) {

	if(params->hpfiltf0 == 0)
		return;

	int i,j,fnumber,use_hmm = 0,filter_update_interval = 1;
	double f0 = 0,lambda = 0,weight;

	/* Filter cut-off frequencies */
	double f[25] = {0,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500};

	/* Denominator filter coefficients */
	double lsfa[24][7] =    {{0.012270374463510,   0.015772995424344,   0.024530703714270,   0.080413887423843,   0.400306289905254,   1.168999898025715,   2.132579893858550},
							 {0.017920282072328,   0.022791623836745,   0.029864786258385,   0.055369341235469,   0.250523957519273,   1.092123690431822,   2.105116976227308},
	   					     {0.021801057749634,   0.028684712190248,   0.034646162981531,   0.057546773586775,   0.246163207189489,   1.093903629649360,   2.108370657386230},
	   					     {0.030330579315104,   0.038234559055245,   0.050121987909179,   0.091732027364528,   0.329220914237973,   1.125584823436989,   2.116231821108477},
	   					     {0.034811643040824,   0.044349408406948,   0.053516371808008,   0.087458975943472,   0.308231588488585,   1.118896822551608,   2.116399287588925},
	   					     {0.043917420002783,   0.054747497330268,   0.074390982686771,   0.143405359634844,   0.436338014468054,   1.182631435560068,   2.136140897442202},
	   					     {0.045038897777332,   0.058258259000167,   0.068963787550062,   0.109221377061591,   0.344220168514159,   1.134851229904599,   2.124295109888979},
	   					     {0.051022535502281,   0.065831379814600,   0.077573712097764,   0.121610480053570,   0.364861933549680,   1.146550846159653,   2.126584406920744},
	   					     {0.057632798735370,   0.073608263773643,   0.086404993167967,   0.134181910858372,   0.384477060556410,   1.154156703313648,   2.134888076326625},
	   					     {0.064123563307627,   0.081389884334947,   0.095130140196939,   0.146329709620344,   0.403020870485622,   1.163653230220413,   2.139554601837831},
	   					     {0.070489123224552,   0.089082609042766,   0.103643756411925,   0.157832024397575,   0.420045784848692,   1.172254221962511,   2.139045944360825},
	   					     {0.076847236579249,   0.096765479627175,   0.112050759011413,   0.168910344402175,   0.436071212884462,   1.181903820775476,   2.140464848942151},
	   					     {0.083670629006476,   0.104607866948277,   0.120599031060558,   0.180042602015452,   0.451542294904885,   1.189290523546819,   2.144017570035968},
	   					     {0.089932593812542,   0.112338083896512,   0.128836675422754,   0.190418798017626,   0.465709447605648,   1.200408605940246,   2.144702419538647},
	   					     {0.097482721230326,   0.120403139593136,   0.137524292415345,   0.201381669315467,   0.480314653853837,   1.205081242188550,   2.149181575561907},
	   					     {0.104600479029988,   0.128368705707085,   0.145929375571240,   0.211617389775297,   0.493657999260862,   1.212296579325581,   2.151596677897734},
	   					     {0.112590768888632,   0.136542902549503,   0.154586373033906,   0.222085490310900,   0.506940958793273,   1.217511577523374,   2.163600455800026},
	   					     {0.120305856730245,   0.144659614898539,   0.163016448320984,   0.231948743065545,   0.519466412869994,   1.224899556795356,   2.172566273838793},
	   					     {0.126196196114728,   0.152429039122201,   0.170619462400073,   0.240153856910034,   0.529226373592494,   1.235527230693612,   2.155469675089261},
	   					     {0.135957950019333,   0.160953029397963,   0.179571203980533,   0.250426168574974,   0.541985366804235,   1.235994160018460,   2.179609273087876},
	   					     {0.144100745276714,   0.169162818316366,   0.187753461107828,   0.259137041396872,   0.552306542274397,   1.240620073961813,   2.183055012971304},
	   					     {0.152457690627618,   0.177406803994366,   0.195849185376123,   0.267429744429561,   0.561963473760917,   1.244674940290760,   2.186392580112432},
	   					     {0.161053099588424,   0.185684023980750,   0.203849642456801,   0.275267011545265,   0.570885488302771,   1.248003394526696,   2.189775517375000},
	   					     {0.169933394492268,   0.193993666339236,   0.211740771952925,   0.282593310976812,   0.579087806020812,   1.250410947768047,   2.193033533970162}};

	/* Numerator filter coefficients */
	double lsfb[24][5] =    {{0.000030707852178,   0.003290687909261,   0.008730452554971,   0.008758594958107,   0.638065301730947},
			   				 {0.000345919101634,   0.007464820245500,   0.013446868902239,   0.013524918943862,   0.252043490985464},
							 {0.009337032323646,   0.012733180521797,   0.018752843006707,   0.018887151908989,   0.124438874030317},
			   			     {0.001750143943771,   0.014695620878432,   0.023413418898011,   0.023932233792984,   0.348241402501227},
							 {0.007451258505131,   0.023588219512075,   0.028822776612533,   0.029638359563087,   0.209598248087015},
			   			     {0.001262844775076,   0.019228486757771,   0.034126483729823,   0.034313458775845,   0.533669297472435},
				   			 {0.003966950786406,   0.027567085337588,   0.030762627018637,   0.041065592365883,   0.041079101614055},
				   			 {0.007300878835979,   0.031740322902411,   0.031817781139977,   0.047113960631324,   0.047413638436212},
					   		 {0.013297800640636,   0.036151790887025,   0.037029955814879,   0.053378560084371,   0.053478127353236},
					   		 {0.004245482517036,   0.040719988302241,   0.041791183344858,   0.059844792535173,   0.059931017064342},
                             {0.013325335393929,   0.045481683243416,   0.046722772464647,   0.066520788321735,   0.066692797195485},
	                         {0.011043680418257,   0.050443635379848,   0.052305089861292,   0.073408015773066,   0.073456841095412},
                             {0.010460216064980,   0.055656179782247,   0.056019990537994,   0.080510475101328,   0.080529142221492},
	                         {0.005042256738375,   0.060986479133597,   0.061700727238748,   0.087811525503098,   0.087832573628210},
                             {0.021201154212601,   0.066716566424458,   0.070923784343375,   0.095350487083348,   0.095396267816331},
	                         {0.021743620358827,   0.072604498913997,   0.080230354225900,   0.103093601812379,   0.103095461338701},
                             {0.010361144083279,   0.078885501870204,   0.081497400090665,   0.111068522558155,   0.111953788509510},
	                         {0.035482770544817,   0.085365388042439,   0.094373216452207,   0.119258245225241,   0.119606909976157},
                             {0.008760915951485,   0.091788399421470,   0.093320007127655,   0.127626103353867,   0.127670078949185},
	                         {0.014865429114325,   0.099220356759925,   0.099842903153314,   0.136297513507208,   0.136872182683673},
                             {0.010181473814981,   0.106654721678083,   0.107407671710693,   0.145160461404689,   0.145544940575618},
	                         {0.013092604335205,   0.114476964121103,   0.116199737975443,   0.154255878946170,   0.154357148101622},
                             {0.019538269162167,   0.122717547053384,   0.127482889764472,   0.163581821500224,   0.164070219003255},
							 {0.006517570557916,   0.131457491434920,   0.135046228721468,   0.173152990621689,   0.173468137485191}};

	/* Allocate filter matrices */
	gsl_matrix *LSFA = gsl_matrix_alloc(fundf->size,7);
	gsl_matrix *LSFB = gsl_matrix_alloc(fundf->size,5);

	/* Find appropriate filter coefficients according to F0 */
	for(i=0;i<fundf->size;i++) {

		/* Fundamental frequency */
		if(gsl_vector_get(fundf,i) != 0)
			f0 = gsl_vector_get(fundf,i);
		else
			f0 = 40;

		/* Search for nearest upper cut-off frequency */
		fnumber = 0;
		while(f[fnumber] < f0 && fnumber != 24)
			fnumber++;

		/* Determine weight */
		if(fnumber == 0)
			weight = 1;
		else if(fnumber == 24)
			weight = 0;
		else
			weight = (f[fnumber]-f0)/(f[fnumber]-f[fnumber-1]);

		/* Interpolate between two filters */
		for(j=0;j<LSFA->size2;j++)
			gsl_matrix_set(LSFA, i, j, weight*lsfa[GSL_MAX(fnumber-2,0)][j] + (1.0-weight)*lsfa[GSL_MAX(fnumber-1,0)][j]);
		for(j=0;j<LSFB->size2;j++)
			gsl_matrix_set(LSFB, i, j, weight*lsfb[GSL_MAX(fnumber-2,0)][j] + (1.0-weight)*lsfb[GSL_MAX(fnumber-1,0)][j]);
	}

	/* Initialize filter */
	int q,mlen;
	long int o;
	double xr,x,ffr,tmpr,Bb;
	double *sigma;
	long int len = signal->size;
	int bdim = LSFA->size2 + 1;
	int adim = LSFB->size2 + 1;
	double *Ar = (double *)calloc(adim,sizeof(double));
	double *Br = (double *)calloc(bdim,sizeof(double));
	double *ynr = (double *)calloc(signal->size,sizeof(double));
	double *rsignal = (double *)calloc(signal->size,sizeof(double));
	double *rmem = (double *)calloc((bdim+2),sizeof(double));
	gsl_vector *A = gsl_vector_alloc(adim);
	gsl_vector *B = gsl_vector_alloc(bdim);
	gsl_matrix *LSFA_i = gsl_matrix_alloc(params->signal_length,LSFA->size2);
	gsl_matrix *LSFB_i = gsl_matrix_alloc(params->signal_length,LSFB->size2);

	/* Smooth and interpolate */
	Smooth_interp_lsf(LSFA_i,LSFA,len,use_hmm,5);
	Smooth_interp_lsf(LSFB_i,LSFB,len,use_hmm,5);

	/* Set signal to array */
	for(i=0;i<len;i++) {
		rsignal[i] = gsl_vector_get(signal,i);
	}

	/* Initialize */
	sigma = NFArray(bdim+2);
	Bb = 0;
	if(adim >= bdim)
		mlen = adim;
	else
		mlen = bdim + 1;

	/* Warped filtering */
	for(o=0;o<len;o++) {

		/* Update filter coefficients */
		if(o%filter_update_interval == 0) {
			lsf2poly(LSFA_i,B,o,params->use_hmm);
			lsf2poly(LSFB_i,A,o,params->use_hmm);
			for(i=0;i<adim;i++)
				Ar[i] = gsl_vector_get(A,i);
			for(i=0;i<bdim;i++)
				Br[i] = gsl_vector_get(B,i);
			alphas2sigmas(Br,sigma,lambda,bdim-1);
			Bb = 1/Br[0];
		}

		xr = rsignal[o]*Bb;

		/* Update feedbackward sum */
		for(q=0;q<bdim;q++) {
			xr -= sigma[q]*rmem[q];
		}
		xr = xr/sigma[bdim];
		x = xr*Ar[0];

		/* Update inner states */
		for(q=0;q<mlen;q++) {
			tmpr = rmem[q] + lambda*(rmem[q+1] - xr);
			rmem[q] = xr;
			xr = tmpr;
		}

		/* Update feedforward sum */
		for(q=0,ffr=0.0;q<adim-1;q++) {
			ffr += Ar[q+1]*rmem[q+1];
		}

	   /* Update output */
	   ynr[o] = x + ffr;
	}

	/* Set output to vector */
	for(i=0;i<len;i++) {
		gsl_vector_set(signal,i,ynr[i]);
	}

	/* Free memory */
	free(ynr);
	free(rsignal);
	free(Ar);
	free(Br);
	free(rmem);
	free(sigma);
	gsl_vector_free(A);
	gsl_vector_free(B);
	gsl_matrix_free(LSFA);
	gsl_matrix_free(LSFB);
	gsl_matrix_free(LSFA_i);
	gsl_matrix_free(LSFB_i);
}











/**
 * Function Evaluate_matrix_std
 *
 * Evaluate standard deviations of matrix
 *
 * @param std vector
 * @param data matrix
 */
void Evaluate_matrix_std(gsl_vector *std, gsl_matrix *data) {

	int i,j;
	gsl_vector *mean = gsl_vector_calloc(data->size2);

    /* Evaluate mean and std of matrix data */
    for(i=0;i<data->size2;i++)
    	for(j=0;j<data->size1;j++)
    		gsl_vector_set(mean,i,gsl_vector_get(mean,i) + gsl_matrix_get(data,j,i));
    for(i=0;i<data->size2;i++)
    	gsl_vector_set(mean,i,gsl_vector_get(mean,i)/data->size1);
    for(i=0;i<data->size2;i++)
		for(j=0;j<data->size1;j++)
			gsl_vector_set(std,i,gsl_vector_get(std,i) + powf(gsl_matrix_get(data,j,i)-gsl_vector_get(mean,i),2));
    for(i=0;i<data->size2;i++)
		gsl_vector_set(std,i,sqrt(gsl_vector_get(std,i)/(data->size1-1)));

    /* Free memory */
    gsl_vector_free(mean);
}



/**
 * Function Evaluate_vector_std
 *
 * Evaluate standard deviation of vector
 *
 * @param std
 * @param data
 */
void Evaluate_vector_std(gsl_vector *std, gsl_vector *data) {

	int i;
	double mean = 0;

    /* Evaluate mean and std of vector data */
    for(i=0;i<data->size;i++)
		mean += gsl_vector_get(data,i);
    mean = mean/data->size;
    for(i=0;i<data->size;i++)
    	gsl_vector_set(std,0,gsl_vector_get(std,0) + powf(gsl_vector_get(data,i)-mean,2));
	gsl_vector_set(std,0,sqrt(gsl_vector_get(std,0)/(data->size-1)));
}










/**
 * Function ReadFileDouble
 *
 * Read double values from file
 *
 * @param name filename
 * @return vector containing the values
 */
gsl_vector *ReadFileDouble(char *name) {

	FILE *file;
	int fileLen,n,i;

	/* Open file */
	file = fopen(name, "rb"); // Open as binary
	if(!file) {
		printf("Error opening file %s: %s\n", name, strerror(errno));
		return NULL;
	}

	/* Get file length */
	fseek(file, 0, SEEK_END);
	fileLen = ftell(file);
	fseek(file, 0, SEEK_SET);

	/* Allocate memory */
	double *buffer = (double *)malloc(fileLen);
	if(!buffer) {
		printf("Memory error!\n");
        fclose(file);
		return NULL;
	}

	/* Read file contents into buffer */
	n = fread(buffer, fileLen, 1, file);
	fclose(file);

	/* Set values to vector */
	gsl_vector *vector = gsl_vector_calloc(fileLen/sizeof(double));
	for(i=0;i<fileLen/sizeof(double);i++) {
		gsl_vector_set(vector,i,buffer[i]);
	}
	free(buffer);
	return vector;
}





/**
 * Function ReadFileFloat
 *
 * Read float values from file
 *
 * @param name filename
 * @return vector containing the values
 */
gsl_vector *ReadFileFloat(char *name) {

	FILE *file;
	int fileLen,n,i;

	/* Open file */
	file = fopen(name, "rb"); // Open as binary
	if(!file) {
		printf("Error opening file %s: %s\n", name, strerror(errno));
		return NULL;
	}

	/* Get file length */
	fseek(file, 0, SEEK_END);
	fileLen = ftell(file);
	fseek(file, 0, SEEK_SET);

	/* Allocate memory */
	float *buffer = (float *)malloc(fileLen);
	if(!buffer) {
		printf("Memory error!\n");
        fclose(file);
		return NULL;
	}

	/* Read file contents into buffer */
	n = fread(buffer, fileLen, 1, file);
	fclose(file);

	/* Set values to vector */
	gsl_vector *vector = gsl_vector_calloc(fileLen/sizeof(float));
	for(i=0;i<fileLen/sizeof(float);i++) {
		gsl_vector_set(vector,i,buffer[i]);
	}
	free(buffer);
	return vector;
}










/**
 * Function Gain_normalization
 *
 * Normalize the gain of a signal.
 * Result is saved in-place to vector signal
 *
 * @param signal input signal
 * @param gain gain vector
 * @param frame_length frame length in samples
 * @param shift shift length in samples
 * @param speed syntesis speed
 * @param pitch synthesis pitch
 * @param gain_threshold threshold for gain normalization
 *
 */
void Gain_normalization(gsl_vector *signal,
						gsl_vector *gain,
						int frame_length,
						int shift,
						double speed,
						double pitch,
						double gain_threshold) {

	int i,j;
	double sum;
	gsl_vector *norm = gsl_vector_calloc(gain->size);
	gsl_vector *norm_i = gsl_vector_alloc(signal->size);

	/* Calculate gain normalization vector */
	for(i=0;i<gain->size;i++) {
		sum = 0;
		for(j=i*shift/speed;j<i*shift/speed+frame_length/speed;j++) {
			sum += pow(gsl_vector_get(signal,GSL_MIN(signal->size-1,j))*HANN(j-i*shift/speed,frame_length),2);
		}
		if(sum < gain_threshold) {sum = gain_threshold;}
		gsl_vector_set(norm,i,sqrt((E_REF*powf(10.0,gsl_vector_get(gain,i)/10.0))/sum));
	}

	/* Linear interpolation */
	Interpolate_lin(norm,norm_i);

	/* Set gain */
	for(i=0;i<signal->size;i++) {
		gsl_vector_set(signal,i,gsl_vector_get(signal,i)*gsl_vector_get(norm_i,i));
	}

	/* Free memory */
	gsl_vector_free(norm);
	gsl_vector_free(norm_i);
}





/**
 * Function LSF_MOD
 *
 * Modify LSF matrix
 *
 * @param matrix matrix to be modified
 *
 */
void LSF_MOD(gsl_matrix *lsf) {

	int i,j;
	double d = 1.4;
	for(i=0;i<lsf->size1;i++) {
		for(j=0;j<lsf->size2;j++) {
			gsl_matrix_set(lsf,i,j, gsl_matrix_get(lsf,i,j)/d);
		}
	}
	LSF_fix_matrix(lsf);
}








/**
 * Function Integrate_matrix
 *
 * Integrate matrix
 *
 * @param matrix matrix to be integrated
 *
 */
void Integrate_matrix(gsl_matrix *matrix) {

	int i,j;
	for(i=0;i<matrix->size1;i++)
		for(j=1;j<matrix->size2;j++)
			gsl_matrix_set(matrix,i,j,gsl_matrix_get(matrix,i,j-1) + gsl_matrix_get(matrix,i,j));
}


/**
 * Function Integrate
 *
 * Integrate vector
 *
 * @param vector to be integrated
 *
 */
void Integrate(gsl_vector *vector, double leak) {

	int i;
	for(i=1;i<vector->size;i++)
		gsl_vector_set(vector,i,gsl_vector_get(vector,i-1)*leak + gsl_vector_get(vector,i));
}



/**
 * Function Free_pulselib_variables
 *
 * Free pulse library variables
 *
 * @param
 */
void Free_pulselib_variables(gsl_matrix *pulses, gsl_matrix *pulses_rs, gsl_matrix *pwaveform, gsl_vector *pulse_lengths,
		gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *pharm, gsl_matrix *phnr, gsl_vector *pgain, gsl_vector *ph1h2,
		gsl_vector *pnaq, gsl_vector *pca_mean, gsl_matrix *pca_pc, gsl_matrix *pca_w_lib, gsl_vector *stoch_env, gsl_vector *stoch_sp, PARAM *params) {

	/* Do not free if pulse library is not in use */
	if(params->use_pulselib == 0)
		return;

	/* Free pulse library */
	gsl_matrix_free(pulses);
	gsl_matrix_free(pulses_rs);
	gsl_vector_free(pulse_lengths);
	gsl_matrix_free(plsf);
	gsl_vector_free(pgain);
	if(params->use_waveform == 1) gsl_matrix_free(pwaveform);
	if(params->use_tilt == 1) gsl_matrix_free(ptilt);
	if(params->use_harmonics == 1) gsl_matrix_free(pharm);
	if(params->use_hnr == 1) gsl_matrix_free(phnr);
	if(params->use_h1h2 == 1) gsl_vector_free(ph1h2);
	if(params->use_naq == 1) gsl_vector_free(pnaq);

	/* Free pulse library pca parameters */
	if(params->use_pulselib_pca == 1) {
		gsl_vector_free(pca_mean);
		gsl_matrix_free(pca_pc);
	}

	/* Free pulse pca parameters */
	if(params->use_pulse_pca == 1)
		gsl_matrix_free(pca_w_lib);
}







/**
 * Function Free_variables
 *
 * Free synthesis variables
 *
 * @param
 */
void Free_variables(gsl_vector *original_pulse, gsl_vector *excitation_voiced, gsl_vector *excitation_unvoiced, gsl_vector *fundf, gsl_vector *gain,
		gsl_vector *gain_new, gsl_matrix *LSF, gsl_matrix *LSF2, gsl_matrix *LSF_interp, gsl_matrix *glflowsp, gsl_matrix *glflowsp_new,
		gsl_matrix *hnr, gsl_matrix *hnr_new, gsl_matrix *harmonics, gsl_matrix *waveform, gsl_vector *h1h2, gsl_vector *naq,
		gsl_vector *resynthesis_pulse_index, gsl_vector *pulse_clus_id, gsl_matrix *pulse_clusters, gsl_matrix *pca_w, 
		gsl_matrix **DNN_W, gsl_vector **input_minmax, gsl_vector *dnnpulseindices, gsl_vector *dnnpulses, PARAM *params) {

	/* Free variables */
	gsl_vector_free(original_pulse);
	gsl_vector_free(excitation_voiced);
	gsl_vector_free(excitation_unvoiced);
	gsl_vector_free(fundf);
	gsl_vector_free(gain);
	gsl_vector_free(gain_new);
	gsl_matrix_free(LSF);
	gsl_matrix_free(glflowsp_new);
	gsl_matrix_free(hnr_new);
	gsl_vector_free(resynthesis_pulse_index);
	if(LSF_interp != NULL) gsl_matrix_free(LSF_interp);
	if(params->use_tilt == 1 && glflowsp != NULL) gsl_matrix_free(glflowsp);
	if(params->use_hnr == 1 && hnr != NULL) gsl_matrix_free(hnr);
	if(params->use_harmonics == 1 && harmonics != NULL) gsl_matrix_free(harmonics);
	if(params->use_waveform == 1 && waveform != NULL) gsl_matrix_free(waveform);
	if(params->use_h1h2 == 1 && h1h2 != NULL) gsl_vector_free(h1h2);
	if(params->use_naq == 1 && naq != NULL) gsl_vector_free(naq);
	if(params->sep_vuv_spectrum == 1 && LSF2 != NULL) gsl_matrix_free(LSF2);
	if(params->pulse_clustering == 1) {
		if(pulse_clus_id != NULL) gsl_vector_free(pulse_clus_id);
		if(pulse_clusters != NULL)gsl_matrix_free(pulse_clusters);
	}

	/* Free pulse library PCA weights */
	if(params->use_pulselib_pca == 1)
		gsl_matrix_free(pca_w);

	/* Free DNN weights if used */
	int i;
	if(params->use_dnn_pulsegen == 1) {
		for(i=0;i<params->dnn_weight_dims->size/2;i++)
			gsl_matrix_free(DNN_W[i]);
		gsl_vector_free(dnnpulseindices);
		gsl_vector_free(dnnpulses);
		if(params->dnn_input_normalized == 1)
			gsl_vector_free(input_minmax[0]);
	}
	free(DNN_W);
	free(input_minmax);

	/* Free variables in param */
	for(i=0;i<params->synlistlen;i++)
		free(params->synlist[i]);
	free(params->synlist);
	free(params->synlist_filename);
	free(params->dnnpath);
	free(params->pulse_filename);
	free(params->pulselibrary_filename);
	gsl_vector_free(params->paramweights);
	gsl_vector_free(params->dnn_weight_dims);
}







/**
 * Function Save_excitation_to_wav
 *
 * Save excitation vector to wav file
 *
 * @param excitation vector
 * @param params parameters
 */
void Save_excitation_to_wav(gsl_vector *excitation, PARAM *params) {

	if(params->write_excitation_to_wav == 1) {
		gsl_vector *temp = gsl_vector_alloc(excitation->size);
		gsl_vector_memcpy(temp,excitation);
		Scale_signal(temp,SCALE_IF_GREATER_THAN_ONE);
		Save_signal_to_file(temp,params,FILENAME_ENDING_EXCITATION);
		gsl_vector_free(temp);
	}
}







/**
 * Function Normalize_pulse_library_var
 *
 * Normalize pulse library parameters according to synthesis parameters (normalize mean and variance)
 *
 * @param pulse lib params
 * @param synthesis params
 */
void Normalize_pulse_library_var(gsl_matrix *plsf,gsl_matrix *ptilt,gsl_matrix *pharm,gsl_matrix *phnr,gsl_matrix *pwaveform,
		gsl_vector *pgain,gsl_vector *ph1h2,gsl_vector *pnaq,gsl_matrix *lsf,gsl_matrix *tilt,gsl_matrix *harm,gsl_matrix *hnr,
		gsl_matrix *waveform,gsl_vector *gain,gsl_vector *h1h2,gsl_vector *naq, gsl_vector *fundf,PARAM *params) {

	/* Exit if pulse library is not used */
	if(params->use_pulselib == 0)
		return;

	/* Exit if normalization is not set on */
	if(params->normalize_pulselib == 0)
		return;

	/* Exit if adaptation to pulse library parameters is set on */
	if(params->adapt_to_pulselib == 1) {
		printf("Warning: Pulse library normalization and adaptation to pulse library parameters cannot be used simultaneously!\n");
		return;
	}

	/* Initialize */
	int i,j,k;
	double mean_lib,mean_par,std_lib,std_par;

	/* Normalize LSF */
	for(i=0;i<lsf->size2;i++) {

		/* Evaluate means */
		mean_par = 0;
		mean_lib = 0;
		for(j=0;j<plsf->size1;j++)
			mean_lib += gsl_matrix_get(plsf,j,i);
		k = 0;
		for(j=0;j<lsf->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				mean_par += gsl_matrix_get(lsf,j,i);
				k++;
			}
		}
		mean_lib = mean_lib/plsf->size1;
		mean_par = mean_par/k;

		/* Evaluate standard deviations */
		std_par = 0;
		std_lib = 0;
		for(j=0;j<plsf->size1;j++)
			std_lib += (gsl_matrix_get(plsf,j,i)-mean_lib)*(gsl_matrix_get(plsf,j,i)-mean_lib);
		std_lib = sqrt(std_lib/(plsf->size1-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<lsf->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_matrix_get(lsf,j,i)-mean_par)*(gsl_matrix_get(lsf,j,i)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;

		/* Normalize */
		for(j=0;j<plsf->size1;j++)
			gsl_matrix_set(plsf,j,i,(gsl_matrix_get(plsf,j,i)-mean_lib)/std_lib*std_par + mean_par);
	}

	/* Normalize LSFsource */
	for(i=0;i<tilt->size2;i++) {

		/* Evaluate means */
		mean_par = 0;
		mean_lib = 0;
		for(j=0;j<ptilt->size1;j++)
			mean_lib += gsl_matrix_get(ptilt,j,i);
		k = 0;
		for(j=0;j<tilt->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				mean_par += gsl_matrix_get(tilt,j,i);
				k++;
			}
		}
		mean_lib = mean_lib/ptilt->size1;
		mean_par = mean_par/k;

		/* Evaluate standard deviations */
		std_par = 0;
		std_lib = 0;
		for(j=0;j<ptilt->size1;j++)
			std_lib += (gsl_matrix_get(ptilt,j,i)-mean_lib)*(gsl_matrix_get(ptilt,j,i)-mean_lib);
		std_lib = sqrt(std_lib/(ptilt->size1-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<tilt->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_matrix_get(tilt,j,i)-mean_par)*(gsl_matrix_get(tilt,j,i)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;

		/* Normalize */
		for(j=0;j<ptilt->size1;j++)
			gsl_matrix_set(ptilt,j,i,(gsl_matrix_get(ptilt,j,i)-mean_lib)/std_lib*std_par + mean_par);
	}

	/* Normalize HNR */
	for(i=0;i<hnr->size2;i++) {

		/* Evaluate means */
		mean_par = 0;
		mean_lib = 0;
		for(j=0;j<phnr->size1;j++)
			mean_lib += gsl_matrix_get(phnr,j,i);
		k = 0;
		for(j=0;j<hnr->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				mean_par += gsl_matrix_get(hnr,j,i);
				k++;
			}
		}
		mean_lib = mean_lib/phnr->size1;
		mean_par = mean_par/k;

		/* Evaluate standard deviations */
		std_par = 0;
		std_lib = 0;
		for(j=0;j<phnr->size1;j++)
			std_lib += (gsl_matrix_get(phnr,j,i)-mean_lib)*(gsl_matrix_get(phnr,j,i)-mean_lib);
		std_lib = sqrt(std_lib/(phnr->size1-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<hnr->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_matrix_get(hnr,j,i)-mean_par)*(gsl_matrix_get(hnr,j,i)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;

		/* Normalize */
		for(j=0;j<phnr->size1;j++)
			gsl_matrix_set(phnr,j,i,(gsl_matrix_get(phnr,j,i)-mean_lib)/std_lib*std_par + mean_par);
	}

	/* Normalize Harmonics */
	if(params->use_harmonics == 1) {
		for(i=0;i<hnr->size2;i++) {

			/* Evaluate means */
			mean_par = 0;
			mean_lib = 0;
			for(j=0;j<pharm->size1;j++)
				mean_lib += gsl_matrix_get(pharm,j,i);
			k = 0;
			for(j=0;j<harm->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					mean_par += gsl_matrix_get(harm,j,i);
					k++;
				}
			}
			mean_lib = mean_lib/pharm->size1;
			mean_par = mean_par/k;

			/* Evaluate standard deviations */
			std_par = 0;
			std_lib = 0;
			for(j=0;j<pharm->size1;j++)
				std_lib += (gsl_matrix_get(pharm,j,i)-mean_lib)*(gsl_matrix_get(pharm,j,i)-mean_lib);
			std_lib = sqrt(std_lib/(pharm->size1-1));
			if(std_lib == 0)
				std_lib = 1;
			k = 0;
			for(j=0;j<harm->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					std_par += (gsl_matrix_get(harm,j,i)-mean_par)*(gsl_matrix_get(harm,j,i)-mean_par);
					k++;
				}
			}
			std_par = sqrt(std_par/(k-1));
			if(std_par == 0)
				std_par = 1;

			/* Normalize */
			for(j=0;j<pharm->size1;j++)
				gsl_matrix_set(pharm,j,i,(gsl_matrix_get(pharm,j,i)-mean_lib)/std_lib*std_par + mean_par);
		}
	}

	/* Normalize Waveform */
	if(params->use_waveform == 1) {
		for(i=0;i<hnr->size2;i++) {

			/* Evaluate means */
			mean_par = 0;
			mean_lib = 0;
			for(j=0;j<pwaveform->size1;j++)
				mean_lib += gsl_matrix_get(pwaveform,j,i);
			k = 0;
			for(j=0;j<waveform->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					mean_par += gsl_matrix_get(waveform,j,i);
					k++;
				}
			}
			mean_lib = mean_lib/pwaveform->size1;
			mean_par = mean_par/k;

			/* Evaluate standard deviations */
			std_par = 0;
			std_lib = 0;
			for(j=0;j<pwaveform->size1;j++)
				std_lib += (gsl_matrix_get(pwaveform,j,i)-mean_lib)*(gsl_matrix_get(pwaveform,j,i)-mean_lib);
			std_lib = sqrt(std_lib/(pwaveform->size1-1));
			if(std_lib == 0)
				std_lib = 1;
			k = 0;
			for(j=0;j<waveform->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					std_par += (gsl_matrix_get(waveform,j,i)-mean_par)*(gsl_matrix_get(waveform,j,i)-mean_par);
					k++;
				}
			}
			std_par = sqrt(std_par/(k-1));
			if(std_par == 0)
				std_par = 1;

			/* Normalize */
			for(j=0;j<pwaveform->size1;j++)
				gsl_matrix_set(pwaveform,j,i,(gsl_matrix_get(pwaveform,j,i)-mean_lib)/std_lib*std_par + mean_par);
		}
	}

	/* Normalize Gain */
	mean_lib = Mean(pgain);
	k = 0;
	mean_par = 0;
	for(i=0;i<gain->size;i++) {
		if(gsl_vector_get(fundf,i) > 0) {
			mean_par += gsl_vector_get(gain,i);
			k++;
		}
	}
	mean_par = mean_par/k;
	std_par = 0;
	std_lib = 0;
	for(j=0;j<pgain->size;j++)
		std_lib += (gsl_vector_get(pgain,j)-mean_lib)*(gsl_vector_get(pgain,j)-mean_lib);
	std_lib = sqrt(std_lib/(pgain->size-1));
	if(std_lib == 0)
		std_lib = 1;
	k = 0;
	for(j=0;j<gain->size;j++) {
		if(gsl_vector_get(fundf,j) > 0) {
			std_par += (gsl_vector_get(gain,j)-mean_par)*(gsl_vector_get(gain,j)-mean_par);
			k++;
		}
	}
	std_par = sqrt(std_par/(k-1));
	if(std_par == 0)
		std_par = 1;
	gsl_vector_add_constant(pgain,-mean_lib);
	gsl_vector_scale(pgain, std_par/std_lib);
	gsl_vector_add_constant(pgain,mean_par);

	/* Normalize H1H2 */
	if(params->use_h1h2 == 1) {
		mean_lib = Mean(ph1h2);
		k = 0;
		mean_par = 0;
		for(i=0;i<h1h2->size;i++) {
			if(gsl_vector_get(fundf,i) > 0) {
				mean_par += gsl_vector_get(h1h2,i);
				k++;
			}
		}
		mean_par = mean_par/k;
		std_par = 0;
		std_lib = 0;
		for(j=0;j<ph1h2->size;j++)
			std_lib += (gsl_vector_get(ph1h2,j)-mean_lib)*(gsl_vector_get(ph1h2,j)-mean_lib);
		std_lib = sqrt(std_lib/(ph1h2->size-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<h1h2->size;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_vector_get(h1h2,j)-mean_par)*(gsl_vector_get(h1h2,j)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;
		gsl_vector_add_constant(ph1h2,-mean_lib);
		gsl_vector_scale(ph1h2, std_par/std_lib);
		gsl_vector_add_constant(ph1h2,mean_par);
	}

	/* Normalize NAQ */
	if(params->use_naq == 1) {
		mean_lib = Mean(pnaq);
		k = 0;
		mean_par = 0;
		for(i=0;i<naq->size;i++) {
			if(gsl_vector_get(fundf,i) > 0) {
				mean_par += gsl_vector_get(naq,i);
				k++;
			}
		}
		mean_par = mean_par/k;
		std_par = 0;
		std_lib = 0;
		for(j=0;j<pnaq->size;j++)
			std_lib += (gsl_vector_get(pnaq,j)-mean_lib)*(gsl_vector_get(pnaq,j)-mean_lib);
		std_lib = sqrt(std_lib/(pnaq->size-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<naq->size;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_vector_get(naq,j)-mean_par)*(gsl_vector_get(naq,j)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;
		gsl_vector_add_constant(pnaq,-mean_lib);
		gsl_vector_scale(pnaq, std_par/std_lib);
		gsl_vector_add_constant(pnaq,mean_par);
	}
}






















/**
 * Function Adapt_synthesis_parameters_var
 *
 * Adapt synthesis parameters according to pulse library parameters (normalize mean and var)
 *
 * @param pulse lib params
 * @param synthesis params
 */
void Adapt_synthesis_parameters_var(gsl_matrix *plsf,gsl_matrix *ptilt,gsl_matrix *pharm,gsl_matrix *phnr,gsl_matrix *pwaveform,
		gsl_vector *pgain,gsl_vector *ph1h2,gsl_vector *pnaq,gsl_matrix *lsf,gsl_matrix *tilt,gsl_matrix *harm,gsl_matrix *hnr,
		gsl_matrix *waveform,gsl_vector *gain,gsl_vector *h1h2,gsl_vector *naq,gsl_vector *fundf,gsl_vector *pulse_lengths,PARAM *params) {


	/* Exit if pulse library is not used */
	if(params->use_pulselib == 0)
		return;

	/* Exit if adaptation is not set on */
	if(params->adapt_to_pulselib == 0)
		return;

	/* Exit if pulse library normalization is set on */
	if(params->normalize_pulselib == 1) {
		printf("Warning: Pulse library normalization and adaptation to pulse library parameters cannot be used simultaneously!\n");
		return;
	}

	/* Initialize */
	int i,j,k;
	double mean_lib,mean_par,std_lib,std_par;

	/* Normalize LSF */
	for(i=0;i<lsf->size2;i++) {

		/* Evaluate means */
		mean_par = 0;
		mean_lib = 0;
		for(j=0;j<plsf->size1;j++)
			mean_lib += gsl_matrix_get(plsf,j,i);
		k = 0;
		for(j=0;j<lsf->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				mean_par += gsl_matrix_get(lsf,j,i);
				k++;
			}
		}
		mean_lib = mean_lib/plsf->size1;
		mean_par = mean_par/k;

		/* Evaluate standard deviations */
		std_par = 0;
		std_lib = 0;
		for(j=0;j<plsf->size1;j++)
			std_lib += (gsl_matrix_get(plsf,j,i)-mean_lib)*(gsl_matrix_get(plsf,j,i)-mean_lib);
		std_lib = sqrt(std_lib/(plsf->size1-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<lsf->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_matrix_get(lsf,j,i)-mean_par)*(gsl_matrix_get(lsf,j,i)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;

		/* Normalize */
		for(j=0;j<lsf->size1;j++)
			if(gsl_vector_get(fundf,j) > 0)
				gsl_matrix_set(lsf,j,i,(gsl_matrix_get(lsf,j,i)-mean_par)/std_par*(std_lib*params->adapt_coeff + std_par*(1.0-params->adapt_coeff)) + params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
	}
	LSF_fix_matrix(lsf);

	/* Normalize LSFsource */
	for(i=0;i<tilt->size2;i++) {

		/* Evaluate means */
		mean_par = 0;
		mean_lib = 0;
		for(j=0;j<ptilt->size1;j++)
			mean_lib += gsl_matrix_get(ptilt,j,i);
		k = 0;
		for(j=0;j<tilt->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				mean_par += gsl_matrix_get(tilt,j,i);
				k++;
			}
		}
		mean_lib = mean_lib/ptilt->size1;
		mean_par = mean_par/k;

		/* Evaluate standard deviations */
		std_par = 0;
		std_lib = 0;
		for(j=0;j<ptilt->size1;j++)
			std_lib += (gsl_matrix_get(ptilt,j,i)-mean_lib)*(gsl_matrix_get(ptilt,j,i)-mean_lib);
		std_lib = sqrt(std_lib/(ptilt->size1-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<tilt->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_matrix_get(tilt,j,i)-mean_par)*(gsl_matrix_get(tilt,j,i)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;

		/* Normalize */
		for(j=0;j<tilt->size1;j++)
			if(gsl_vector_get(fundf,j) > 0)
				gsl_matrix_set(tilt,j,i,(gsl_matrix_get(tilt,j,i)-mean_par)/std_par*(std_lib*params->adapt_coeff + std_par*(1.0-params->adapt_coeff)) + params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
	}
	LSF_fix_matrix(tilt);

	/* Normalize HNR */
	for(i=0;i<hnr->size2;i++) {

		/* Evaluate means */
		mean_par = 0;
		mean_lib = 0;
		for(j=0;j<phnr->size1;j++)
			mean_lib += gsl_matrix_get(phnr,j,i);
		k = 0;
		for(j=0;j<hnr->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				mean_par += gsl_matrix_get(hnr,j,i);
				k++;
			}
		}
		mean_lib = mean_lib/phnr->size1;
		mean_par = mean_par/k;

		/* Evaluate standard deviations */
		std_par = 0;
		std_lib = 0;
		for(j=0;j<phnr->size1;j++)
			std_lib += (gsl_matrix_get(phnr,j,i)-mean_lib)*(gsl_matrix_get(phnr,j,i)-mean_lib);
		std_lib = sqrt(std_lib/(phnr->size1-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<hnr->size1;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_matrix_get(hnr,j,i)-mean_par)*(gsl_matrix_get(hnr,j,i)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;

		/* Normalize */
		for(j=0;j<hnr->size1;j++)
			gsl_matrix_set(hnr,j,i,(gsl_matrix_get(hnr,j,i)-mean_par)/std_par*(std_lib*params->adapt_coeff + std_par*(1.0-params->adapt_coeff)) + params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
	}

	/* Normalize Harmonics */
	if(params->use_harmonics == 1) {
		for(i=0;i<hnr->size2;i++) {

			/* Evaluate means */
			mean_par = 0;
			mean_lib = 0;
			for(j=0;j<pharm->size1;j++)
				mean_lib += gsl_matrix_get(pharm,j,i);
			k = 0;
			for(j=0;j<harm->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					mean_par += gsl_matrix_get(harm,j,i);
					k++;
				}
			}
			mean_lib = mean_lib/pharm->size1;
			mean_par = mean_par/k;

			/* Evaluate standard deviations */
			std_par = 0;
			std_lib = 0;
			for(j=0;j<pharm->size1;j++)
				std_lib += (gsl_matrix_get(pharm,j,i)-mean_lib)*(gsl_matrix_get(pharm,j,i)-mean_lib);
			std_lib = sqrt(std_lib/(pharm->size1-1));
			if(std_lib == 0)
				std_lib = 1;
			k = 0;
			for(j=0;j<harm->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					std_par += (gsl_matrix_get(harm,j,i)-mean_par)*(gsl_matrix_get(harm,j,i)-mean_par);
					k++;
				}
			}
			std_par = sqrt(std_par/(k-1));
			if(std_par == 0)
				std_par = 1;

			/* Normalize */
			for(j=0;j<harm->size1;j++)
				gsl_matrix_set(harm,j,i,(gsl_matrix_get(harm,j,i)-mean_par)/std_par*(std_lib*params->adapt_coeff + std_par*(1.0-params->adapt_coeff)) + params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
		}
	}

	/* Normalize Waveform */
	if(params->use_waveform == 1) {
		for(i=0;i<hnr->size2;i++) {

			/* Evaluate means */
			mean_par = 0;
			mean_lib = 0;
			for(j=0;j<pwaveform->size1;j++)
				mean_lib += gsl_matrix_get(pwaveform,j,i);
			k = 0;
			for(j=0;j<waveform->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					mean_par += gsl_matrix_get(waveform,j,i);
					k++;
				}
			}
			mean_lib = mean_lib/pwaveform->size1;
			mean_par = mean_par/k;

			/* Evaluate standard deviations */
			std_par = 0;
			std_lib = 0;
			for(j=0;j<pwaveform->size1;j++)
				std_lib += (gsl_matrix_get(pwaveform,j,i)-mean_lib)*(gsl_matrix_get(pwaveform,j,i)-mean_lib);
			std_lib = sqrt(std_lib/(pwaveform->size1-1));
			if(std_lib == 0)
				std_lib = 1;
			k = 0;
			for(j=0;j<waveform->size1;j++) {
				if(gsl_vector_get(fundf,j) > 0) {
					std_par += (gsl_matrix_get(waveform,j,i)-mean_par)*(gsl_matrix_get(waveform,j,i)-mean_par);
					k++;
				}
			}
			std_par = sqrt(std_par/(k-1));
			if(std_par == 0)
				std_par = 1;

			/* Normalize */
			for(j=0;j<waveform->size1;j++)
				gsl_matrix_set(waveform,j,i,(gsl_matrix_get(waveform,j,i)-mean_par)/std_par*(std_lib*params->adapt_coeff + std_par*(1.0-params->adapt_coeff)) + params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
		}
	}

	/* Normalize Gain */
	mean_lib = Mean(pgain);
	k = 0;
	mean_par = 0;
	for(i=0;i<gain->size;i++) {
		if(gsl_vector_get(fundf,i) > 0) {
			mean_par += gsl_vector_get(gain,i);
			k++;
		}
	}
	mean_par = mean_par/k;
	std_par = 0;
	std_lib = 0;
	for(j=0;j<pgain->size;j++)
		std_lib += (gsl_vector_get(pgain,j)-mean_lib)*(gsl_vector_get(pgain,j)-mean_lib);
	std_lib = sqrt(std_lib/(pgain->size-1));
	if(std_lib == 0)
		std_lib = 1;
	k = 0;
	for(j=0;j<gain->size;j++) {
		if(gsl_vector_get(fundf,j) > 0) {
			std_par += (gsl_vector_get(gain,j)-mean_par)*(gsl_vector_get(gain,j)-mean_par);
			k++;
		}
	}
	std_par = sqrt(std_par/(k-1));
	if(std_par == 0)
		std_par = 1;
	for(i=0;i<gain->size;i++)
		if(gsl_vector_get(fundf,i) > 0) // only voiced
			gsl_vector_set(gain,i,(gsl_vector_get(gain,i)-mean_par)/std_par*(std_lib*params->adapt_coeff + std_par*(1.0-params->adapt_coeff)) + params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
	//gsl_vector_add_constant(gain,-params->adapt_coeff*mean_par);
	//gsl_vector_scale(gain, (params->adapt_coeff*std_lib + (1.0-params->adapt_coeff)*std_par)/std_par);
	//gsl_vector_add_constant(gain,params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);

	/* Normalize H1H2 */
	if(params->use_h1h2 == 1) {
		mean_lib = Mean(ph1h2);
		k = 0;
		mean_par = 0;
		for(i=0;i<h1h2->size;i++) {
			if(gsl_vector_get(fundf,i) > 0) {
				mean_par += gsl_vector_get(h1h2,i);
				k++;
			}
		}
		mean_par = mean_par/k;
		std_par = 0;
		std_lib = 0;
		for(j=0;j<ph1h2->size;j++)
			std_lib += (gsl_vector_get(ph1h2,j)-mean_lib)*(gsl_vector_get(ph1h2,j)-mean_lib);
		std_lib = sqrt(std_lib/(ph1h2->size-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<h1h2->size;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_vector_get(h1h2,j)-mean_par)*(gsl_vector_get(h1h2,j)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;
		gsl_vector_add_constant(h1h2,-params->adapt_coeff*mean_par);
		gsl_vector_scale(h1h2, (params->adapt_coeff*std_lib + (1.0-params->adapt_coeff)*std_par)/std_par);
		gsl_vector_add_constant(h1h2,params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
	}

	/* Normalize NAQ */
	if(params->use_naq == 1) {
		mean_lib = Mean(pnaq);
		k = 0;
		mean_par = 0;
		for(i=0;i<naq->size;i++) {
			if(gsl_vector_get(fundf,i) > 0) {
				mean_par += gsl_vector_get(naq,i);
				k++;
			}
		}
		mean_par = mean_par/k;
		std_par = 0;
		std_lib = 0;
		for(j=0;j<pnaq->size;j++)
			std_lib += (gsl_vector_get(pnaq,j)-mean_lib)*(gsl_vector_get(pnaq,j)-mean_lib);
		std_lib = sqrt(std_lib/(pnaq->size-1));
		if(std_lib == 0)
			std_lib = 1;
		k = 0;
		for(j=0;j<naq->size;j++) {
			if(gsl_vector_get(fundf,j) > 0) {
				std_par += (gsl_vector_get(naq,j)-mean_par)*(gsl_vector_get(naq,j)-mean_par);
				k++;
			}
		}
		std_par = sqrt(std_par/(k-1));
		if(std_par == 0)
			std_par = 1;
		gsl_vector_add_constant(naq,-params->adapt_coeff*mean_par);
		gsl_vector_scale(naq, (params->adapt_coeff*std_lib + (1.0-params->adapt_coeff)*std_par)/std_par);
		gsl_vector_add_constant(naq,params->adapt_coeff*mean_lib + (1.0-params->adapt_coeff)*mean_par);
	}

	/* Normalize f0 */
	mean_lib = 0;
	for(i=0;i<pulse_lengths->size;i++)
		mean_lib += log10(2.0*params->FS/gsl_vector_get(pulse_lengths,i));
	mean_lib = mean_lib/pulse_lengths->size;
	k = 0;
	mean_par = 0;
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) > 0) {
			mean_par += log10(gsl_vector_get(fundf,i));
			k++;
		}
	}
	mean_par = mean_par/k;
	for(i=0;i<fundf->size;i++)
		gsl_vector_set(fundf,i,pow(10,log10(gsl_vector_get(fundf,i)) + params->adapt_coeff*(-mean_par+mean_lib)));
}





/**
 * Select_LSFs_from_pulse_library
 *
 * Replace vocal tract LSFs with the closes LSFs from the pulse library
 *
 * @param lsf original LSFs
 * @param plsf pulse library LSFs
 * @param fundf fundamental frequency
 * @param params
 */
void Select_LSFs_from_pulse_library(gsl_matrix *lsf, gsl_matrix *plsf, gsl_vector *fundf, int win, PARAM *params) {

	if(params->use_pulselib == 0)
		return;
	if(params->use_pulselib_lsf == 0)
		return;

	int i,j,k,n,minind;
	double mean_par, mean_lib;
	gsl_vector *error = gsl_vector_alloc(plsf->size1);
	gsl_matrix *lsf_tmp = gsl_matrix_alloc(lsf->size1,lsf->size2);

	/* Replace vocal tract LSFs with the closest LSFs from the pulse library */
	for(i=0;i<lsf->size1;i++) {
		if(gsl_vector_get(fundf,i) > 0) {
			gsl_vector_set_zero(error);
			for(n=0;n<plsf->size1;n++)
				for(j=0;j<lsf->size2;j++)
					gsl_vector_set(error,n,gsl_vector_get(error,n) + powf(gsl_matrix_get(lsf,i,j)-gsl_matrix_get(plsf,n,j),2));
			minind = gsl_vector_min_index(error);
			for(j=0;j<lsf->size2;j++)
				gsl_matrix_set(lsf_tmp,i,j,gsl_matrix_get(plsf,minind,j));
		}
	}

	/* Average between original LSFs and new LSFs from pulse library */
	for(i=0;i<lsf->size2;i++) {
		for(j=0;j<lsf->size1;j++) {
			if(gsl_vector_get(fundf,j) == 0)
				continue;
			n = 0;
			mean_par = 0;
			mean_lib = 0;
			for(k=j-win;k<j+win;k++) {
				if(k > 0 && k < lsf->size1 && gsl_vector_get(fundf,k) > 0) {
					mean_lib += gsl_matrix_get(lsf_tmp,k,i);
					mean_par += gsl_matrix_get(lsf,k,i);
					n++;
				}
			}
			if(n > 0) {
				mean_lib = mean_lib/n;
				mean_par = mean_par/n;
				gsl_matrix_set(lsf,j,i,gsl_matrix_get(lsf,j,i)-mean_par + mean_lib);
			}
		}
	}

	/* Free memory */
	gsl_vector_free(error);
	gsl_matrix_free(lsf_tmp);
}







/**
 * Function Convert_logF0_to_lin
 *
 * Convert fundamental frequency vector from natural logarithm to linear scale
 *
 * @param fundf
 * @param params
 */
void Convert_logF0_to_lin(gsl_vector *fundf, PARAM *params) {

	if(params->logf0 == 0)
		return;

	int i;
	for(i=0;i<fundf->size;i++)
		gsl_vector_set(fundf,i,exp(gsl_vector_get(fundf,i)));
}






























































































/*****************************************************************************/
/*                          FUNCTIONS NOT IN USE                             */
/*****************************************************************************/


/**
 * Function Wrap
 *
 * Wrap phase angles between [-pi,pi].
 *
 * @param phase vector of phase values
 *
 */
void Wrap(gsl_vector *phase) {

	int i;
	for(i=0;i<phase->size;i++) {
		while(gsl_vector_get(phase,i) < -M_PI)
			gsl_vector_set(phase,i,gsl_vector_get(phase,i) + 2*M_PI);
		while(gsl_vector_get(phase,i) > M_PI)
			gsl_vector_set(phase,i,gsl_vector_get(phase,i) - 2*M_PI);
	}
}




/**
 * Function Unwrap
 *
 * Unwrap phase angles. Algorithm minimizes the incremental phase variation
 * by constraining it to the range [-pi,pi]
 *
 * @param phase vector of phase values
 *
 */
void Unwrap(gsl_vector *phase) {

	int i;
	gsl_vector *dphase = gsl_vector_calloc(phase->size);
	gsl_vector *dphases = gsl_vector_calloc(phase->size);

	/* Evaluate incremental phase variation */
	for(i=0;i<phase->size-1;i++)
		gsl_vector_set(dphase,i,gsl_vector_get(phase,i+1)-gsl_vector_get(phase,i));

	/* Evaluate equivalent phase variations in [-pi,pi) */
	for(i=0;i<phase->size;i++)
		gsl_vector_set(dphases,i,((gsl_vector_get(dphase,i)+M_PI)/(2*M_PI)-floor((gsl_vector_get(dphase,i)+M_PI)/(2*M_PI)))*2*M_PI - M_PI);

	/* Preserve variation sign for pi vs. -pi */
	for(i=0;i<phase->size;i++) {
		if(gsl_vector_get(dphases,i) == -M_PI && gsl_vector_get(dphase,i) > 0)
			gsl_vector_set(dphases,i,M_PI);
	}

	/* Incremental phase corrections */
	for(i=0;i<phase->size;i++)
		gsl_vector_set(dphases,i,gsl_vector_get(dphases,i)-gsl_vector_get(dphase,i));

	/* Ignore correction when incr. variation is < CUTOFF, CUTOFF = M_PI */
	for(i=0;i<phase->size;i++) {
		if(fabs(gsl_vector_get(dphase,i)) < M_PI)
			gsl_vector_set(dphases,i,0);
	}

	/* Integrate corrections */
	for(i=1;i<phase->size;i++)
		gsl_vector_set(dphases,i,gsl_vector_get(dphases,i-1) + gsl_vector_get(dphases,i));

	/* Add correction to phase */
	for(i=0;i<phase->size-1;i++)
		gsl_vector_set(phase,i+1,gsl_vector_get(phase,i+1) + gsl_vector_get(dphases,i));

	/* Free memory */
	gsl_vector_free(dphase);
	gsl_vector_free(dphases);
}





/**
 * Function Unwrap2
 *
 * Fine-tune unwrapped phase angles.
 *
 * @param phase vector of phase values
 * @param tol tolerance
 */
void Unwrap2(gsl_vector *phase, double tol) {

	int i,j,sign;
	double diff;

	if(gsl_vector_get(phase,0) < gsl_vector_get(phase,phase->size-1))
		sign = -1;
	else
		sign = 1;

	for(i=1;i<phase->size-1;i++) {
		diff = gsl_vector_get(phase,i) - gsl_vector_get(phase,i+1);
		if(fabs(diff) > tol) {
			if(sign*diff > 0) {
				for(j=i+1;j<phase->size;j++)
					gsl_vector_set(phase,j,gsl_vector_get(phase,j) + sign*M_PI);
			}
		}
	}
}









// TEST
void Create_highband_excitation(gsl_vector *excitation_highband, gsl_matrix *erbgain, PARAM *params) {

	/////////////////////////
	int FS_high = 44100;
	//int FS_high = 16000;
	/////////////////////////

	int i,j,frame_index,sample_index = 0,erb_channels = erbgain->size2,signal_length = excitation_highband->size;
	double r = FS_high/(double)params->FS;
	gsl_vector *noise;

	/* Start loop */
	while(sample_index < signal_length) {

		frame_index = floor(params->n_frames*(sample_index/(double)signal_length));
		if(frame_index > params->n_frames-1)
			frame_index = params->n_frames-1;

		/* Allocate noise vector */
		noise = gsl_vector_calloc(rint(2.0*r*params->shift/params->speed));




		/**********************************/
		/*           Create noise         */
		/**********************************/

		/* Initialize */
		int n = noise->size;
		int n2 = ceil((n-1)/2.0);
		gsl_vector *real = gsl_vector_calloc(n);
		gsl_vector *imag = gsl_vector_calloc(n);
		gsl_vector *gain_values = gsl_vector_alloc(n2);
		gsl_vector *erb_inds = gsl_vector_alloc(n2);

		/* Evaluate ERB scale indices */
		for(i=0;i<n2;i++)
			gsl_vector_set(erb_inds,i,log10(0.00437*(i/(n2-1.0)*(FS_high/2.0))+1.0)/log10(0.00437*(FS_high/2.0)+1.0)*(erb_channels-SMALL_NUMBER));

		/* Evaluate values according to ERB rate */
		for(i=0;i<n2;i++) {
			j = floor(gsl_vector_get(erb_inds,i));
			gsl_vector_set(gain_values,i,gsl_matrix_get(erbgain,frame_index,j));
		}

		/* Convert HNR values from logarithmic scale to linear scale (actual noise amplitudes) */
		for(i=0;i<n2;i++)
			gsl_vector_set(gain_values,i,pow(10,(gsl_vector_get(gain_values,i)/20.0)));

		/* Modify both the magnitude and the phase of the spectrum  */
		double lfl = 0.3628;
		for(i=rint(lfl*n2);i<n2;i++) {
			gsl_vector_set(imag,i+1, RAND()*gsl_vector_get(gain_values,i));
			gsl_vector_set(real,i+1, RAND()*gsl_vector_get(gain_values,i));
		}

		/* Copy the noise to the folded spectrum as well */
		if(imag->size%2 == 0) {
			for(i=0;i<n2;i++) {
				gsl_vector_set(imag,n2+i,-gsl_vector_get(imag,n2-i));
				gsl_vector_set(real,n2+i,gsl_vector_get(real,n2-i));
			}
		} else {
			for(i=0;i<n2;i++) {
				gsl_vector_set(imag,n2+1+i,-gsl_vector_get(imag,n2-i));
				gsl_vector_set(real,n2+1+i,gsl_vector_get(real,n2-i));
			}
		}

		/* Create halfcomplex data */
		double data[n];
		data[0] = gsl_vector_get(real,0);
		for(i=1;i<n;i++) {
			if(i%2 == 1)
				data[i] = gsl_vector_get(real,(i+1)/2);
			else
				data[i] = gsl_vector_get(imag,i/2);
		}

		/* Inverse FFT */
		gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(n);
		gsl_fft_halfcomplex_wavetable *whc = gsl_fft_halfcomplex_wavetable_alloc(n);
		gsl_fft_halfcomplex_inverse(data,1,n,whc,work);
		gsl_fft_halfcomplex_wavetable_free(whc);
		gsl_fft_real_workspace_free(work);

		/* Set data from array "data" to vector "noise" */
		for(i=0;i<n;i++)
			gsl_vector_set(noise,i,data[i]);

		/* Free memory */
		gsl_vector_free(real);
		gsl_vector_free(imag);
		gsl_vector_free(gain_values);
		gsl_vector_free(erb_inds);

		/**********************************/
		/*      End of create noise       */
		/**********************************/


		/* Windown noise */
		for(i=0;i<n;i++)
			gsl_vector_set(noise,i,gsl_vector_get(noise,i)*HANN(i,n));

		/* Set noise to excitation */
		int q = rint(0.25*n);
		int ind;
		for(i=0;i<noise->size;i++) {
			ind = GSL_MIN(GSL_MAX(sample_index+i-q,0),excitation_highband->size-1);
			gsl_vector_set(excitation_highband,ind,gsl_vector_get(excitation_highband,ind) + gsl_vector_get(noise,i));
		}

		/* Free memory */
		gsl_vector_free(noise);

		/* Increment sample index */
		sample_index += rint(r*params->shift/params->speed);
	}
}








































/*********************************************************************/
/*                       TEST FUNCTIONS                              */
/*********************************************************************/



// TEST FUNCTION
// PRINT VECTOR TO FILE: p1.dat
void VPrint1(gsl_vector *vector) {
	FILE *f = fopen("p1.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p2.dat
void VPrint2(gsl_vector *vector) {
	FILE *f = fopen("p2.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p3.dat
void VPrint3(gsl_vector *vector) {
	FILE *f = fopen("p3.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p4.dat
void VPrint4(gsl_vector *vector) {
	FILE *f = fopen("p4.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p5.dat
void VPrint5(gsl_vector *vector) {
	FILE *f = fopen("p5.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p6.dat
void VPrint6(gsl_vector *vector) {
	FILE *f = fopen("p6.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p7.dat
void VPrint7(gsl_vector *vector) {
	FILE *f = fopen("p7.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p8.dat
void VPrint8(gsl_vector *vector) {
	FILE *f = fopen("p8.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p9.dat
void VPrint9(gsl_vector *vector) {
	FILE *f = fopen("p9.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p10.dat
void VPrint10(gsl_vector *vector) {
	FILE *f = fopen("p10.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p11.dat
void VPrint11(gsl_vector *vector) {
	FILE *f = fopen("p11.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p12.dat
void VPrint12(gsl_vector *vector) {
	FILE *f = fopen("p12.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: m1.dat
void MPrint1(gsl_matrix *matrix) {
	FILE *f = fopen("m1.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: m2.dat
void MPrint2(gsl_matrix *matrix) {
	FILE *f = fopen("m2.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: m3.dat
void MPrint3(gsl_matrix *matrix) {
	FILE *f = fopen("m3.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: m4.dat
void MPrint4(gsl_matrix *matrix) {
	FILE *f = fopen("m4.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: m5.dat
void MPrint5(gsl_matrix *matrix) {
	FILE *f = fopen("m5.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.15f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: m6.dat
void MPrint6(gsl_matrix *matrix) {
	FILE *f = fopen("m6.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.15f");
	fclose(f);
}


// TEST FUNCTION
// PRINT ARRAY TO FILE: a1.dat
void APrint1(double *array, int size) {
	int i;
	gsl_vector *vector = gsl_vector_alloc(size);
	for(i=0;i<size;i++)
		gsl_vector_set(vector,i,array[i]);
	FILE *f = fopen("a1.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
	gsl_vector_free(vector);
}

// TEST FUNCTION
// PRINT ARRAY TO FILE: a2.dat
void APrint2(double *array, int size) {
	int i;
	gsl_vector *vector = gsl_vector_alloc(size);
	for(i=0;i<size;i++)
		gsl_vector_set(vector,i,array[i]);
	FILE *f = fopen("a2.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
	gsl_vector_free(vector);
}

// TEST FUNCTION
// PRINT ARRAY TO FILE: a3.dat
void APrint3(double *array, int size) {
	int i;
	gsl_vector *vector = gsl_vector_alloc(size);
	for(i=0;i<size;i++)
		gsl_vector_set(vector,i,array[i]);
	FILE *f = fopen("a3.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
	gsl_vector_free(vector);
}

// TEST FUNCTION
// PRINT ARRAY TO FILE: a4.dat
void APrint4(double *array, int size) {
	int i;
	gsl_vector *vector = gsl_vector_alloc(size);
	for(i=0;i<size;i++)
		gsl_vector_set(vector,i,array[i]);
	FILE *f = fopen("a4.dat", "w");
	gsl_vector_fprintf(f, vector, "%.15f");
	fclose(f);
	gsl_vector_free(vector);
}

// TEST FUNCTION
// PRINT ELAPSED TIME
void TimePrint(double start_time) {
	printf("\n\nElapsed time: %1.2lf seconds.\n\n",((double)clock()-start_time)/(double)CLOCKS_PER_SEC);
}

// TEST FUNCTION
// PAUSE
void pause(int print) {
	if(print == 1)
		printf("\n\nPAUSED - PRESS ENTER TO CONTINUE\n");
	getchar();
}

// TEST FUNCTION
// TEST FOR NANS AND INFS (MATRIX)
int Find_matrix_NaN_Inf(gsl_matrix *m) {

	int i,j,err = 0;
	for(i=0;i<m->size1;i++) {
		for(j=0;j<m->size2;j++) {
			if(isnan(gsl_matrix_get(m,i,j)) == 1) {
				printf("Warning: NaN value found!\n");
				err = 1;
			}
			if(isinf(gsl_matrix_get(m,i,j)) == 1) {
				printf("Warning: Inf value found!\n");
				err = 1;
			}
		}
	}
	return err;
}

// TEST FUNCTION
// TEST FOR NANS AND INFS (VECTOR)
int Find_vector_NaN_Inf(gsl_vector *v) {

	int i,err = 0;
	for(i=0;i<v->size;i++) {
		if(isnan(gsl_vector_get(v,i)) == 1) {
			printf("Warning: NaN value found!\n");
			err = 1;
		}
		if(isinf(gsl_vector_get(v,i)) == 1) {
			printf("Warning: Inf value found!\n");
			err = 1;
		}
	}
	return err;
}

