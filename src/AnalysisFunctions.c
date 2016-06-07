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
 *              GlottHMM Speech Parameter Extractor
 * <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
 *
 * This program reads a speech file and extracts speech
 * parameters using glottal inverse filtering.
 *
 * This program has been written in Aalto University,
 * Department of Signal Processign and Acoustics, Espoo, Finland
 *
 * Author: Tuomo Raitio
 * Acknowledgements: Antti Suni, Paavo Alku, Martti Vainio
 *
 * File AnalysisFunctions.c
 * Version: 1.1
 *





/***********************************************/
/*                 INCLUDE                     */
/***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sndfile.h> 				/* Read and write wav */
#include <gsl/gsl_vector.h>		/* GNU, Vector */
#include <gsl/gsl_matrix.h>		/* GNU, Matrix */
#include <gsl/gsl_fft_real.h>		/* GSL, FFT */
#include <gsl/gsl_permutation.h>	/* GSL, Permutations */
#include <gsl/gsl_linalg.h>		/* GSL, Linear algebra */
#include <gsl/gsl_spline.h>		/* GSL, Interpolation */
#include <gsl/gsl_errno.h>			/* GSL, Error handling */
#include <gsl/gsl_poly.h>			/* GSL, Polynomials */
#include <gsl/gsl_sort_double.h>	/* GSL, Sort double */
#include <gsl/gsl_sort_vector.h>	/* GSL, Sort vector */
#include <gsl/gsl_complex.h>		/* GSL, Complex numbers */
#include <gsl/gsl_complex_math.h>	/* GSL, Arithmetic operations for complex numbers */
#include <gsl/gsl_fft_halfcomplex.h>/* GSL, FFT halfcomplex */
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fit.h>			/* GSL, Linear regression */
#include <gsl/gsl_multifit.h>		/* GSL, Higher order linear fitting */
#include <libconfig.h>				/* Configuration file */
#include "AnalysisFunctions.h"










/**
 * Function Check_command_line
 *
 * Check command line format and print instructions
 *
 * @param argc number of input arguments
 */
int Check_command_line(int argc) {

	/* Check command line format */
	if (argc < 3 || argc > 4) {

		printf("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
		printf("            GlottHMM - Speech Parameter Extractor (%s)\n",VERSION);
		printf("<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n\n");
		printf("Description:\n\n");
		printf("    Extraction of speech signal into vocal tract filter and voice\n");
		printf("    source parameters using glottal inverse filtering.\n\n");
		printf("Usage:\n\n");
		printf("    Analysis wav_file config_default config_user\n\n");
		printf("	wav_file        - Name of the audio file to be analysed\n");
		printf("	config_default  - Name of the default config file\n");
		printf("	config_user     - Name of the user config file (OPTIONAL) \n\n");
		printf("Version:\n\n");
		printf("    %s (%s)\n\n",VERSION,DATE);
		return EXIT_FAILURE;
	} else {

		/* Print program description */
		printf("\n<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>\n");
		printf("              Speech Parameter Extractor (%s)\n",VERSION);
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
 * @conf_type type of config file (def/usr)
 */
int Assign_config_parameters(struct config_t *conf, PARAM *params, int conf_type) {

	long int ival;
	double fval;
	const char *charval;
	int bool;

	/* Assign config parameters */
	if (config_lookup(conf, SAMPLING_FREQUENCY) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", SAMPLING_FREQUENCY);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, SAMPLING_FREQUENCY) != NULL) {
		config_lookup_int(conf, SAMPLING_FREQUENCY,&ival);
		params->FS = (int)ival;
	}
	if (config_lookup(conf, FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, FRAME_LENGTH) != NULL) {
		config_lookup_float(conf, FRAME_LENGTH,&fval);
		params->frame_length_ms = fval;
	}
	if (config_lookup(conf, UNVOICED_FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", UNVOICED_FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, UNVOICED_FRAME_LENGTH) != NULL) {
		config_lookup_float(conf, UNVOICED_FRAME_LENGTH,&fval);
		params->unvoiced_frame_length_ms = fval;
	}
	if (config_lookup(conf, FRAME_SHIFT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", FRAME_SHIFT);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, FRAME_SHIFT) != NULL) {
		config_lookup_float(conf, FRAME_SHIFT,&fval);
		params->shift_ms = fval;
	}
	if (config_lookup(conf, LPC_ORDER) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", LPC_ORDER);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LPC_ORDER) != NULL) {
		config_lookup_int(conf, LPC_ORDER,&ival);
		params->lpc_order_vt = (int)ival;
	}
	if (config_lookup(conf, LPC_ORDER_SOURCE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", LPC_ORDER_SOURCE);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LPC_ORDER_SOURCE) != NULL) {
		config_lookup_int(conf, LPC_ORDER_SOURCE,&ival);
		params->lpc_order_gl = (int)ival;
	}
	if (config_lookup(conf, DIFFERENTIAL_LSF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", DIFFERENTIAL_LSF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, DIFFERENTIAL_LSF) != NULL) {
		config_lookup_bool(conf, DIFFERENTIAL_LSF,&bool);
		params->differential_lsf = bool;
	}
	if (config_lookup(conf, WARPING_VT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", WARPING_VT);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, WARPING_VT) != NULL) {
		config_lookup_float(conf, WARPING_VT,&fval);
		params->lambda_vt = fval;
	}
	if (config_lookup(conf, WARPING_GL) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", WARPING_GL);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, WARPING_GL) != NULL) {
		config_lookup_float(conf, WARPING_GL,&fval);
		params->lambda_gl = fval;
	}
	if (config_lookup(conf, VOICING_THRESHOLD) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", VOICING_THRESHOLD);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, VOICING_THRESHOLD) != NULL) {
		config_lookup_float(conf, VOICING_THRESHOLD,&fval);
		params->voicing_threshold = fval;
	}
	if (config_lookup(conf, ZCR_THRESHOLD) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", ZCR_THRESHOLD);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, ZCR_THRESHOLD) != NULL) {
		config_lookup_float(conf, ZCR_THRESHOLD,&fval);
		params->ZCR = fval;
	}
	if (config_lookup(conf, F0_MIN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", F0_MIN);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, F0_MIN) != NULL) {
		config_lookup_float(conf, F0_MIN,&fval);
		params->fmin = fval;
	}
	if (config_lookup(conf, F0_MAX) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", F0_MAX);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, F0_MAX) != NULL) {
		config_lookup_float(conf, F0_MAX,&fval);
		params->fmax = fval;
	}
	if (config_lookup(conf, USE_F0_POSTPROCESSING) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_F0_POSTPROCESSING);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_F0_POSTPROCESSING) != NULL) {
		config_lookup_bool(conf, USE_F0_POSTPROCESSING,&bool);
		params->f0_postprocessing = bool;
	}
	if (config_lookup(conf, MAX_NUMBER_OF_PULSES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", MAX_NUMBER_OF_PULSES);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, MAX_NUMBER_OF_PULSES) != NULL) {
		config_lookup_int(conf, MAX_NUMBER_OF_PULSES,&ival);
		params->maxnumberofpulses = (int)ival;
	}
	if (config_lookup(conf, PULSEMAXLEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", PULSEMAXLEN);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, PULSEMAXLEN) != NULL) {
		config_lookup_float(conf, PULSEMAXLEN,&fval);
		params->pulsemaxlen_ms = fval;
	}
	if (config_lookup(conf, RESAMPLED_PULSELEN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", RESAMPLED_PULSELEN);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, RESAMPLED_PULSELEN) != NULL) {
		config_lookup_float(conf, RESAMPLED_PULSELEN,&fval);
		params->rspulsemaxlen_ms = fval;
	}
	if (config_lookup(conf, WAVEFORM_SAMPLES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", WAVEFORM_SAMPLES);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, WAVEFORM_SAMPLES) != NULL) {
		config_lookup_int(conf, WAVEFORM_SAMPLES,&ival);
		params->waveform_samples = (int)ival;
	}
	if (config_lookup(conf, MAX_PULSE_LEN_DIFF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", MAX_PULSE_LEN_DIFF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, MAX_PULSE_LEN_DIFF) != NULL) {
		config_lookup_float(conf, MAX_PULSE_LEN_DIFF,&fval);
		params->max_pulse_len_diff = fval;
	}
	if (config_lookup(conf, INVERT_SIGNAL) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", INVERT_SIGNAL);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, INVERT_SIGNAL) != NULL) {
		config_lookup_bool(conf, INVERT_SIGNAL,&bool);
		params->invert_signal = bool;
	}
	if (config_lookup(conf, EXTRACT_PULSELIB) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_PULSELIB);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_PULSELIB) != NULL) {
		config_lookup_bool(conf, EXTRACT_PULSELIB,&bool);
		params->extract_pulselib_params = bool;
	}
	if (config_lookup(conf, PITCH_SYNCHRONOUS_ANALYSIS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", PITCH_SYNCHRONOUS_ANALYSIS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, PITCH_SYNCHRONOUS_ANALYSIS) != NULL) {
		config_lookup_bool(conf, PITCH_SYNCHRONOUS_ANALYSIS,&bool);
		params->pitch_synchronous_analysis = bool;
	}
	if (config_lookup(conf, F0_CHECK_RANGE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", F0_CHECK_RANGE);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, F0_CHECK_RANGE) != NULL) {
		config_lookup_int(conf, F0_CHECK_RANGE,&ival);
		params->f0_check_range = (int)ival;
	}
	if (config_lookup(conf, RELATIVE_F0_THRESHOLD) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", RELATIVE_F0_THRESHOLD);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, RELATIVE_F0_THRESHOLD) != NULL) {
		config_lookup_float(conf, RELATIVE_F0_THRESHOLD,&fval);
		params->relative_f0_threshold = fval;
	}
	if (config_lookup(conf, F0_FRAME_LENGTH) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", F0_FRAME_LENGTH);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, F0_FRAME_LENGTH) != NULL) {
		config_lookup_float(conf, F0_FRAME_LENGTH,&fval);
		params->f0_frame_length_ms = fval;
	}
	if (config_lookup(conf, LP_STABILIZED) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", LP_STABILIZED);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LP_STABILIZED) != NULL) {
		config_lookup_bool(conf, LP_STABILIZED,&bool);
		params->lp_stabilized = bool;
	}
	if (config_lookup(conf, USE_IAIF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_IAIF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_IAIF) != NULL) {
		config_lookup_bool(conf, USE_IAIF,&bool);
		params->use_iaif = bool;
	}
	if (config_lookup(conf, LPC_ORDER_GL_IAIF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", LPC_ORDER_GL_IAIF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LPC_ORDER_GL_IAIF) != NULL) {
		config_lookup_int(conf, LPC_ORDER_GL_IAIF,&ival);
		params->lpc_order_gl_iaif = (int)ival;
	}
	if (config_lookup(conf, USE_MOD_IAIF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_MOD_IAIF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_MOD_IAIF) != NULL) {
		config_lookup_bool(conf, USE_MOD_IAIF,&bool);
		params->use_mod_iaif = bool;
	}
	if (config_lookup(conf, HNR_CHANNELS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", HNR_CHANNELS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, HNR_CHANNELS) != NULL) {
		config_lookup_int(conf, HNR_CHANNELS,&ival);
		params->hnr_channels = (int)ival;
	}
	if (config_lookup(conf, NUMBER_OF_HARMONICS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", NUMBER_OF_HARMONICS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, NUMBER_OF_HARMONICS) != NULL) {
		config_lookup_int(conf, NUMBER_OF_HARMONICS,&ival);
		params->number_of_harmonics = (int)ival;
	}
	if (config_lookup(conf, FORMANT_PRE_ENH_LPC_DELTA) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", FORMANT_PRE_ENH_LPC_DELTA);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, FORMANT_PRE_ENH_LPC_DELTA) != NULL) {
		config_lookup_float(conf, FORMANT_PRE_ENH_LPC_DELTA,&fval);
		params->formant_enh_lpc_delta = fval;
	}
	if (config_lookup(conf, FORMANT_PRE_ENH_COEFF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", FORMANT_PRE_ENH_COEFF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, FORMANT_PRE_ENH_COEFF) != NULL) {
		config_lookup_float(conf, FORMANT_PRE_ENH_COEFF,&fval);
		params->formant_enh_coeff = fval;
	}
	if (config_lookup(conf, SEPARATE_VOICED_UNVOICED_SPECTRUM) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", SEPARATE_VOICED_UNVOICED_SPECTRUM);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, SEPARATE_VOICED_UNVOICED_SPECTRUM) != NULL) {
		config_lookup_bool(conf, SEPARATE_VOICED_UNVOICED_SPECTRUM,&bool);
		params->sep_vuv_spectrum = bool;
	}
	if (config_lookup(conf, EXTRACT_SOURCE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_SOURCE);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_SOURCE) != NULL) {
		config_lookup_bool(conf, EXTRACT_SOURCE,&bool);
		params->extract_source = bool;
	}
	if (config_lookup(conf, HP_FILTERING) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", HP_FILTERING);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, HP_FILTERING) != NULL) {
		config_lookup_bool(conf, HP_FILTERING,&bool);
		params->hp_filtering = bool;
	}
	if (config_lookup(conf, EXTRACT_F0) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_F0);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_F0) != NULL) {
		config_lookup_bool(conf, EXTRACT_F0,&bool);
		params->extract_f0 = bool;
	}
	if (config_lookup(conf, EXTRACT_GAIN) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_GAIN);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_GAIN) != NULL) {
		config_lookup_bool(conf, EXTRACT_GAIN,&bool);
		params->extract_gain = bool;
	}
	if (config_lookup(conf, EXTRACT_LSF) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_LSF);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_LSF) != NULL) {
		config_lookup_bool(conf, EXTRACT_LSF,&bool);
		params->extract_lsf = bool;
	}
	if (config_lookup(conf, EXTRACT_LSFSOURCE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_LSFSOURCE);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_LSFSOURCE) != NULL) {
		config_lookup_bool(conf, EXTRACT_LSFSOURCE,&bool);
		params->extract_tilt = bool;
	}
	if (config_lookup(conf, EXTRACT_HNR) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_HNR);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_HNR) != NULL) {
		config_lookup_bool(conf, EXTRACT_HNR,&bool);
		params->extract_hnr = bool;
	}
	if (config_lookup(conf, EXTRACT_HARMONICS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_HARMONICS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_HARMONICS) != NULL) {
		config_lookup_bool(conf, EXTRACT_HARMONICS,&bool);
		params->extract_harmonics = bool;
	}
	if (config_lookup(conf, EXTRACT_H1H2) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_H1H2);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_H1H2) != NULL) {
		config_lookup_bool(conf, EXTRACT_H1H2,&bool);
		params->extract_h1h2 = bool;
	}
	if (config_lookup(conf, EXTRACT_NAQ) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_NAQ);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_NAQ) != NULL) {
		config_lookup_bool(conf, EXTRACT_NAQ,&bool);
		params->extract_naq = bool;
	}
	if (config_lookup(conf, EXTRACT_WAVEFORM) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_WAVEFORM);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_WAVEFORM) != NULL) {
		config_lookup_bool(conf, EXTRACT_WAVEFORM,&bool);
		params->extract_waveform = bool;
	}
	if (config_lookup(conf, EXTRACT_INFOFILE) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_INFOFILE);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_INFOFILE) != NULL) {
		config_lookup_bool(conf, EXTRACT_INFOFILE,&bool);
		params->extract_info = bool;
	}
	if (config_lookup(conf, USE_EXTERNAL_F0) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", USE_EXTERNAL_F0);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, USE_EXTERNAL_F0) != NULL) {
		config_lookup_bool(conf, USE_EXTERNAL_F0,&bool);
		params->use_external_f0 = bool;
	}
	if (config_lookup(conf, NOISE_REDUCTION_ANALYSIS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", NOISE_REDUCTION_ANALYSIS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, NOISE_REDUCTION_ANALYSIS) != NULL) {
		config_lookup_bool(conf, NOISE_REDUCTION_ANALYSIS,&bool);
		params->noise_reduction_analysis = bool;
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
	if (config_lookup(conf, LOG_F0) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", LOG_F0);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LOG_F0) != NULL) {
		config_lookup_bool(conf, LOG_F0,&bool);
		params->logf0 = bool;
	}
	if (config_lookup(conf, EXTRACT_ONLY_UNIQUE_PULSES) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_ONLY_UNIQUE_PULSES);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_ONLY_UNIQUE_PULSES) != NULL) {
		config_lookup_bool(conf, EXTRACT_ONLY_UNIQUE_PULSES,&bool);
		params->extract_only_unique_pulses = bool;
	}
	if (config_lookup(conf, EXTRACT_ONE_PULSE_PER_FRAME) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_ONE_PULSE_PER_FRAME);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_ONE_PULSE_PER_FRAME) != NULL) {
		config_lookup_bool(conf, EXTRACT_ONE_PULSE_PER_FRAME,&bool);
		params->extract_one_pulse_per_frame = bool;
	}
	if (config_lookup(conf, EXTRACT_FFT_SPECTRA) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", EXTRACT_FFT_SPECTRA);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTRACT_FFT_SPECTRA) != NULL) {
		config_lookup_bool(conf, EXTRACT_FFT_SPECTRA,&bool);
		params->write_fftspectra_to_file = bool;
	}
	if (config_lookup(conf, UNVOICED_PRE_EMPHASIS) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n", UNVOICED_PRE_EMPHASIS);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, UNVOICED_PRE_EMPHASIS) != NULL) {
		config_lookup_bool(conf, UNVOICED_PRE_EMPHASIS,&bool);
		params->unvoiced_pre_emphasis = bool;
	}

	/* Data format */
	if (config_lookup(conf, DATA_FORMAT) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n",DATA_FORMAT);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, DATA_FORMAT) != NULL) {
		config_lookup_string(conf,DATA_FORMAT,&charval);
		if(strcmp(charval,DATA_FORMAT_ASCII) == 0)
			params->data_format = DATA_FORMAT_ID_ASCII;
		else if(strcmp(charval,DATA_FORMAT_BINARY) == 0)
			params->data_format = DATA_FORMAT_ID_BINARY;
		else {
			printf("\nError: Invalid configuration value \"%s\".\n", DATA_FORMAT);
			return EXIT_FAILURE;
		}
	}

	/* LP method */
	if (config_lookup(conf, LP_METHOD) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n",LP_METHOD);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LP_METHOD) != NULL) {
		config_lookup_string(conf,LP_METHOD,&charval);
		if(strcmp(charval,LP_METHOD_LPC) == 0)
			params->lp_method = LP_METHOD_ID_LPC;
		else if(strcmp(charval,LP_METHOD_WLP) == 0)
			params->lp_method = LP_METHOD_ID_WLP;
		else if(strcmp(charval,LP_METHOD_XLP) == 0)
			params->lp_method = LP_METHOD_ID_XLP;
		else {
			printf("\nError: Invalid configuration value \"%s\".\n", LP_METHOD);
			return EXIT_FAILURE;
		}
	}

	/* LP weighting method */
	if (config_lookup(conf, LP_WEIGHTING) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n",LP_WEIGHTING);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, LP_WEIGHTING) != NULL) {
		config_lookup_string(conf,LP_WEIGHTING,&charval);
		if(strcmp(charval,LP_WEIGHTING_STE) == 0)
			params->lp_weighting = LP_WEIGHTING_ID_STE;
		else if(strcmp(charval,LP_WEIGHTING_GCI) == 0)
			params->lp_weighting = LP_WEIGHTING_ID_GCI;
		else {
			printf("\nError: Invalid configuration value \"%s\".\n", LP_WEIGHTING);
			return EXIT_FAILURE;
		}
	}

	/* Formant enhancement method */
	if (config_lookup(conf, FORMANT_PRE_ENH_METHOD) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n",FORMANT_PRE_ENH_METHOD);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, FORMANT_PRE_ENH_METHOD) != NULL) {
		config_lookup_string(conf,FORMANT_PRE_ENH_METHOD,&charval);
		if(strcmp(charval,FORMANT_ENH_LSF) == 0)
			params->formant_enh_method = FORMANT_ENH_ID_LSF;
		else if(strcmp(charval,FORMANT_ENH_LPC) == 0)
			params->formant_enh_method = FORMANT_ENH_ID_LPC;
		else if(strcmp(charval,FORMANT_ENH_NONE) == 0)
			params->formant_enh_method = FORMANT_ENH_ID_NONE;
		else {
			printf("\nError: Invalid configuration value \"%s\".\n", FORMANT_PRE_ENH_METHOD);
			return EXIT_FAILURE;
		}
	}

	/* High-pass filter filename */
	if (config_lookup(conf, HPFILTER_FILENAME) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n",HPFILTER_FILENAME);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, HPFILTER_FILENAME) != NULL) {
		if(conf_type == USR_CONF) // Free default
			free(params->hpfilter_filename);
		config_lookup_string(conf,HPFILTER_FILENAME,&charval);
		params->hpfilter_filename = (char *)malloc((strlen(charval)+1)*sizeof(char));
		strcpy(params->hpfilter_filename,charval);
	}

	/* External F0 filename */
	if (config_lookup(conf, EXTERNAL_F0_FILENAME) == NULL && conf_type == DEF_CONF) {
		printf("\nError: Could not read default configuration \"%s\".\n",EXTERNAL_F0_FILENAME);
		return EXIT_FAILURE;
	} else if (config_lookup(conf, EXTERNAL_F0_FILENAME) != NULL) {
		if(conf_type == USR_CONF) // Free default
			free(params->external_f0_filename);
		config_lookup_string(conf,EXTERNAL_F0_FILENAME,&charval);
		params->external_f0_filename = (char *)malloc((strlen(charval)+1)*sizeof(char));
		strcpy(params->external_f0_filename,charval);
	}

	/* Convert milliseconds to samples */
	params->frame_length = rint(params->FS*params->frame_length_ms/1000);
	params->shift = rint(params->FS*params->shift_ms/1000);
	params->pulsemaxlen = rint(params->FS*params->pulsemaxlen_ms/1000);
	params->rspulsemaxlen = rint(params->FS*params->rspulsemaxlen_ms/1000);
	params->f0_frame_length = rint(params->FS*params->f0_frame_length_ms/1000);
	params->unvoiced_frame_length = rint(params->FS*params->unvoiced_frame_length_ms/1000);

	/* Free memory */
	config_destroy(conf);
	free(conf);

	return EXIT_SUCCESS;
}

















/**
 * Function Check_parameter_validity
 *
 * Check the validity of parameters
 *
 * @param filename name of the configuration file
 */
int Check_parameter_validity(PARAM *params) {

	if(params->lpc_order_gl < 3) {
		printf("\nError: The degree of the source spectral model is too low!\n\n");return EXIT_FAILURE;}
	if(params->f0_frame_length_ms < 1000/params->fmin) {
		printf("\nError: F0 frame length is too short with respect to the lowest F0 value!\n\n");return EXIT_FAILURE;}
	if(params->lambda_vt < -1 || params->lambda_vt > 1) {
		printf("\nError: Vocal tract warping coefficient must be a value in the range [-1,1]!\n\n");return EXIT_FAILURE;}
	if(params->lambda_gl < -1 || params->lambda_gl > 1) {
		printf("\nError: Voice source warping coefficient must be a value in the range [-1,1]!\n\n");return EXIT_FAILURE;}
	if(params->f0_frame_length_ms/1000*params->FS > FFT_LENGTH_LONG) {
		printf("\nError: F0 frame length greater than FFT length!\n\n");return EXIT_FAILURE;}
	if(params->frame_length_ms/1000*params->FS > FFT_LENGTH) {
		printf("\nError: Frame length greater than FFT length!\n\n");return EXIT_FAILURE;}
	if(params->unvoiced_frame_length > params->frame_length) {
		printf("\nError: Unvoiced frame length is greater than overall frame length!\n\n");return EXIT_FAILURE;}
	return EXIT_SUCCESS;
}







/**
 * Function Read_soundfile
 *
 * Read soundfile to vector. Also specify sampling frequency, number of frames, and signal length.
 * Invert signal if required.
 *
 * @param filename
 * @return signal
 */
gsl_vector *Read_soundfile(char *filename, PARAM *params) {

	/* Initialize */
	int i;
	SNDFILE *soundfile;
	SF_INFO sfinfo;
	sfinfo.frames = 0;
	sfinfo.samplerate = 0;
	sfinfo.channels = 0;
	sfinfo.format = 0;
	sfinfo.sections = 0;
	sfinfo.seekable = 0;

	/* Open soundfile */
	soundfile = sf_open(filename, SFM_READ, &sfinfo);
	if(soundfile == NULL) {
		fprintf(stderr, "\nError: Cannot open file \"%s\".\n", filename);
		return NULL;
	}

	/* Check sampling frequency */
	if(params->FS != sfinfo.samplerate) {
		printf("\nError: Sampling frequency of the sound file is different than indicated in config file (%i Hz vs %i Hz)\n\n",(int)sfinfo.samplerate,params->FS);
		return NULL;
	}

	/* Count the number of frames that fits to the signal (round to next integer) */
	params->n_frames = ceil(sfinfo.frames/(double)params->shift-(double)params->frame_length/(double)params->shift)+1;
	params->signal_length = (params->n_frames + params->frame_length/(double)params->shift - 1)*(double)params->shift;

	/* Read sound file into vector */
	gsl_vector *signal = gsl_vector_calloc(params->signal_length);
	double *temp = (double *)calloc(sfinfo.frames, sizeof(double));
	sf_read_double(soundfile, temp, sfinfo.frames);
	for(i=0; i<sfinfo.frames; i++)
		gsl_vector_set(signal, i, temp[i]);
	sf_close(soundfile);
	free(temp);

	/* Invert signal if required */
	if(params->invert_signal == 1) {
		for(i=0;i<signal->size;i++)
			gsl_vector_set(signal,i,-gsl_vector_get(signal,i));
	}

	/* Return signal vector */
	return signal;
}











/**
 * Function HP_Filter
 *
 * High-pass filter signal
 *
 * @param signal pointer to the samples
 */
int HighPassFilter(gsl_vector *signal, PARAM *params) {

	/* Do not filter if not requested */
	if(params->hp_filtering == 0) {
		free(params->hpfilter_filename);
		return EXIT_SUCCESS;
	}

	/* Evaluate filter coefficient file length */
	int n = EvalFileLength(params->hpfilter_filename);
	if(n == -1)
		return EXIT_FAILURE;

	/* Open file */
	FILE *filter = fopen(params->hpfilter_filename, "r");
	if(!filter) {
		printf("Error opening file \"%s\": %s\n", params->hpfilter_filename, strerror(errno));
		return EXIT_FAILURE;
	}

	/* Initialize */
	int i,j;
	double sum;
	gsl_vector *temp = gsl_vector_calloc(signal->size);
	gsl_vector *coeffs = gsl_vector_alloc(n);

	/* Read coefficients from file */
	gsl_vector_fscanf(filter,coeffs);

	/* Filter signal */
	for(i=0; i<signal->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, n-1); j++) {
			sum += gsl_vector_get(signal,i-j)*gsl_vector_get(coeffs,j);
		}
		gsl_vector_set(temp, GSL_MAX(i-(n-1)/2,0), sum);
	}

	/* Copy "temp" samples to "signal" */
	for(i=0; i<signal->size; i++) {
		gsl_vector_set(signal, i, gsl_vector_get(temp, i));
	}

	/* Free memory */
	gsl_vector_free(temp);
	gsl_vector_free(coeffs);
	fclose(filter);
	free(params->hpfilter_filename);

	return EXIT_SUCCESS;
}





/**
 * Function Read_external_f0
 *
 * Read external F0 file.
 *
 * -Interpolate if required.
 * -Remove some samples in the beginning and in the end of the vector
 * due to the windowing scheme used in Analysis.
 *
 * @param fundf pointer to f0 vector
 * @param params
 */
int Read_external_f0(gsl_vector **fundf, PARAM *params) {

	/* Do not read if not requested, allocate normal f0 vector */
	if(params->use_external_f0 == 0) {
		*(fundf) = gsl_vector_alloc(params->n_frames);
		free(params->external_f0_filename);
		return EXIT_SUCCESS;
	}

	/* Evaluate external f0 file length */
	int n = EvalFileLength(params->external_f0_filename);
	if(n == -1)
		return EXIT_FAILURE;

	/* Open file */
	FILE *f0_file = fopen(params->external_f0_filename, "r");
	if(!f0_file) {
		printf("Error opening file \"%s\": %s\n", params->external_f0_filename, strerror(errno));
		return EXIT_FAILURE;
	}

	/* Read f0 from file */
	gsl_vector *f0_raw = gsl_vector_alloc(n);
	gsl_vector_fscanf(f0_file,f0_raw);

	/* Interpolate to a specific length (n_frames + 2*empty_frames */
	int empty_frames = rint((params->frame_length/(double)params->shift - 1)/2);
	gsl_vector *f0_new = gsl_vector_alloc(params->n_frames + 2*empty_frames);
	if(n != params->n_frames + 2*empty_frames) {
		printf("External F0 is not the same length (%i) as required (%i).\n",n,params->n_frames + 2*empty_frames);
		printf("Apply interpolation to get correct length!\n\n");
		Interpolate(f0_raw,f0_new);
		int i;
		for(i=0;i<f0_new->size;i++)  // Remove possible bad values
			if(gsl_vector_get(f0_new,i) < params->fmin)
				gsl_vector_set(f0_new,i,0);
	} else {
		gsl_vector_memcpy(f0_new,f0_raw);
	}

	/* Remove "empty frames" in the beginning and in the end of the vector */
	gsl_vector *f0 = gsl_vector_alloc(params->n_frames);
	int i;
	for(i=0;i<f0->size;i++)
		gsl_vector_set(f0,i,gsl_vector_get(f0_new,i+empty_frames));

	/* Set pointer */
	*(fundf) = f0;

	/* Free memory */
	fclose(f0_file);
	free(params->external_f0_filename);
	gsl_vector_free(f0_raw);
	gsl_vector_free(f0_new);

	return EXIT_SUCCESS;
}







/**
 * Function EvalFileLength
 *
 * Read file and count the number of lines
 *
 * @param name filename
 * @return number of parameters or -1 if file could not be opened
 */
int EvalFileLength(const char *name) {

	FILE *file;
	char s[DEF_STRING_LEN];
	int ind;

	/* Open file */
	file = fopen(name, "r");
	if(!file) {
		printf("Error opening file \"%s\": %s\n", name, strerror(errno));
		return -1;
	}

	/* Read lines until EOF */
	ind = 0;
	while(fscanf(file,"%s",s) != EOF)
		ind++;

	fclose(file);
	return ind;
}




/**
 * Function Print_progress
 *
 * Print progress of the analysis
 *
 * @param index time index
 * @param n_frames total number of frames
 * @param audio file name
 * @param FS sampling frequency
 */
void Print_progress(int index, int n_frames, char *filename, int FS) {

	int idx = 0;
	if(index == 0)
		printf("Extracting parameters from \"%s\" (%i Hz)\n\n", filename, FS);
	if(index % GSL_MAX((n_frames-1)/PROGRESS_UPDATE_INTERVAL,1) == 0) {
		idx = rint(100.0*index/(double)(n_frames-1));
		if(idx != 100)
			printf("%i %%\n", idx);
	}
	if(index == n_frames-1)
		printf("100 %%\n");
}










/**
 * Function Allocate_variables
 *
 * Allocate memory for variables
 *
 * @param ...
 */
void Allocate_variables(gsl_vector **frame, gsl_vector **frame0, gsl_vector **glottal, gsl_vector **gain, gsl_vector **uvgain,
		gsl_vector **f0_frame, gsl_vector **glottal_f0, gsl_vector **f0_frame0, gsl_vector **glottsig, gsl_vector **glottsig_f0, gsl_vector **source_signal,
		gsl_matrix **fundf_candidates, gsl_matrix **LSF, gsl_matrix **LSF2, gsl_matrix **bp_gain, gsl_matrix **spectral_tilt, gsl_matrix **HNR,
		gsl_matrix **waveform, gsl_matrix **harmonics, gsl_vector **h1h2, gsl_vector **naq, gsl_matrix **fftmatrix_vt, gsl_matrix **fftmatrix_src,
		gsl_matrix **fftmatrix_uv, PARAM *params) {

	/* Allocate vector */
	*(source_signal) = gsl_vector_calloc(params->signal_length);
	*(frame) = gsl_vector_alloc(params->frame_length);
	*(frame0) = gsl_vector_calloc(params->lpc_order_vt);
	*(glottal) = gsl_vector_alloc(params->frame_length);
	*(glottsig) = gsl_vector_alloc(params->frame_length);
	*(glottsig_f0) = gsl_vector_alloc(params->f0_frame_length);
	*(gain) = gsl_vector_alloc(params->n_frames);
	*(uvgain) = gsl_vector_alloc(params->n_frames);
	*(h1h2) = gsl_vector_calloc(params->n_frames);
	*(naq) = gsl_vector_calloc(params->n_frames);
	*(f0_frame) = gsl_vector_calloc(params->f0_frame_length);
	*(f0_frame0) = gsl_vector_calloc(params->lpc_order_vt);
	*(glottal_f0) = gsl_vector_alloc(params->f0_frame_length);

	/* Allocate matrices */
	*(fundf_candidates) = gsl_matrix_calloc(params->n_frames, NUMBER_OF_F0_CANDIDATES);
	*(spectral_tilt) = gsl_matrix_calloc(params->n_frames, params->lpc_order_gl);
	*(LSF) = gsl_matrix_calloc(params->n_frames, params->lpc_order_vt);
	*(LSF2) = gsl_matrix_calloc(params->n_frames, params->lpc_order_vt);
	*(HNR) = gsl_matrix_calloc(params->n_frames, params->hnr_channels);
	*(bp_gain) = gsl_matrix_calloc(params->n_frames, GAIN_FREQUENCY_BANDS);
	*(waveform) = gsl_matrix_calloc(params->n_frames,params->waveform_samples);
	*(harmonics) = gsl_matrix_calloc(params->n_frames,params->number_of_harmonics);
	*(fftmatrix_vt) = gsl_matrix_calloc(params->n_frames,OUTPUT_FFT_LENGTH/2);
	*(fftmatrix_src) = gsl_matrix_calloc(params->n_frames,OUTPUT_FFT_LENGTH/2);
	*(fftmatrix_uv) = gsl_matrix_calloc(params->n_frames,OUTPUT_FFT_LENGTH/2);
}









/**
 * Function Allocate_pulselib_variables
 *
 * Allocate memory for pulse library variables
 *
 * @param ...
 */
void Allocate_pulselib_variables(gsl_matrix **gpulses, gsl_matrix **gpulses_rs, gsl_vector **pulse_inds, gsl_vector **pulse_pos, gsl_vector **pulse_lengths, gsl_matrix **plsf,
		gsl_matrix **ptilt, gsl_matrix **pharm, gsl_matrix **phnr, gsl_matrix **pwaveform, gsl_vector **pgain, gsl_vector **ph1h2, gsl_vector **pnaq, PARAM *params) {

	/* Allocate parameters if used */
	if(params->extract_pulselib_params == 1) {

		/* If pulse library parameters are used, allocate memory */
		*(gpulses) = gsl_matrix_calloc(params->maxnumberofpulses,params->pulsemaxlen);
		*(gpulses_rs) = gsl_matrix_calloc(params->maxnumberofpulses,params->rspulsemaxlen);
		*(pulse_inds) = gsl_vector_calloc(params->maxnumberofpulses);
		*(pulse_pos) = gsl_vector_calloc(params->maxnumberofpulses);
		*(pulse_lengths) = gsl_vector_calloc(params->maxnumberofpulses);
		*(plsf) = gsl_matrix_calloc(params->maxnumberofpulses,params->lpc_order_vt);
		*(ptilt) = gsl_matrix_calloc(params->maxnumberofpulses,params->lpc_order_gl);
		*(pharm) = gsl_matrix_calloc(params->maxnumberofpulses,params->number_of_harmonics);
		*(phnr) = gsl_matrix_calloc(params->maxnumberofpulses,params->hnr_channels);
		*(pwaveform) = gsl_matrix_calloc(params->maxnumberofpulses,params->waveform_samples);
		*(pgain) = gsl_vector_calloc(params->maxnumberofpulses);
		*(ph1h2) = gsl_vector_calloc(params->maxnumberofpulses);
		*(pnaq) = gsl_vector_calloc(params->maxnumberofpulses);
		params->number_of_pulses = 0;

	} else {

		/* Set NULL if not used */
		params->number_of_pulses = 0;
		*(gpulses) = NULL;
		*(gpulses_rs) = NULL;
		*(pulse_inds) = NULL;
		*(pulse_pos) = NULL;
		*(pulse_lengths) = NULL;
		*(plsf) = NULL;
		*(ptilt) = NULL;
		*(pharm) = NULL;
		*(phnr) = NULL;
		*(pwaveform) = NULL;
		*(pgain) = NULL;
		*(ph1h2) = NULL;
		*(pnaq) = NULL;
	}
}









/**
 * Function Get_samples_to_frames
 *
 * Get samples to frames
 *
 * @param frame pointer to signal vector
 * @param frame pointer to main frame
 * @param frame0 pointer to pre-frame
 * @param frame pointer to f0 frame
 * @param frame0 pointer to f0 pre-frame
 * @param shift size
 * @parame index time
 */
void Get_samples_to_frames(gsl_vector *signal, gsl_vector *frame, gsl_vector *frame0, gsl_vector *f0_frame, gsl_vector *f0_frame0, int shift, int index) {

	/* Initialize */
	int i;
	int frame_length = frame->size;
	int f0_frame_length = f0_frame->size;
	int lpc_order_vt = frame0->size;
	int signal_length = signal->size;

	/* Get samples to frame */
	for(i=0; i<frame_length; i++)
		gsl_vector_set(frame, i, gsl_vector_get(signal, index*shift+i));

	/* Get pre-frame samples for smooth filtering */
	for(i=0; i<lpc_order_vt; i++) {
		if(index*shift+i-lpc_order_vt >= 0)
			gsl_vector_set(frame0, i, gsl_vector_get(signal, index*shift+i-lpc_order_vt));
		else
			gsl_vector_set(frame0, i, 0);
	}

	/* Get f0_frame and f0_frame0 (longer frame) */
	/* In the beginning of the signal */
	if(index*shift - ceil(f0_frame_length/2.0) + ceil(frame_length/2.0) < 0) {
		for(i=0; i<f0_frame_length; i++)
			gsl_vector_set(f0_frame, i, gsl_vector_get(signal, i));
		gsl_vector_set_zero(f0_frame0);
	/* In the end of the signal */
	} else if(index*shift + ceil(f0_frame_length/2.0) + ceil(frame_length/2.0) > signal_length-1) {
		for(i=0; i<f0_frame_length; i++)
			gsl_vector_set(f0_frame, i, gsl_vector_get(signal, signal->size-f0_frame_length+i));
		for(i=0; i<lpc_order_vt; i++)
			gsl_vector_set(f0_frame0, i, gsl_vector_get(signal, signal->size-f0_frame_length+i-lpc_order_vt));
	/* Between */
	} else {
		for(i=0; i<f0_frame_length; i++)
			gsl_vector_set(f0_frame, i, gsl_vector_get(signal, index*shift-ceil(f0_frame_length/2.0)+ceil(frame_length/2.0)+i));
		for(i=0; i<lpc_order_vt; i++)
			gsl_vector_set(f0_frame0, i, gsl_vector_get(signal, GSL_MAX(index*shift-ceil(f0_frame_length/2.0)+ceil(frame_length/2.0)+i-lpc_order_vt,0)));
	}
}








/**
 * Function Gain
 *
 * Calculate energy
 *
 * @param frame pointer to the samples
 * @param gain vector for results
 * @param index time index
 * @param windowing switch for using windowing
 */
void Gain(gsl_vector *frame, gsl_vector *gain, int index, int windowing) {

	int i;
	double sum;

	/* Windowing switch */
	if(windowing == 1) {

		/* Windowing */
		gsl_vector *wframe = gsl_vector_alloc(frame->size);
		for(i=0;i<frame->size;i++)
			gsl_vector_set(wframe,i,gsl_vector_get(frame,i)*HANN(i,frame->size));

		/* Evaluate gain of frame, normalize energy per sample basis */
		sum = 0;
		for(i=0;i<frame->size;i++) {
			sum = sum + gsl_vector_get(wframe,i)*gsl_vector_get(wframe,i);
		}
		gsl_vector_set(gain, index, 10.0*log10((8.0/3.0)*sum/E_REF/((double)(frame->size))));

		/* Free memory */
		gsl_vector_free(wframe);

	} else {

		/* Evaluate gain of frame, normalize energy per sample basis */
		sum = 0;
		for(i=0;i<frame->size;i++) {
			sum = sum + gsl_vector_get(frame,i)*gsl_vector_get(frame,i);
		}
		gsl_vector_set(gain, index, 10.0*log10(sum/E_REF/((double)(frame->size))));
	}

	/* Ensure non-infinity values */
	if(isinf(gsl_vector_get(gain,index)) != 0)
		gsl_vector_set(gain,index,MIN_LOG_POWER);
}





/**
 * Function UnvoicedGain
 *
 * Calculate energy from potentially unvoiced frames (shorter window)
 *
 * @param frame pointer to the samples
 * @param gain vector for results
 * @param index time index
 * @param windowing switch for using windowing
 * @param unvoiced_frame_length
 */
void UnvoicedGain(gsl_vector *frame, gsl_vector *gain, int index, int windowing, int unvoiced_frame_length) {

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
 * Function BandPassGain
 *
 * Calculate bandpass gains from Fourier transform
 *
 * @param frame pointer to samples
 * @param bp_gain values to be saved
 * @param params parameter structure
 * @param index frame index
 *
 */
void BandPassGain(gsl_vector *frame, gsl_matrix *bp_gain, PARAM *params, int index) {

	int i;
	double data[FFT_LENGTH] = {0};
	gsl_vector *fft = gsl_vector_calloc(FFT_LENGTH/2);

	/* FFT */
	for (i=0; i<frame->size; i++) {
		data[i] = gsl_vector_get(frame,i); // No windowing here
	}
	gsl_fft_real_radix2_transform(data, 1, FFT_LENGTH);
	for(i=1; i<FFT_LENGTH/2; i++) {
		gsl_vector_set(fft, i, sqrt(pow(data[i], 2) + pow(data[FFT_LENGTH-i], 2)));
	}
	gsl_vector_set(fft, 0, data[0]);

	/* Calculate gain for each band: 0--1, 1--2, 2--4, 4--6, 6--FS/2 kHz */
	int k,freq_lims[6] = {0,GSL_MIN(1000,params->FS/2),GSL_MIN(2000,params->FS/2),GSL_MIN(4000,params->FS/2),GSL_MIN(6000,params->FS/2),params->FS/2};
	double weights[5] = {1,1,0.5,0.5,1000/(double)(params->FS/2-freq_lims[4])};
	for(k=0; k<5; k++) {
		for(i=rint(freq_lims[k]/(double)params->FS*(double)FFT_LENGTH); i<rint(freq_lims[k+1]/(double)params->FS*(double)FFT_LENGTH); i++)
			gsl_matrix_set(bp_gain, index, k, gsl_matrix_get(bp_gain, index, k) + gsl_vector_get(fft, i));
		gsl_matrix_set(bp_gain,index,k,weights[k]*gsl_matrix_get(bp_gain,index,k));
	}
	gsl_vector_free(fft);
}







/**
 * Function InverseFiltering
 *
 * Glottal inverse filtering:
 * Estimate source signal and parametrize spectrum
 * Iterative Adaptive Inverse Filtering (IAIF)
 *
 * @param frame pointer to the samples
 * @param frame0 pointer to pre-frame samples
 * @param glottal pointer to the glottal source signal
 * @param LSF pointer to the voiced LSFs (glottal inverse filtered)
 * @param LSF2 pointer to the unvoiced LSFs (normal LPC)
 * @param spectral_tilt spectral tilt vector
 * @param index time index
 * @param lpc_order_vt LPC degree 1
 * @param LPC_g LPC degree 2
 * @param rho leaky integrator constant
 * @param lpc_order_gl spectral tilt degree
 */
void InverseFiltering(gsl_vector *frame,
						gsl_vector *frame0,
						gsl_vector *glottal,
						gsl_matrix *LSF,
						gsl_matrix *LSF2,
						gsl_matrix *spectral_tilt,
						gsl_vector *fundf,
						gsl_vector *glottsig,
						gsl_matrix *fftmatrix_vt,
						gsl_matrix *fftmatrix_src,
						gsl_matrix *fftmatrix_uv,
						int index,
						PARAM *params) {

	gsl_vector *a_l = gsl_vector_calloc(2);
	gsl_vector *a_p = gsl_vector_calloc(params->lpc_order_vt+1);
	gsl_vector *a_g = gsl_vector_calloc(params->lpc_order_gl+1);
	gsl_vector *a_g_iaif = gsl_vector_calloc(params->lpc_order_gl_iaif+1);
	gsl_vector *a_tilt = gsl_vector_calloc(params->lpc_order_gl+1);
	gsl_vector *frame_concat = gsl_vector_alloc(params->lpc_order_vt+frame->size);
	gsl_vector *B = gsl_vector_alloc(1);
	gsl_vector_set(B,0,1);
	int i;

	/* Concatenate frame and frame0 */
	for(i=0; i<params->lpc_order_vt; i++)
		gsl_vector_set(frame_concat, i, gsl_vector_get(frame0, i));
	for(i=params->lpc_order_vt; i<frame->size+params->lpc_order_vt; i++)
		gsl_vector_set(frame_concat, i, gsl_vector_get(frame, i-params->lpc_order_vt));

	/* Use glottal inverse filtering */
	if(params->use_iaif == 1) {

		/* 2. LPC-analysis (order l = 1) */
		/* 3. Inverse filtering */
		WLPC(frame, a_l, params->lambda_gl);
		WFilter(frame_concat, glottal, a_l, B, params->lambda_gl);
		Remove_mean(glottal);

		/* FFT of first VT estimate */
		EvalFFTSpectrum(glottal, fftmatrix_vt, index, params);

		/* 4. Spectral analysis, select (W)LPC/(S)WLP/(S)XLP */
		/* 5. Inverse filtering, warped only if WLPS is used */
		if(params->lp_method == LP_METHOD_ID_WLP) {			// SWLP-analysis (order p, M = p, lag = 1)
			SWLP(glottal,a_p,params->lpc_order_vt,1,params->lp_weighting,fundf,params->FS,index,glottsig,params->lp_stabilized);
			Filter(frame_concat,glottal,a_p);
		} else if(params->lp_method == LP_METHOD_ID_XLP) {	// SXLP-analysis (order p, w = p)
			SXLP(glottal,a_p,params->lpc_order_vt,params->lp_stabilized);
			Filter(frame_concat,glottal,a_p);
		} else {											// LPC-analysis (order p)
			WLPC(glottal,a_p,params->lambda_vt);
			WFilter(frame_concat,glottal,a_p,B,params->lambda_vt);
		}

		/* Convert LPC to LSF */
		if(params->use_mod_iaif == 1)
			Convert_matrix_to_LSF(LSF, a_p, index);

		/* If full IAIF is used */
		if(params->use_mod_iaif == 0) {

			/* 6. Integration */
			Integrator(glottal, LIP_RADIATION);
			Remove_mean(glottal);

			/* 7. WLPC-analysis (order g = lpc_order_gl_iaif) */
			/* 8. Inverse filtering */
			/* 9. Integration */
			WLPC(glottal, a_g_iaif, params->lambda_gl);
			WFilter(frame_concat, glottal, a_g_iaif, B, params->lambda_gl);
			Integrator(glottal, LIP_RADIATION);
			Remove_mean(glottal);

			/* FFT of final VT estimate */
			EvalFFTSpectrum(glottal, fftmatrix_vt, index, params);

			/* 10. Spectral analysis, select (W)LPC/(S)WLP/(S)XLP */
			/* 11. Inverse filtering, warped only if WLPC is used */
			if(params->lp_method == LP_METHOD_ID_WLP) {			// SWLP-analysis (order p, M = p, lag = 1)
				SWLP(glottal,a_p,params->lpc_order_vt,1,params->lp_weighting,fundf,params->FS,index,glottsig,params->lp_stabilized);
				Filter(frame_concat,glottal,a_p);
			} else if(params->lp_method == LP_METHOD_ID_XLP) {	// SXLP-analysis (order p, w = p)
				SXLP(glottal,a_p,params->lpc_order_vt,1);
				Filter(frame_concat,glottal,a_p);
			} else {											// LPC-analysis (order r = p), WLPC
				WLPC(glottal,a_p,params->lambda_vt);
				WFilter(frame_concat,glottal,a_p,B,params->lambda_vt);
			}

			/* Convert LPC to LSF */
			Convert_matrix_to_LSF(LSF, a_p, index);

		}

		/* 12. Integration */
		/* No integration if pre-emphasis is used */
		if(USE_PRE_EMPH == 0)
			Integrator(glottal, LIP_RADIATION);
		Remove_mean(glottal);

	} else {

		/* Pre-emphasis */
		Differentiate(frame,LIP_RADIATION);

		/* FFT of VT estimate */
		EvalFFTSpectrum(frame, fftmatrix_vt, index, params);

		/* If glottal inverse filtering is not used, use (W)LPC/(S)WLP/(S)XLP only once */
		if(params->lp_method == LP_METHOD_ID_WLP) {
			SWLP(frame,a_p,params->lpc_order_vt,1,params->lp_weighting,fundf,params->FS,index,glottsig,params->lp_stabilized);
			Filter(frame_concat,glottal,a_p);
		} else if(params->lp_method == LP_METHOD_ID_XLP) {
			SXLP(frame,a_p,params->lpc_order_vt,1);
			Filter(frame_concat,glottal,a_p);
		} else {
			WLPC(frame, a_p, params->lambda_vt);
			WFilter(frame_concat,glottal,a_p,B,params->lambda_vt);
		}

		/* De-emphasis */
		Integrator(frame, LIP_RADIATION);

		/* Convert LPC to LSF */
		Convert_matrix_to_LSF(LSF, a_p, index);
	}

	/* Define frame for unvoiced segments */
	gsl_vector *uvframe = gsl_vector_alloc(params->unvoiced_frame_length);
	int cntr = rint((frame->size-params->unvoiced_frame_length)/2.0);
	for(i=0;i<params->unvoiced_frame_length;i++)
		gsl_vector_set(uvframe,i,gsl_vector_get(frame,i+cntr));

	/* LPC for unvoiced segments, LPC/WLPC (possibility to use pre-emphasis, might make unvoiced
	 * excitation too high-pass and must be consistent with synthesis settings) */
	if(params->unvoiced_pre_emphasis == 1)
		Differentiate(uvframe,LIP_RADIATION);
	WLPC(uvframe, a_p, params->lambda_vt);
	if(params->unvoiced_pre_emphasis == 1)
		Integrator(uvframe, LIP_RADIATION);
	Convert_matrix_to_LSF(LSF2, a_p, index);

	/* Evaluate spectral tilt, LPC/WLPC */
	WLPC(glottal, a_tilt, params->lambda_gl);
	Convert_matrix_to_LSF(spectral_tilt, a_tilt, index);

	/* FFT of final SRC estimate and unvoiced frame */
	EvalFFTSpectrum(glottal, fftmatrix_src, index, params);
	EvalFFTSpectrum(uvframe, fftmatrix_uv, index, params);

	/* Free memory */
	gsl_vector_free(uvframe);
	gsl_vector_free(a_l);
	gsl_vector_free(a_p);
	gsl_vector_free(a_g);
	gsl_vector_free(a_g_iaif);
	gsl_vector_free(a_tilt);
	gsl_vector_free(frame_concat);
	gsl_vector_free(B);
}






/**
 * Function InverseFiltering_long
 *
 * Glottal inverse filtering for long frame:
 * Estimate source signal with Iterative Adaptive Inverse Filtering (IAIF)
 *
 * @param frame pointer to the samples
 * @param frame0 pointer to pre-frame samples
 * @param glottal pointer to the glottal source signal
 * @param lpc_order_vt LPC degree 1
 * @param LPC_g LPC degree 2
 * @param rho leaky integrator constant
 */
void InverseFiltering_long(gsl_vector *frame,
					gsl_vector *frame0,
					gsl_vector *glottal,
					gsl_vector *fundf,
					gsl_vector *glottsig,
					int index,
					PARAM *params) {

	// TODO:
	/* Parameters for this function will discard PARAMS below */
	int lp_method = 0;
	int lp_weighting = 0;

	/* Initialize */
	gsl_vector *a_l = gsl_vector_calloc(2);
	gsl_vector *a_p = gsl_vector_calloc(params->lpc_order_vt+1);
	gsl_vector *a_g = gsl_vector_calloc(params->lpc_order_gl+1);
	gsl_vector *a_g_iaif = gsl_vector_calloc(params->lpc_order_gl_iaif+1);
	gsl_vector *frame_concat = gsl_vector_alloc(params->lpc_order_vt+frame->size);
	gsl_vector *B = gsl_vector_alloc(1);
	gsl_vector_set(B,0,1);
	int i;

	/* Concatenate frame and frame0 */
	for(i=0; i<params->lpc_order_vt; i++)
		gsl_vector_set(frame_concat, i, gsl_vector_get(frame0, i));
	for(i=params->lpc_order_vt; i<frame->size+params->lpc_order_vt; i++)
		gsl_vector_set(frame_concat, i, gsl_vector_get(frame, i-params->lpc_order_vt));

	/* Use glottal inverse filtering */
	if(params->use_iaif == 1) {

		/* 2. WLPC-analysis (order 1) */
		/* 3. Inverse filtering */
		WLPC(frame, a_l, params->lambda_gl);
		WFilter(frame_concat, glottal, a_l, B, params->lambda_gl);
		Remove_mean(glottal);

		/* 4. Spectral analysis, select (W)LPC/(S)WLP/(S)XLP */
		/* 5. Inverse filtering, warped only if WLPC is used */
		if(lp_method == 1) {				//  SWLP-analysis (order p, M = p, lag = 1)
			SWLP(glottal,a_p,params->lpc_order_vt,1,lp_weighting,fundf,params->FS,index,glottsig,params->lp_stabilized);
			Filter(frame_concat, glottal, a_p);
		} else if(lp_method == 2) {			// SXLP-analysis (oder p, w = p)
			SXLP(glottal,a_p,params->lpc_order_vt,params->lp_stabilized);
			Filter(frame_concat, glottal, a_p);
		} else {							// WLPC-analysis (order p)
			WLPC(glottal, a_p, params->lambda_vt);
			WFilter(frame_concat, glottal, a_p,B,params->lambda_vt);
		}

		/* 6. Integration */
		Integrator(glottal, LIP_RADIATION);

		/* If full IAIF is used */
		if(params->use_mod_iaif == 0) {

			/* 7. WLPC-analysis (order g) */
			/* 8. Inverse filtering */
			/* 9. Integration */
			WLPC(glottal, a_g_iaif, params->lambda_gl);
			WFilter(frame_concat, glottal, a_g_iaif, B, params->lambda_gl);
			Integrator(glottal, LIP_RADIATION);

			/* 10. Spectral analysis, select (W)LPC/(S)WLP/(S)XLP,
			 * 11. Inverse filtering, warped only if WLPC is used */
			if(lp_method == 1) {				// SWLP-analysis (order p, M = p, lag = 1)
				SWLP(glottal,a_p,params->lpc_order_vt,1,lp_weighting,fundf,params->FS,index,glottsig,params->lp_stabilized);
				Filter(frame_concat, glottal, a_p);
			} else if(lp_method == 2) {			// SXLP-analysis (oder p, w = p)
				SXLP(glottal,a_p,params->lpc_order_vt,params->lp_stabilized);
				Filter(frame_concat, glottal, a_p);
			} else {							// LPC-analysis (order p), warping can be used
				WLPC(glottal, a_p, params->lambda_vt);
				WFilter(frame_concat, glottal, a_p, B, params->lambda_vt);
			}

			/* 12. Integration */
			Integrator(glottal, LIP_RADIATION);
		}

		/* If glottal inverse filtering is not used */
	} else {

		/* Pre-emphasis */
		Differentiate(frame,LIP_RADIATION);

		/* LPC */
		WLPC(frame, a_p, params->lambda_vt);
		WFilter(frame_concat, glottal, a_p, B, params->lambda_vt);
		Integrator(glottal, LIP_RADIATION);

		/* De-emphasis */
		Integrator(frame, LIP_RADIATION);
	}

	/* Free memory */
	gsl_vector_free(a_l);
	gsl_vector_free(a_p);
	gsl_vector_free(a_g);
	gsl_vector_free(a_g_iaif);
	gsl_vector_free(frame_concat);
	gsl_vector_free(B);
}


























/**
 * Function Define_current_f0
 *
 * Define current f0 value: If index is zero, use default value. If f0 is zero, use previous non-zero value
 *
 * @param fundf f0 vector
 * @param f0 previous f0 value
 * @param index
 */
double Define_current_f0(gsl_vector *fundf, double f0, int index) {
	if(index == 0)
		f0 = DEFAULT_F0;
	if(gsl_vector_get(fundf, index) != 0)
		f0 = gsl_vector_get(fundf, index);
	return f0;
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
 * Function Differentiate_matrix_plus_dpi_sqrt
 *
 * Differentiate matrix and add one value that is the distance of the last LSF to PI. Convert to SQRT.
 *
 * @param m1 matrix to be processed
 *
 */
void Differentiate_LSFs(gsl_matrix **m1) {

	/* Initialize, allocate new matrix (size+1) */
	int i,k;
	gsl_matrix *m2 = gsl_matrix_calloc((*m1)->size1,(*m1)->size2+1);

	/* Process: 1:orig 2-N:diff N+1:diff-pi */
	for(k=0;k<(*m1)->size1;k++) {

		/* Set first LSF as is */
		gsl_matrix_set(m2,k,0,gsl_matrix_get(*m1,k,0));

		/* Set the following differentials */
		for(i=1;i<(*m1)->size2;i++)
			gsl_matrix_set(m2,k,i,gsl_matrix_get((*m1),k,i)-gsl_matrix_get((*m1),k,i-1));

		/* Set the distance to PI */
		gsl_matrix_set(m2,k,m2->size2-1,M_PI-gsl_matrix_get((*m1),k,(*m1)->size2-1));
	}

	/* Take square root */
	for(k=0;k<m2->size1;k++)
		for(i=0;i<m2->size2;i++)
			gsl_matrix_set(m2,k,i,sqrt(gsl_matrix_get(m2,k,i)));

	/* Free old matrix and set pointer to new matrix */
	gsl_matrix_free(*(m1));
	*(m1) = m2;
}






















/**
 * Function Extract_pulses
 *
 * Extract pulses for pulse library
 *
 * @param ...
 *
 */
void Extract_pulses(gsl_vector *speech_frame, gsl_vector *frame_orig, gsl_vector *fundf, gsl_vector *naq, gsl_matrix *gpulses,
		gsl_matrix *gpulses_rs, gsl_vector *pulse_pos, gsl_vector *pulse_inds, gsl_vector *pulse_lengths, gsl_matrix *plsf,
		gsl_matrix *ptilt, gsl_matrix *pharm, gsl_matrix *phnr, gsl_vector *pgain, gsl_vector *ph1h2, gsl_vector *pnaq,
		gsl_matrix *lsf, gsl_matrix *spectral_tilt,	gsl_matrix *harmonics, gsl_vector *h1h2, gsl_matrix *hnr_i,gsl_vector *gain,
		gsl_matrix *pwaveform, gsl_matrix *waveform, int index, PARAM *params) {

	/* Stop if unvoiced */
	if(gsl_vector_get(fundf,index) == 0)
		return;

	/* Find glottal closure instants (GCIs). Stop if not found */
	gsl_vector *indices = Find_GCIs(frame_orig,fundf,params->FS,index);
	if(indices == NULL)
		return;

	/* Initialize */
	int i,j,original_pulse_ind,t0,n_added_pulses = 0,p = lsf->size2,pg = spectral_tilt->size2;
	int wav_start_ind = floor(params->waveform_samples/2.0);
	int wav_end_ind = wav_start_ind + params->waveform_samples;
	double temp;
	gsl_vector *pulse_rs = gsl_vector_alloc(params->rspulsemaxlen);
	gsl_vector *pulse_wf = gsl_vector_alloc(2*waveform->size2);
	gsl_matrix *frame_lsf = gsl_matrix_calloc(MAX_NUMBER_OF_PULSES_IN_FRAME,lsf->size2);
	gsl_matrix *frame_source_lsf = gsl_matrix_calloc(MAX_NUMBER_OF_PULSES_IN_FRAME,spectral_tilt->size2);
	gsl_vector *frame = gsl_vector_calloc(frame_orig->size);
	gsl_vector_memcpy(frame,frame_orig);
	gsl_vector *pulse_temp;
	gsl_vector *pulse_nowin;
	gsl_vector *t0frame,*t0frame_vt,*t0frame_concat,*lsf_tmp1,*lsf_tmp2,*lsf_g1,*lsf_g2,*a,*ag,*a1,*B;
	if(params->pitch_synchronous_analysis == 1) {
		lsf_tmp1 = gsl_vector_calloc(p);
		lsf_tmp2 = gsl_vector_calloc(p);
		lsf_g1 = gsl_vector_calloc(pg);
		lsf_g2 = gsl_vector_calloc(pg);
		ag = gsl_vector_calloc(pg+1);
		a = gsl_vector_calloc(p+1);
		a1 = gsl_vector_calloc(2);
		B = gsl_vector_alloc(1);
		gsl_vector_set(B,0,1);
	} else {
		lsf_tmp1 = NULL;
		lsf_tmp2 = NULL;
		lsf_g1 = NULL;
		lsf_g2 = NULL;
		ag = NULL;
		a1 = NULL;
		a = NULL;
		B = NULL;
	}

	/* Differentiate speech frame */
	Differentiate(frame,LEAK);

	/* Extract each complete two-period glottal flow waveform */
	i = 0;
	while(indices != NULL && i<indices->size) {
		original_pulse_ind = GSL_MAX(index*params->shift-ceil(params->f0_frame_length/2.0)+ceil(params->frame_length/2.0),0) + gsl_vector_get(indices,i);
		if(gsl_vector_get(indices,i) != 0) {
			if(i+2 > indices->size-1)
				break;
			if(gsl_vector_get(indices,i+2) > frame->size-2)
				break;
			pulse_temp = gsl_vector_alloc(gsl_vector_get(indices,i+2)-gsl_vector_get(indices,i)+1);
			pulse_nowin = gsl_vector_alloc(pulse_temp->size);

			/* Get samples and windowing */
			for(j=0;j<pulse_temp->size;j++) {
				gsl_vector_set(pulse_temp,j,gsl_vector_get(frame,gsl_vector_get(indices,i)+j)*HANN(j,pulse_temp->size));
				gsl_vector_set(pulse_nowin,j,gsl_vector_get(frame,gsl_vector_get(indices,i)+j));
			}

			/* If f0 is significantly different from the pulse period, the pulse may not be correctly estimated -> discard */
			if(fabs(pulse_temp->size-rint(2.0*params->FS/gsl_vector_get(fundf,index)))/rint(2.0*params->FS/gsl_vector_get(fundf,index)) > params->max_pulse_len_diff) {
				gsl_vector_free(pulse_temp);
				gsl_vector_free(pulse_nowin);
				i++;
				continue;
			}

			/* Give warning if pulses do not fit to the matrix */
			if(pulse_temp->size > params->pulsemaxlen)
				printf("ERROR: The pulse is longer (%i) than PULSEMAXLEN (%i).\n",(int)pulse_temp->size,(int)gpulses->size2);


			/******************************/
			/* Pitch-synchronous analysis */
			if(params->pitch_synchronous_analysis == 1) {

				/* Initialize pulse inverse filtering */
				t0 = pulse_temp->size;
				t0frame = gsl_vector_alloc(t0);
				t0frame_vt = gsl_vector_alloc(t0);
				t0frame_concat = gsl_vector_calloc(p+t0);

				/* Copy T0 segment from speech frame */
				for(j=0;j<t0;j++)
					gsl_vector_set(t0frame,j,gsl_vector_get(speech_frame,gsl_vector_get(indices,i)+j));

				/* Copy frame0 from speech frame */
				for(j=0; j<p; j++)
					if(gsl_vector_get(indices,i)+j-p >= 0)
						gsl_vector_set(t0frame_concat, j, gsl_vector_get(speech_frame, gsl_vector_get(indices,i)+j-p));
				for(j=p; j<t0+p; j++)
					gsl_vector_set(t0frame_concat, j, gsl_vector_get(t0frame, j-p));


				/**********************************/
				/* IAIF glottal inverse filtering */

				/* If IAIF is used */
				if(params->use_iaif == 1) {

					/* Evaluate glottal spectrum */
					Remove_mean(t0frame);
					WLPC(t0frame,a1,params->lambda_gl);

					/* Filter out glottal spectrum */
					WFilter(t0frame_concat, t0frame_vt, a1, B, params->lambda_gl);

					/* Evaluate VT spectrum */
					Remove_mean(t0frame_vt);
					WLPC(t0frame_vt,a,params->lambda_vt);

					/* Filter VT out */
					WFilter(t0frame_concat,t0frame_vt,a,B,params->lambda_vt);
					Remove_mean(t0frame_vt);

					/* If full IAIF is used */
					if(params->use_mod_iaif == 0) {

						/* Evaluate glottal spectrum with higher order */
						WLPC(t0frame_vt,ag,params->lambda_gl);

						/* Filter out glottal spectrum */
						WFilter(t0frame_concat, t0frame_vt, ag, B, params->lambda_gl);

						/* Evaluate VT spectrum */
						Remove_mean(t0frame_vt);
						WLPC(t0frame_vt,a,params->lambda_vt);

						/* Filter VT out */
						WFilter(t0frame_concat,t0frame_vt,a,B,params->lambda_vt);
						Remove_mean(t0frame_vt);
					}
				} else {

					/* Use basic LPC-based inverse filtering */
					WLPC(t0frame,a,params->lambda_vt);
					WFilter(t0frame_concat,t0frame_vt,a,B,params->lambda_vt);
					Remove_mean(t0frame_vt);
				}

				/* Evaluate glottal spectrum */
				WLPC(t0frame_vt,ag,params->lambda_gl);

				/* Convert to LSF, and copy LSFs to vector for averaging */
				Convert_vector_to_LSF(a, lsf_tmp1);
				Convert_vector_to_LSF(ag, lsf_g1);

				/* Set vocal tract and glottal spectrum */
				for(j=0;j<lsf_tmp2->size;j++)
					gsl_vector_set(lsf_tmp2,j,gsl_vector_get(lsf_tmp2,j) + gsl_vector_get(lsf_tmp1,j));
				for(j=0;j<lsf_g1->size;j++)
					gsl_vector_set(lsf_g2,j,gsl_vector_get(lsf_g2,j) + gsl_vector_get(lsf_g1,j));

				/* Set new pulse to pulse_temp (windowing) and to pulse_nowin */
				for(j=0;j<t0;j++) {
					gsl_vector_set(pulse_temp,j,gsl_vector_get(t0frame_vt,j)*HANN(j,t0));
					gsl_vector_set(pulse_nowin,j,gsl_vector_get(t0frame_vt,j));
				}

				/* Free memory */
				gsl_vector_free(t0frame_concat);
				gsl_vector_free(t0frame);
				gsl_vector_free(t0frame_vt);
			}

			/* Evaluate normalized amplitude quotient (NAQ) */
			double dpeak = BIG_POS_NUMBER;
			for(j=rint(pulse_nowin->size/4.0);j<rint(3.0*pulse_temp->size/4.0);j++)
				if(gsl_vector_get(pulse_nowin,j) < dpeak)
					dpeak = gsl_vector_get(pulse_nowin,j);
			dpeak = fabs(dpeak);
			Integrator(pulse_nowin,LEAK);
			double fac = BIG_NEG_NUMBER;
			for(j=rint(pulse_nowin->size/5.0);j<rint(4.0*pulse_temp->size/5.0);j++)
				if(gsl_vector_get(pulse_nowin,j) > fac)
					fac = gsl_vector_get(pulse_nowin,j);
			fac = fabs(fac);
			gsl_vector_set(naq,index,fac/dpeak/(pulse_nowin->size/2.0));

			/* Interpolate (rs) */
			Interpolate(pulse_temp,pulse_rs);

			/* Normalize the energy of the pulse */
			double e = 0;
			for(j=0;j<pulse_rs->size;j++)
				e += gsl_vector_get(pulse_rs,j)*gsl_vector_get(pulse_rs,j);
			e = sqrt(e);
			for(j=0;j<pulse_rs->size;j++)
				gsl_vector_set(pulse_rs,j,gsl_vector_get(pulse_rs,j)/e);

			/* Interpolate waveform to twice the number of samples in waveform,
			 * but take only half of the later in the center of the pulse */
			Interpolate(pulse_temp,pulse_wf);

			/* Normalize pulse waveform minimum to -1 */
			temp = 1;
			for(j=0;j<pulse_wf->size;j++)
				if(gsl_vector_get(pulse_wf,j) < temp)
					temp = gsl_vector_get(pulse_wf,j);
			temp = fabs(temp);
			if(temp >= 0)
				temp = 1;
			for(j=0;j<pulse_wf->size;j++)
				gsl_vector_set(pulse_wf,j,gsl_vector_get(pulse_wf,j)/temp);
			for(j=wav_start_ind;j<wav_end_ind;j++)
				gsl_matrix_set(waveform,index,j-wav_start_ind,gsl_matrix_get(waveform,index,j-wav_start_ind) + gsl_vector_get(pulse_wf,j)); // Sum pulses and divide by N in the end

			/* Set pulse parameters for pulse library */
			if(params->extract_pulselib_params == 1 && params->number_of_pulses < params->maxnumberofpulses) {

				/* Add gpulses */
				for(j=0;j<pulse_temp->size;j++)
					gsl_matrix_set(gpulses,params->number_of_pulses,j,gsl_vector_get(pulse_temp,j));

				/* Add gpulses_rs */
				for(j=0;j<pulse_rs->size;j++)
					gsl_matrix_set(gpulses_rs,params->number_of_pulses,j,gsl_vector_get(pulse_rs,j));

				/* Set length, position, and pulse index */
				gsl_vector_set(pulse_lengths,params->number_of_pulses,pulse_temp->size);
				gsl_vector_set(pulse_pos,params->number_of_pulses,index);
				gsl_vector_set(pulse_inds,params->number_of_pulses,original_pulse_ind);

				/* Set frame-wise pulse parameters */
				for(j=0;j<harmonics->size2;j++)
					gsl_matrix_set(pharm,params->number_of_pulses,j,gsl_matrix_get(harmonics,index,j));
				for(j=0;j<hnr_i->size2;j++)
					gsl_matrix_set(phnr,params->number_of_pulses,j,gsl_matrix_get(hnr_i,index,j));
				gsl_vector_set(pgain,params->number_of_pulses,gsl_vector_get(gain,index));
				gsl_vector_set(ph1h2,params->number_of_pulses,gsl_vector_get(h1h2,index));

				/* Set individual pulse parameters */
				for(j=wav_start_ind;j<wav_end_ind;j++)
					gsl_matrix_set(pwaveform,params->number_of_pulses,j-wav_start_ind,gsl_vector_get(pulse_wf,j));
				gsl_vector_set(pnaq,params->number_of_pulses,gsl_vector_get(naq,index));

				/* Set either individual or frame-wise parameters */
				if(params->pitch_synchronous_analysis == 1) {
					for(j=0;j<lsf->size2;j++)
						gsl_matrix_set(plsf,params->number_of_pulses,j,gsl_vector_get(lsf_tmp1,j));
					for(j=0;j<spectral_tilt->size2;j++)
						gsl_matrix_set(ptilt,params->number_of_pulses,j,gsl_vector_get(lsf_g1,j));
				} else {
					for(j=0;j<lsf->size2;j++)
						gsl_matrix_set(plsf,params->number_of_pulses,j,gsl_matrix_get(lsf,index,j));
					for(j=0;j<spectral_tilt->size2;j++)
						gsl_matrix_set(ptilt,params->number_of_pulses,j,gsl_matrix_get(spectral_tilt,index,j));
				}
			}

			/* Save vocal tract and source LSFs of each pulse (from which the closest to the average is selected).
			 * This is only used in pitch-synchronous analysis. */
			if(params->pitch_synchronous_analysis == 1) {
				for(j=0;j<lsf->size2;j++)
					gsl_matrix_set(frame_lsf,n_added_pulses,j,gsl_vector_get(lsf_tmp1,j));
				for(j=0;j<spectral_tilt->size2;j++)
					gsl_matrix_set(frame_source_lsf,n_added_pulses,j,gsl_vector_get(lsf_g1,j));
			}

			/* Increment pulse index */
			params->number_of_pulses++;
			n_added_pulses++;

			/* Give warning if pulses do not fit to the matrix (if pulses are stored) */
			if(params->extract_pulselib_params == 1 && params->number_of_pulses > params->maxnumberofpulses)
				printf("WARNING: The number of pulses has exceeded the value MAX_NUMBER_OF_PULSES (%i).\n",params->maxnumberofpulses);

			/* Free memory */
			gsl_vector_free(pulse_temp);
			gsl_vector_free(pulse_nowin);
		}
		i++;
	}

	/* Replace vocal tract and source LSFs from pulse library closest to the average of the frame */
	if(params->pitch_synchronous_analysis == 1 && n_added_pulses > 0) {
		gsl_vector *error = gsl_vector_calloc(n_added_pulses);
		for(i=0;i<n_added_pulses;i++)
			for(j=0;j<lsf->size2;j++)
				gsl_vector_set(error,i, gsl_vector_get(error,i) + fabs(gsl_vector_get(lsf_tmp2,j)/n_added_pulses - gsl_matrix_get(frame_lsf,i,j)));
		int minind = gsl_vector_min_index(error);
		for(j=0;j<lsf->size2;j++)
			gsl_matrix_set(lsf,index,j,gsl_matrix_get(frame_lsf,minind,j));
		for(j=0;j<spectral_tilt->size2;j++)
			gsl_matrix_set(spectral_tilt,index,j,gsl_matrix_get(frame_source_lsf,minind,j));
		gsl_vector_free(error);
	}

	/* Normalize the averaged pulse waveform so that the minimum is scaled to -1 */
	if(n_added_pulses > 0) {
		temp = 1;
		for(j=0;j<waveform->size2;j++) {
			if(gsl_matrix_get(waveform,index,j) < temp)
				temp = gsl_matrix_get(waveform,index,j);
		}
		temp = fabs(temp);
		for(j=0;j<waveform->size2;j++)
			gsl_matrix_set(waveform,index,j,gsl_matrix_get(waveform,index,j)/temp);
	}

	/* Free memory */
	gsl_matrix_free(frame_lsf);
	gsl_matrix_free(frame_source_lsf);
	gsl_vector_free(indices);
	gsl_vector_free(pulse_rs);
	gsl_vector_free(pulse_wf);
	gsl_vector_free(frame);
	if(params->pitch_synchronous_analysis == 1) {
		gsl_vector_free(lsf_tmp1);
		gsl_vector_free(lsf_tmp2);
		gsl_vector_free(lsf_g1);
		gsl_vector_free(lsf_g2);
		gsl_vector_free(a);
		gsl_vector_free(ag);
		gsl_vector_free(a1);
		gsl_vector_free(B);
	}
}









/**
 * Function Find_GCIs
 *
 * Find glottal closure instants (GCIs)
 *
 * @param ...
 *
 */
gsl_vector *Find_GCIs(gsl_vector *frame_orig, gsl_vector *fundf, int FS, int index) {

	int i,j,min_ind,ind_ind,t0_tmp;
	double t0,min_val,temp;
	gsl_vector *indices = gsl_vector_calloc(100);
	gsl_vector *final_inds;
	gsl_vector *frame = gsl_vector_calloc(frame_orig->size);

	/* Differentiate frame */
	gsl_vector_memcpy(frame,frame_orig);
	Differentiate(frame,LEAK);

	/* Find t0 minima of the glottal waveform in order to find GCIs,
	 * start evaluating from the mimina, first backward, then forward */
	min_ind = gsl_vector_min_index(frame);
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
		i++;
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
 * Function WFilter
 *
 * Warped FIR/IIR filter.
 *
 * This function is more or less adapted from WarpTB (MATLAB toolbox for frequency-warped signal processing)
 * authored by Aki H\E4rm\E4 and Matti Karjalainen.
 *
 * @param signal signal vector
 * @param A filter numerator
 * @param B filter denominator
 * @param lambda warping parameter
 *
 */
void WFilter(gsl_vector *signal, gsl_vector *result, gsl_vector *A, gsl_vector *B, double lambda) {

	int i,q,mlen;
    long int o;
    double xr,x,ffr,tmpr,Bb;
    double *sigma;
    long int len = signal->size;
    int adim = A->size;
    int bdim = B->size;
    double *Ar = (double *)calloc(adim,sizeof(double));
    double *Br = (double *)calloc(bdim,sizeof(double));
    double *ynr = (double *)calloc(signal->size,sizeof(double));
    double *rsignal = (double *)calloc(signal->size,sizeof(double));
    double *rmem = (double *)calloc(GSL_MAX(adim,bdim)+2,sizeof(double));

    /* Set signal to array */
    for(i=0;i<len;i++) {
		rsignal[i] = gsl_vector_get(signal,i);
	}

    /* Set A and B to arrays */
	for(i=0;i<adim;i++) {
		Ar[i] = gsl_vector_get(A,i);
	}
	for(i=0;i<bdim;i++) {
		Br[i] = gsl_vector_get(B,i);
	}

    /* Initialize */
    sigma = NFArray(bdim+2);
    alphas2sigmas(Br,sigma,lambda,bdim-1);
    if(adim >= bdim)
    	mlen = adim;
    else
    	mlen = bdim + 1;
    Bb = 1/Br[0];

    /* Warped filtering */
    for(o=0;o<len;o++) {

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

    /* Set output to result */
    int order = signal->size-result->size;
    for(i=order;i<signal->size;i++) {
    	gsl_vector_set(result,i-order,ynr[i]);
    }

    /* Free memory */
	free(ynr);
	free(rsignal);
	free(Ar);
	free(Br);
	free(rmem);
	free(sigma);
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
 * @return pointer to array
 *
 */
double *NFArray(int size) {

	double *p;
	p = (double *)calloc(sizeof(*p),size);
	return p;
}











/**
 * Function FundF
 *
 * Estimate fundamental frequency of a frame
 *
 * @param ...
 */
void FundF(gsl_vector *frame,
			gsl_vector *signal,
			gsl_matrix *fundf_candidates,
			gsl_vector *fundf,
			gsl_matrix *bp_gain,
			int index,
			PARAM *params) {

	/* Exit if external F0 file is used */
	if(params->use_external_f0 == 1)
		return;

	int i,n,zero_crossings,ind;
	double r_max;
	gsl_vector *max_inds = gsl_vector_calloc(NUMBER_OF_F0_CANDIDATES);
	gsl_vector *max_inds_interp = gsl_vector_calloc(NUMBER_OF_F0_CANDIDATES);
	gsl_vector *r = gsl_vector_calloc(frame->size);
	gsl_vector *r_copy = gsl_vector_calloc(frame->size);
	gsl_vector *r_norm = gsl_vector_alloc(frame->size);

	/* Count the number of zero crossings */
	zero_crossings = 0;
	for(i=0; i<signal->size-1; i++)
		if(GSL_SIGN(gsl_vector_get(signal, i)) != GSL_SIGN(gsl_vector_get(signal, i+1)))
			zero_crossings++;

	/* Autocorrelation sequence */
	for(i=0; i<frame->size; i++)
	    for(n=i; n<frame->size; n++)
	        gsl_vector_set(r, i, gsl_vector_get(r, i)+gsl_vector_get(frame, n)*gsl_vector_get(frame, n-i));

	/* Normalize r */
	r_max = gsl_vector_max(r);
	for(i=0; i<frame->size; i++)
	   gsl_vector_set(r_norm,i,gsl_vector_get(r,i)/r_max);

	/* Copy vector r for interpolation */
	gsl_vector_memcpy(r_copy,r);

	/* Clear samples when the index exceeds the fundamental frequency limits */
	for(i=0; i<rint(params->FS/params->fmax); i++)
		gsl_vector_set(r,i,NOTMAX);
	for(i=rint(params->FS/params->fmin); i<r->size; i++)
		gsl_vector_set(r,i,NOTMAX);

	/* Clear samples descending from the end-points */
	ind = rint(params->FS/params->fmax);
	while(gsl_vector_get(r,ind)-gsl_vector_get(r,ind+1) > 0) {
		gsl_vector_set(r,ind,NOTMAX);
		ind++;
		if(ind+1>r->size-1)
			break;
	}
	ind = rint(params->FS/params->fmin)-1;
	while(gsl_vector_get(r,ind)-gsl_vector_get(r,ind-1) > 0) {
		gsl_vector_set(r,ind,NOTMAX);
		ind--;
		if(ind-1<0) {
			break;
		}
	}

	/* Get T0 and F0 candidates */
	for(i=0;i<NUMBER_OF_F0_CANDIDATES;i++) {

		/* Get the i:th T0 index estimate */
		gsl_vector_set(max_inds,i,gsl_vector_max_index(r));

		/* Fit quadratic function to the peak of ACF to cancel the effect of the sampling period
		 * (i.e. parabolic interpolation) */
		gsl_vector_set(max_inds_interp,i,Parabolic_interpolation(r_copy,gsl_vector_get(max_inds,i),params));

		/* Set the F0 candidate */
		if(gsl_vector_get(max_inds_interp,i) <= 0) {
			gsl_matrix_set(fundf_candidates,index,i,0);
			break;
		} else {
			gsl_matrix_set(fundf_candidates,index,i,params->FS/gsl_vector_get(max_inds_interp,i));
				if(gsl_matrix_get(fundf_candidates,index,i) > params->fmax || gsl_matrix_get(fundf_candidates,index,i) < params->fmin || gsl_matrix_get(bp_gain, index, 0) < params->voicing_threshold/2.0 || zero_crossings > params->ZCR*2.0)
					gsl_matrix_set(fundf_candidates,index,i,0);
		}

		/* Clear the descending samples from the i:th maximum */
		ind = rint(gsl_vector_get(max_inds,i));
		while(gsl_vector_get(r,ind)-gsl_vector_get(r,ind+1) > 0) {
			gsl_vector_set(r,ind,NOTMAX);
			ind++;
			if(ind+1>r->size-1) {
				break;
			}
		}
		ind = GSL_MAX(rint(gsl_vector_get(max_inds,i)-1),1);
		while(gsl_vector_get(r,ind)-gsl_vector_get(r,ind-1) > 0) {
			gsl_vector_set(r,ind,NOTMAX);
			ind--;
			if(ind-1<0) {
				break;
			}
		}
	}

	/* Decide voiced/unvoiced. If voiced, set F0, otherwise set zero */
	if(gsl_matrix_get(bp_gain, index, 0) < params->voicing_threshold || zero_crossings > params->ZCR || gsl_vector_get(max_inds_interp,0) == 0)
		gsl_vector_set(fundf, index, 0);
	else
		gsl_vector_set(fundf, index, params->FS/gsl_vector_get(max_inds_interp,0));

    /* Free memory */
	gsl_vector_free(r);
	gsl_vector_free(r_copy);
	gsl_vector_free(max_inds);
	gsl_vector_free(max_inds_interp);
	gsl_vector_free(r_norm);
}









/**
 * Function Parabolic_interpolation
 *
 * Fit quadratic function to the peak of autocorrelation function (ACF) to cancel the effect of the sampling period
 * (i.e. parabolic interpolation)
 *
 * @param r ACF
 * @param maxind index at maximum value
 * @param params
 * @return T0 value
 */
double Parabolic_interpolation(gsl_vector *r, int maxind, PARAM *params) {

	/* Allocate variables */
	int i;
	double xi, chisq, T0;
	gsl_matrix *X, *cov;
	gsl_vector *y, *w, *c;
	gsl_multifit_linear_workspace *work;
	X = gsl_matrix_alloc(F0_INTERP_SAMPLES, 3);
	y = gsl_vector_alloc(F0_INTERP_SAMPLES);
	w = gsl_vector_alloc(F0_INTERP_SAMPLES);
	c = gsl_vector_alloc(3);
	cov = gsl_matrix_alloc(3, 3);
	work = gsl_multifit_linear_alloc(F0_INTERP_SAMPLES, 3);

	/* Set data */
	for(i=0;i<F0_INTERP_SAMPLES;i++) {
		xi = maxind-(F0_INTERP_SAMPLES-1)/2 + i;
		gsl_matrix_set(X, i, 0, 1.0);
		gsl_matrix_set(X, i, 1, xi);
		gsl_matrix_set(X, i, 2, xi*xi);
		gsl_vector_set(y, i, gsl_vector_get(r, GSL_MIN(GSL_MAX(maxind-(F0_INTERP_SAMPLES-1)/2 + i,0),r->size-1)));
		gsl_vector_set(w, i, 1.0);
	}

	/* Quadratic fitting */
	/* Evaluate the value of the function at the zero of the derivative.
	 * This is the interpolated length of the fundamental period in samples. */
	gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);
	T0 = -gsl_vector_get(c,1)/(2.0*gsl_vector_get(c,2));

	/* Free memory */
	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(w);
	gsl_vector_free(c);
	gsl_multifit_linear_free(work);

	/* Make sure the parabolic interpolation did perform succesfully */
	if(isnan(T0) == 0 && params->FS/T0 < params->fmax && params->FS/T0 > params->fmin)
		return T0;
	else
		return 0;
}












/**
 * Function F0_postprocess
 *
 * Various postprocessing for F0:
 *  -Median filtering
 *  -Filling small gaps
 *  -Tajectory estimation for discontinuities
 *
 * @param fundf F0
 * @param fundf_candidates candidate F0 values
 * @param params
 */
void F0_postprocess(gsl_vector *fundf, gsl_matrix *fundf_candidates, PARAM *params) {

	/* Exit if external F0 file is used */
	if(params->use_external_f0 == 1)
		return;

	/* Copy original F0 */
	gsl_vector *fundf_orig = gsl_vector_alloc(fundf->size);
	gsl_vector_memcpy(fundf_orig,fundf);

	/* Process */
	MedFilt3(fundf);
	Fill_f0_gaps(fundf, params);
	Fundf_postprocessing(fundf, fundf_orig, fundf_candidates, params);
	MedFilt3(fundf);
	Fill_f0_gaps(fundf, params);
	Fundf_postprocessing(fundf, fundf_orig, fundf_candidates, params);
	MedFilt3(fundf);

	/* Free memory */
	gsl_vector_free(fundf_orig);
}



/**
 * Fill_f0_gaps
 *
 * Fill small gaps in F0
 *
 * @param fundf F0
 */
void Fill_f0_gaps(gsl_vector *fundf, PARAM *params) {

	int i,j,voiced;
	double fundf_est,mu,std,sum,n,lim,ave;
	double f0jump01,f0jump45,f0jump02,f0jump03,f0jump12,f0jump13,f0jump42,f0jump43,f0jump52,f0jump53,f0jump05,f0jump14;
	gsl_vector *fundf_fill = gsl_vector_alloc(F0_FILL_RANGE);

	/* Estimate mean (mu) and standard deviation (std) of voiced parts */
	sum = 0;
	n = 0;
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) != 0) {
			sum += gsl_vector_get(fundf,i);
			n++;
		}
	}
	mu = sum/n;
	sum = 0;
	n = 0;
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) != 0) {
			sum += pow(mu-gsl_vector_get(fundf,i),2);
			n++;
		}
	}
	std = sqrt(sum/(n-1));

	/* Go throught all F0 values and fill small gaps (even if voiced) */
	for(i=0;i<fundf->size-F0_FILL_RANGE;i++) {
		fundf_est = 0;
		voiced = 0;
		for(j=0;j<F0_FILL_RANGE;j++) {
			gsl_vector_set(fundf_fill,j,gsl_vector_get(fundf,i+j));
			if(gsl_vector_get(fundf_fill,j) > 0) {
				voiced++;
				fundf_est += gsl_vector_get(fundf_fill,j);
			}
		}
		if(gsl_vector_get(fundf_fill,0) > 0 && gsl_vector_get(fundf_fill,1) > 0 && gsl_vector_get(fundf_fill,5) > 0 && gsl_vector_get(fundf_fill,4) > 0) {
			if(gsl_vector_get(fundf_fill,2) == 0 && gsl_vector_get(fundf_fill,3) == 0) {
				gsl_vector_set(fundf,i+2,fundf_est/4.0);
				gsl_vector_set(fundf,i+3,fundf_est/4.0);
			} else if(gsl_vector_get(fundf_fill,2) == 0) {
				gsl_vector_set(fundf,i+2,fundf_est/5.0);
			} else if(gsl_vector_get(fundf_fill,3) == 0) {
				gsl_vector_set(fundf,i+3,fundf_est/5.0);
			}
		}

		/* If all values are voiced, replace gaps of size two of which values differ significantly from average value */
		if(voiced == F0_FILL_RANGE) {
			f0jump01 = fabs(gsl_vector_get(fundf_fill,0)-gsl_vector_get(fundf_fill,1));
			f0jump45 = fabs(gsl_vector_get(fundf_fill,4)-gsl_vector_get(fundf_fill,5));
			f0jump02 = fabs(gsl_vector_get(fundf_fill,0)-gsl_vector_get(fundf_fill,2));
			f0jump03 = fabs(gsl_vector_get(fundf_fill,0)-gsl_vector_get(fundf_fill,3));
			f0jump12 = fabs(gsl_vector_get(fundf_fill,1)-gsl_vector_get(fundf_fill,2));
			f0jump13 = fabs(gsl_vector_get(fundf_fill,1)-gsl_vector_get(fundf_fill,3));
			f0jump42 = fabs(gsl_vector_get(fundf_fill,4)-gsl_vector_get(fundf_fill,2));
			f0jump43 = fabs(gsl_vector_get(fundf_fill,4)-gsl_vector_get(fundf_fill,3));
			f0jump52 = fabs(gsl_vector_get(fundf_fill,5)-gsl_vector_get(fundf_fill,2));
			f0jump53 = fabs(gsl_vector_get(fundf_fill,5)-gsl_vector_get(fundf_fill,3));
			f0jump05 = fabs(gsl_vector_get(fundf_fill,0)-gsl_vector_get(fundf_fill,5));
			f0jump14 = fabs(gsl_vector_get(fundf_fill,1)-gsl_vector_get(fundf_fill,4));
			lim = params->relative_f0_threshold*std;
			ave = (gsl_vector_get(fundf_fill,0) + gsl_vector_get(fundf_fill,1) + gsl_vector_get(fundf_fill,4) + gsl_vector_get(fundf_fill,5))/4.0;
			if(f0jump01 < lim && f0jump45 < lim && f0jump02 > lim && f0jump03 > lim && f0jump12 > lim && f0jump13 > lim && f0jump42 > lim && f0jump43 > lim && f0jump52 > lim && f0jump53 > lim && f0jump05 < lim && f0jump14 < lim) {
				gsl_vector_set(fundf,i+2,ave);
				gsl_vector_set(fundf,i+3,ave);
			}
		}
	}

	/* Go throught all F0 values and eliminate small voiced regions */
	for(i=0;i<fundf->size-F0_FILL_RANGE;i++) {
		for(j=0;j<F0_FILL_RANGE;j++) {
			gsl_vector_set(fundf_fill,j,gsl_vector_get(fundf,i+j));
		}
		if(gsl_vector_get(fundf_fill,0) == 0 && gsl_vector_get(fundf_fill,1) == 0 && gsl_vector_get(fundf_fill,5) == 0 && gsl_vector_get(fundf_fill,4) == 0) {
			if(gsl_vector_get(fundf_fill,2) > 0 || gsl_vector_get(fundf_fill,3) > 0) {
				gsl_vector_set(fundf,i+2,0);
				gsl_vector_set(fundf,i+3,0);
			}
		}
	}
	gsl_vector_free(fundf_fill);
}




/**
 * Fundf_postprocessing
 *
 * Refine f0-trajectory
 *
 * @param fundf fundamental frequency values
 * @param fundf_orig original fundamental frequency values
 * @param fundf_candidates candidates for F0
 * @param params parameter structure
 */
void Fundf_postprocessing(gsl_vector *fundf, gsl_vector *fundf_orig, gsl_matrix *fundf_candidates, PARAM *params) {

	int i,j,ind,voiced_ind_b,voiced_ind_f,unvoiced_ind,n;
	double f0_forward,f0_backward,x[params->f0_check_range],y[params->f0_check_range],w[params->f0_check_range],f0jump_b,f0jump_f;
	double c0,c1,cov00,cov01,cov11,chisq_f,chisq_b,mu,std,sum;

	/* Estimate mean (mu) and standard deviation (std) of voiced parts */
	sum = 0;
	n = 0;
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) != 0) {
			sum += gsl_vector_get(fundf,i);
			n++;
		}
	}
	mu = sum/n;
	sum = 0;
	n = 0;
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) != 0) {
			sum += pow(mu-gsl_vector_get(fundf,i),2);
			n++;
		}
	}
	std = sqrt(sum/(n-1));

	/* Go throught all F0 values and create weighted linear estimates for all F0 samples both from backward and forward */
	// TODO: Nonlinear fitting!
	if(params->f0_postprocessing == 1) {

		/* Start looping F0 */
		for(i=params->f0_check_range;i<fundf->size-params->f0_check_range;i++) {

			/* Check values backward (left to right) */
			ind = 0;
			voiced_ind_b = 0;
			unvoiced_ind = 0;
			for(j=params->f0_check_range;j>0;j--) {
				if(gsl_vector_get(fundf,i-j) == 0) {
					w[ind] = 0;
					unvoiced_ind++;
				} else {
					w[ind] = 1;
					voiced_ind_b++;
				}
				x[ind] = ind;
				y[ind] = gsl_vector_get(fundf,i-j);
				ind++;
			}
			if(unvoiced_ind == params->f0_check_range || voiced_ind_b < 3) {
				f0_backward = 0;
			} else {

				/* Weighted linear fitting, estimate new value */
				gsl_fit_wlinear(x,1,w,1,y,1,(double)params->f0_check_range,&c0,&c1,&cov00,&cov01,&cov11,&chisq_b);
				f0_backward = c0 + c1*(double)params->f0_check_range;
			}

			/* Check values forward (right to left) */
			ind = 0;
			voiced_ind_f = 0;
			unvoiced_ind = 0;
			for(j=1;j<params->f0_check_range+1;j++) {
				if(gsl_vector_get(fundf,i+j) == 0) {
					w[ind] = 0;
					unvoiced_ind++;
				} else {
					w[ind] = 1;
					voiced_ind_f++;
				}
				x[ind] = ind;
				y[ind] = gsl_vector_get(fundf,i+j);
				ind++;
			}
			if(unvoiced_ind == params->f0_check_range || voiced_ind_f < 3) {
				f0_forward = 0;
			} else {

				/* Weighted linear fitting, estimate new value */
				gsl_fit_wlinear(x,1,w,1,y,1,(double)params->f0_check_range,&c0,&c1,&cov00,&cov01,&cov11,&chisq_f);
				f0_forward = c0 - c1;
			}

			/* Evaluate relative jump in F0 from both directions */
			if(gsl_vector_get(fundf,i) != 0) {
				f0jump_b = fabs(gsl_vector_get(fundf,i)-gsl_vector_get(fundf,i-1))/gsl_vector_get(fundf,i);
				f0jump_f = fabs(gsl_vector_get(fundf,i)-gsl_vector_get(fundf,i+1))/gsl_vector_get(fundf,i);
			} else {
				f0jump_b = 0;
				f0jump_f = 0;
			}

			/* Set estimate if the relative jump is high enough,
			 * take into account the number of used voiced values
			 * in the fitting (voiced_ind) and the magnitude of the residual (chisq) */
			if(gsl_vector_get(fundf,i) != 0) {
				double new_f0 = gsl_vector_get(fundf,i);
				double f0jump = 0;
				if(voiced_ind_b > voiced_ind_f) { // Compare number of voiced frames
					if(f0_backward > params->fmin && f0_backward < params->fmax && f0jump_b >= params->relative_f0_threshold) {
						new_f0 = f0_backward;
						f0jump = f0jump_b;
					}
				} else if(voiced_ind_b < voiced_ind_f) {
					if(f0_forward > params->fmin && f0_forward < params->fmax && f0jump_f >= params->relative_f0_threshold) {
						new_f0 = f0_forward;
						f0jump = f0jump_f;
					}
				} else {
					if(chisq_b < chisq_f) { // Compare goodness of fit
						if(f0_backward > params->fmin && f0_backward < params->fmax && f0jump_b >= params->relative_f0_threshold) {
							new_f0 = f0_backward;
							f0jump = f0jump_b;
						}
					} else {
						if(f0_forward > params->fmin && f0_forward < params->fmax && f0jump_f >= params->relative_f0_threshold) {
							new_f0 = f0_forward;
							f0jump = f0jump_f;
						}
					}
				}

				/* Check if second F0 candidate is close. If so, set that value instead if estimated value */
				if(f0jump > params->relative_f0_threshold*std) {
					if(fabs(gsl_matrix_get(fundf_candidates,i,1) - new_f0) < f0jump)
						gsl_vector_set(fundf,i,gsl_matrix_get(fundf_candidates,i,1));
					else
						gsl_vector_set(fundf,i,new_f0);
				}
			}

			/* Make sure the correction does not lead the F0 trajectory to a wrong tract by estimating
			 * the cumulative difference between the original and new F0 curves.
			 * This is just a precaution, leading to a wrong F0 track should not normally happen */
			double cum_error = 0;
			for(j=GSL_MAX(i-F0_CUM_ERROR_RANGE,0);j<i;j++)
				cum_error += fabs(gsl_vector_get(fundf_orig,j)-gsl_vector_get(fundf,j))/mu;
			if(cum_error > F0_CUM_ERROR_LIM)
				gsl_vector_set(fundf,i,gsl_vector_get(fundf_orig,i));
		}
	}
}












/**
 * Function Integrator
 *
 * Leaky integrator
 *
 * @param frame pointer to the samples
 * @param rho leaky integrator constant
 */
void Integrator(gsl_vector *frame, double leak) {

	int i;
	double temp = 0;
	for(i=0;i<frame->size;i++) {
		gsl_vector_set(frame, i, gsl_vector_get(frame, i)+leak*temp);
		temp = gsl_vector_get(frame, i);
	}
}








/**
 * Function Add_missing_frames
 *
 * Add missing frames to the beginning and the end of parameters.
 * This is due to the analysis scheme, where the analysis starts
 * and ends with whole frames instead of zero-padding.
 *
 */
void Add_missing_frames(gsl_vector **pfundf,gsl_matrix **pfundf_candidates,gsl_vector **pgain,gsl_vector **puvgain,gsl_matrix **phnr,
		gsl_matrix **pspectral_tilt,gsl_matrix **pLSF,gsl_matrix **pLSF2,gsl_matrix **pharmonics,gsl_matrix **pwaveform,
		gsl_vector **ph1h2,gsl_vector **pnaq, gsl_matrix **pfftmatrix_vt, gsl_matrix **pfftmatrix_src, gsl_matrix **pfftmatrix_uv, PARAM *params) {

	int i,j;
	int empty_frames = rint((params->frame_length/(double)params->shift - 1)/2);
	int n_frames = (*(pfundf))->size;
	int n_frames_new = n_frames + 2*empty_frames;

	/* Set pointers */
	gsl_vector *fundf = *pfundf;
	gsl_vector *gain = *pgain;
	gsl_vector *uvgain = *puvgain;
	gsl_vector *h1h2 = *ph1h2;
	gsl_vector *naq = *pnaq;
	gsl_matrix *hnr = *phnr;
	gsl_matrix *spectral_tilt = *pspectral_tilt;
	gsl_matrix *LSF = *pLSF;
	gsl_matrix *LSF2 = *pLSF2;
	gsl_matrix *harmonics = *pharmonics;
	gsl_matrix *waveform = *pwaveform;
	gsl_matrix *fundf_candidates = *pfundf_candidates;
	gsl_matrix *fftmatrix_vt = *pfftmatrix_vt;
	gsl_matrix *fftmatrix_src = *pfftmatrix_src;
	gsl_matrix *fftmatrix_uv = *pfftmatrix_uv;

	/* Allocate new variables */
	gsl_vector *fundf_new = gsl_vector_alloc(n_frames_new);
	gsl_vector *gain_new = gsl_vector_alloc(n_frames_new);
	gsl_vector *uvgain_new = gsl_vector_alloc(n_frames_new);
	gsl_vector *h1h2_new = gsl_vector_alloc(n_frames_new);
	gsl_vector *naq_new = gsl_vector_alloc(n_frames_new);
	gsl_matrix *hnr_new = gsl_matrix_alloc(n_frames_new,hnr->size2);
	gsl_matrix *spectral_tilt_new = gsl_matrix_alloc(n_frames_new,spectral_tilt->size2);
	gsl_matrix *LSF_new = gsl_matrix_alloc(n_frames_new,LSF->size2);
	gsl_matrix *LSF2_new = gsl_matrix_alloc(n_frames_new,LSF2->size2);
	gsl_matrix *harmonics_new = gsl_matrix_alloc(n_frames_new,harmonics->size2);
	gsl_matrix *waveform_new = gsl_matrix_alloc(n_frames_new,waveform->size2);
	gsl_matrix *fundf_candidates_new = gsl_matrix_alloc(n_frames_new,fundf_candidates->size2);
	gsl_matrix *fftmatrix_vt_new = gsl_matrix_alloc(n_frames_new,fftmatrix_vt->size2);
	gsl_matrix *fftmatrix_src_new = gsl_matrix_alloc(n_frames_new,fftmatrix_src->size2);
	gsl_matrix *fftmatrix_uv_new = gsl_matrix_alloc(n_frames_new,fftmatrix_uv->size2);

	/* Fill in values */
	for(i=0;i<n_frames;i++) {
		gsl_vector_set(fundf_new,i+empty_frames,gsl_vector_get(fundf,i));
		gsl_vector_set(gain_new,i+empty_frames,gsl_vector_get(gain,i));
		gsl_vector_set(uvgain_new,i+empty_frames,gsl_vector_get(uvgain,i));
		gsl_vector_set(h1h2_new,i+empty_frames,gsl_vector_get(h1h2,i));
		gsl_vector_set(naq_new,i+empty_frames,gsl_vector_get(naq,i));
		for(j=0;j<hnr->size2;j++)
			gsl_matrix_set(hnr_new,i+empty_frames,j,gsl_matrix_get(hnr,i,j));
		for(j=0;j<spectral_tilt->size2;j++)
			gsl_matrix_set(spectral_tilt_new,i+empty_frames,j,gsl_matrix_get(spectral_tilt,i,j));
		for(j=0;j<LSF->size2;j++) {
			gsl_matrix_set(LSF_new,i+empty_frames,j,gsl_matrix_get(LSF,i,j));
			gsl_matrix_set(LSF2_new,i+empty_frames,j,gsl_matrix_get(LSF2,i,j));
		}
		for(j=0;j<harmonics->size2;j++)
			gsl_matrix_set(harmonics_new,i+empty_frames,j,gsl_matrix_get(harmonics,i,j));
		for(j=0;j<waveform->size2;j++)
			gsl_matrix_set(waveform_new,i+empty_frames,j,gsl_matrix_get(waveform,i,j));
		for(j=0;j<fundf_candidates->size2;j++)
			gsl_matrix_set(fundf_candidates_new,i+empty_frames,j,gsl_matrix_get(fundf_candidates,i,j));
		for(j=0;j<fftmatrix_vt->size2;j++) {
			gsl_matrix_set(fftmatrix_vt_new,i+empty_frames,j,gsl_matrix_get(fftmatrix_vt,i,j));
			gsl_matrix_set(fftmatrix_src_new,i+empty_frames,j,gsl_matrix_get(fftmatrix_src,i,j));
			gsl_matrix_set(fftmatrix_uv_new,i+empty_frames,j,gsl_matrix_get(fftmatrix_uv,i,j));
		}
	}

	/* Add empty values at the beginning */
	for(i=0;i<empty_frames;i++) {
		gsl_vector_set(fundf_new,i,gsl_vector_get(fundf,empty_frames));
		gsl_vector_set(gain_new,i,gsl_vector_get(gain,empty_frames));
		gsl_vector_set(uvgain_new,i,gsl_vector_get(uvgain,empty_frames));
		gsl_vector_set(h1h2_new,i,gsl_vector_get(h1h2,empty_frames));
		gsl_vector_set(naq_new,i,gsl_vector_get(naq,empty_frames));
		for(j=0;j<hnr->size2;j++)
			gsl_matrix_set(hnr_new,i,j,gsl_matrix_get(hnr,empty_frames,j));
		for(j=0;j<spectral_tilt->size2;j++)
			gsl_matrix_set(spectral_tilt_new,i,j,gsl_matrix_get(spectral_tilt,empty_frames,j));
		for(j=0;j<LSF->size2;j++) {
			gsl_matrix_set(LSF_new,i,j,gsl_matrix_get(LSF,empty_frames,j));
			gsl_matrix_set(LSF2_new,i,j,gsl_matrix_get(LSF2,empty_frames,j));
		}
		for(j=0;j<harmonics->size2;j++)
			gsl_matrix_set(harmonics_new,i,j,gsl_matrix_get(harmonics,empty_frames,j));
		for(j=0;j<waveform->size2;j++)
			gsl_matrix_set(waveform_new,i,j,gsl_matrix_get(waveform,empty_frames,j));
		for(j=0;j<fundf_candidates->size2;j++)
			gsl_matrix_set(fundf_candidates_new,i,j,gsl_matrix_get(fundf_candidates,empty_frames,j));
		for(j=0;j<fftmatrix_vt->size2;j++) {
			gsl_matrix_set(fftmatrix_vt_new,i,j,gsl_matrix_get(fftmatrix_vt,empty_frames,j));
			gsl_matrix_set(fftmatrix_src_new,i,j,gsl_matrix_get(fftmatrix_src,empty_frames,j));
			gsl_matrix_set(fftmatrix_uv_new,i,j,gsl_matrix_get(fftmatrix_uv,empty_frames,j));
		}
	}

	/* Add empty values in the end */
	for(i=n_frames_new-empty_frames-1;i<n_frames_new;i++) {
		gsl_vector_set(fundf_new,i,gsl_vector_get(fundf,n_frames-1));
		gsl_vector_set(gain_new,i,gsl_vector_get(gain,n_frames-1));
		gsl_vector_set(uvgain_new,i,gsl_vector_get(uvgain,n_frames-1));
		gsl_vector_set(h1h2_new,i,gsl_vector_get(h1h2,n_frames-1));
		gsl_vector_set(naq_new,i,gsl_vector_get(naq,n_frames-1));
		for(j=0;j<hnr->size2;j++)
			gsl_matrix_set(hnr_new,i,j,gsl_matrix_get(hnr,n_frames-1,j));
		for(j=0;j<spectral_tilt->size2;j++)
			gsl_matrix_set(spectral_tilt_new,i,j,gsl_matrix_get(spectral_tilt,n_frames-1,j));
		for(j=0;j<LSF->size2;j++) {
			gsl_matrix_set(LSF_new,i,j,gsl_matrix_get(LSF,n_frames-1,j));
			gsl_matrix_set(LSF2_new,i,j,gsl_matrix_get(LSF2,n_frames-1,j));
		}
		for(j=0;j<harmonics->size2;j++)
			gsl_matrix_set(harmonics_new,i,j,gsl_matrix_get(harmonics,n_frames-1,j));
		for(j=0;j<waveform->size2;j++)
			gsl_matrix_set(waveform_new,i,j,gsl_matrix_get(waveform,n_frames-1,j));
		for(j=0;j<fundf_candidates->size2;j++)
			gsl_matrix_set(fundf_candidates_new,i,j,gsl_matrix_get(fundf_candidates,n_frames-1,j));
		for(j=0;j<fftmatrix_vt->size2;j++) {
			gsl_matrix_set(fftmatrix_vt_new,i,j,gsl_matrix_get(fftmatrix_vt,n_frames-1,j));
			gsl_matrix_set(fftmatrix_src_new,i,j,gsl_matrix_get(fftmatrix_src,n_frames-1,j));
			gsl_matrix_set(fftmatrix_uv_new,i,j,gsl_matrix_get(fftmatrix_uv,n_frames-1,j));
		}
	}

	/* Free old variables */
	gsl_vector_free(fundf);
	gsl_vector_free(gain);
	gsl_vector_free(uvgain);
	gsl_vector_free(h1h2);
	gsl_vector_free(naq);
	gsl_matrix_free(hnr);
	gsl_matrix_free(spectral_tilt);
	gsl_matrix_free(LSF);
	gsl_matrix_free(LSF2);
	gsl_matrix_free(harmonics);
	gsl_matrix_free(waveform);
	gsl_matrix_free(fundf_candidates);
	gsl_matrix_free(fftmatrix_vt);
	gsl_matrix_free(fftmatrix_src);
	gsl_matrix_free(fftmatrix_uv);

	/* Set new variables to pointers */
	*(pfundf) = fundf_new;
	*(pgain) = gain_new;
	*(puvgain) = uvgain_new;
	*(ph1h2) = h1h2_new;
	*(pnaq) = naq_new;
	*(phnr) = hnr_new;
	*(pspectral_tilt) = spectral_tilt_new;
	*(pLSF) = LSF_new;
	*(pLSF2) = LSF2_new;
	*(pharmonics) = harmonics_new;
	*(pwaveform) = waveform_new;
	*(pfundf_candidates) = fundf_candidates_new;
	*(pfftmatrix_vt) = fftmatrix_vt_new;
	*(pfftmatrix_src) = fftmatrix_src_new;
	*(pfftmatrix_uv) = fftmatrix_uv_new;
}












/**
 * Function LPC
 *
 * Calculate Linear Prediction (LP) coefficients using
 * autocorrelation method.
 *
 * @param frame pointer to the samples
 * @param a pointer to LPC coefficient vector
 */
void LPC(gsl_vector *frame, gsl_vector *a) {

	int i,j,n,s,p = a->size-1;
	double win = 0,sum = 0;
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

	/* Autocorrelation sequence */
	for(i=0; i<p+1;i++) {
		for(n=i; n<wframe->size; n++) {
			gsl_vector_set(r, i, gsl_vector_get(r, i)+gsl_vector_get(wframe, n)*gsl_vector_get(wframe, n-i));
		}
	}

	/* Autocorrelation matrix (Toeplitz) */
	for(i=0; i<p;i++) {
	    for(j=0; j<p; j++) {
	        gsl_matrix_set(R, i, j, gsl_vector_get(r, abs(i-j)));
	        sum += gsl_vector_get(r, abs(i-j));
	    }
	}

	/* Vector b */
	for(i=1; i<p+1; i++) {
		gsl_vector_set(b, i-1, gsl_vector_get(r, i));
		sum += gsl_vector_get(r, i);
	}

	/* Ra=r solver (LU-decomposition) (Do not evaluate LU if sum = 0) */
	if(sum != 0) {
		gsl_linalg_LU_decomp(R, perm, &s);
		gsl_linalg_LU_solve(R, perm, b, a_temp);
	}

	/* Set LP-coefficients to vector "a" */
	for(i=1; i<a->size; i++) {
		gsl_vector_set(a, i, (-1)*gsl_vector_get(a_temp, i-1));
	}
	gsl_vector_set(a, 0, 1);

	/* Remove real roots */
	if(ROOT_SCALING == 1) {
		if(a->size > ROOT_SCALE_MIN_DEGREE) {
			RealRootScale(a);
		}
	}

	/* Replace NaN-values with zeros in case of all-zero frames */
	for(i=0;i<a->size;i++) {
		if(gsl_isnan(gsl_vector_get(a,i)))
			gsl_vector_set(a,i,0);
	}

	/* Free memory */
	gsl_vector_free(wframe);
	gsl_vector_free(a_temp);
	gsl_vector_free(r);
	gsl_vector_free(b);
	gsl_matrix_free(R);
	gsl_permutation_free(perm);
}













/**
 * Function LPC_Postfilter
 *
 * Enhance Formants by modifying the re-evaluated the LPC power spectrum,
 * and evaluating the LPC-coefficients again
 *
 * @param LSF pointer to the LSF matrix
 * @param params parameter structure
 */
void LPC_Postfilter(gsl_matrix *LSF, PARAM *params) {

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

	/* Loop for every index */
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
		ModPowerSpectrum(S,formants,nf,params->formant_enh_coeff,params->formant_enh_lpc_delta);

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
		Convert_matrix_to_LSF(LSF, a, fi);
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
void ModPowerSpectrum(gsl_vector *s, gsl_vector *formants, int n, double gamma, double delta) {

	int i,j;

	/* Modify spectrum at frequencies 0 and fs/2 within a constant number of bins from the formant peak */
	for(i=0;i<gsl_vector_get(formants,0) - delta;i++)
		gsl_vector_set(s,i,gsl_vector_get(s,i)*gamma);
	for(i=gsl_vector_get(formants,n-1) + delta+1;i<s->size/2+1;i++)
		gsl_vector_set(s,i,gsl_vector_get(s,i)*gamma);

	/* Modify spectrum within a constant number of bins from the formant peaks */
	for(i=0;i<n-1;i++) {
		for(j=gsl_vector_get(formants,i) + delta+1;j<gsl_vector_get(formants,i+1) - delta;j++)
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
 * Function WLPC
 *
 * Calculate Warped Linear Prediction (Warped LP) coefficients using
 * autocorrelation method.
 *
 * @param frame pointer to the samples
 * @param a pointer to coefficiets
 * @param p LPC degree
 * @param lambda warping coefficient
 */
void WLPC(gsl_vector *frame, gsl_vector *a, double lambda) {

	int i,j,s,p = a->size-1;
	double win = 0,sum = 0;
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
	        sum += gsl_vector_get(r, abs(i-j));
	    }
	}

	/* Vector b */
	for(i=1; i<p+1; i++) {
		gsl_vector_set(b, i-1, gsl_vector_get(r, i));
		sum += gsl_vector_get(r, i);
	}

	/* Ra=r solver (LU-decomposition) (Do not evaluate LU if sum = 0) */
	if(sum != 0) {
		gsl_linalg_LU_decomp(R, perm, &s);
		gsl_linalg_LU_solve(R, perm, b, a_temp);
	}

	/* Set LP-coefficients to vector "a" */
	for(i=1; i<a->size; i++) {
		gsl_vector_set(a, i, (-1)*gsl_vector_get(a_temp, i-1));
	}
	gsl_vector_set(a, 0, 1);

	/* Remove real roots */
	if(ROOT_SCALING == 1) {
		if(a->size > ROOT_SCALE_MIN_DEGREE) {
			RealRootScale(a);
		}
	}

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
	if(weighting == LP_WEIGHTING_ID_STE)
		Eval_STE_weight(wframe,weight,M,lag);
	else if(weighting == LP_WEIGHTING_ID_GCI) {
		//Eval_GCI_weight1(glottsig,weight,fundf,FS,index);  // Original
		Eval_GCI_weight2(glottsig,weight,fundf,FS,index);    // Manu's work
		//Eval_GCI_weight3(glottsig,weight,fundf,FS,index);  // Manu's work with F0 adaptation
	}

	/* Make sure weights are positive */
	for(i=0;i<weight->size;i++)
		if(gsl_vector_get(weight,i) <= 0)
			gsl_vector_set(weight,i,0.0001);

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

	/* Stabilize unstable filter by scaling the poles along the unit circle */
	if(stabilized == 0)
		Pole_stabilize(a);

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
 * Function Eval_GCI_weight1
 *
 * Find GCIs (Glottal Closure Instants) and construct weighting for de-emphasizing GCIs in (S)WLP
 *
 * @param ...
 *
 */
void Eval_GCI_weight1(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index) {

	/* Estimate glottal closure instants */
	gsl_vector *inds = Find_GCIs(glottsig, fundf, FS, index);

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
 * Function Eval_GCI_weight2
 *
 * Find GCIs (Glottal Closure Instants) and construct weighting for de-emphasizing GCIs in (S)WLP
 *
 * @param ...
 *
 */
void Eval_GCI_weight2(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index) {

	/* Estimate glottal closure instants */
	gsl_vector *inds = Find_GCIs(glottsig, fundf, FS, index);

	/* If unvoiced or GCIs not found, simply set weight to 1 */
	if(inds == NULL || gsl_vector_get(fundf,index) == 0) {
		if(inds != NULL)
			gsl_vector_free(inds);
		gsl_vector_set_all(weight,1);
		return;
	}

	/* Algorithm parameters */
	double pq = 0.05;
	double dq = 0.7;
	double d = 0.00001;
	int nramp = DEFAULT_NRAMP;

	/* Sanity check */
	if(dq + pq > 1.0)
	    dq = 1.0 - pq;

	/* Initialize */
	int i,j,t,t1 = 0,t2 = 0;

	/* Set weight according to GCIs */
	gsl_vector_set_all(weight,d);
	for(i=0;i<inds->size-1;i++) {
		t = gsl_vector_get(inds,i+1)-gsl_vector_get(inds,i);
		t1 = round(dq*t);
		t2 = round(pq*t);
		while(t1+t2 > t)
			t1 = t1-1;
		for(j=gsl_vector_get(inds,i)+t2;j<gsl_vector_get(inds,i)+t2+t1;j++)
			gsl_vector_set(weight,j,1);
		if(nramp > 0) {
			for(j=gsl_vector_get(inds,i)+t2;j<gsl_vector_get(inds,i)+t2+nramp;j++)
				gsl_vector_set(weight,j,(j-gsl_vector_get(inds,i)-t2+1)/(nramp+1));
			if(gsl_vector_get(inds,i)+t2+t1-nramp >= 0)
				for(j=gsl_vector_get(inds,i)+t2+t1-nramp;j<gsl_vector_get(inds,i)+t2+t1;j++)
					gsl_vector_set(weight,j,1.0-(j-gsl_vector_get(inds,i)-t2-t1+nramp+1)/(nramp+1));
		}
	}

	/* Set last weighting */
	i = i-1;
	int Nend = glottsig->size-(t2+gsl_vector_get(inds,i+1));
	if(t2+gsl_vector_get(inds,i+1) < glottsig->size) {
		if(t1+t2 < Nend) {
			for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+t1;j++)
				gsl_vector_set(weight,j,1);
			if(nramp > 0) {
				for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+nramp;j++)
					gsl_vector_set(weight,j,(j-gsl_vector_get(inds,i+1)-t2+1)/(nramp+1));
				for(j=gsl_vector_get(inds,i+1)+t2+t1-nramp;j<gsl_vector_get(inds,i+1)+t2+t1;j++)
					gsl_vector_set(weight,j,1.0-(j-gsl_vector_get(inds,i+1)-t2-t1+nramp+1)/(nramp+1));
			}
		} else {
			t1 = Nend-t2;
			for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+t1-1;j++)
				gsl_vector_set(weight,j,1);
			if(nramp > 0)
				for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+nramp;j++)
					gsl_vector_set(weight,j,(j-gsl_vector_get(inds,i+1)-t2+1)/(nramp+1));
		}
	}

	/* Free memory */
	gsl_vector_free(inds);
}








/**
 * Function Eval_GCI_weight3
 *
 * Find GCIs (Glottal Closure Instants) and construct weighting for de-emphasizing GCIs in (S)WLP
 *
 * @param ...
 *
 */
void Eval_GCI_weight3(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index) {

	/* Estimate glottal closure instants */
	gsl_vector *inds = Find_GCIs(glottsig, fundf, FS, index);

	/* If unvoiced or GCIs not found, simply set weight to 1 */
	if(inds == NULL || gsl_vector_get(fundf,index) == 0) {
		if(inds != NULL)
			gsl_vector_free(inds);
		gsl_vector_set_all(weight,1);
		return;
	}

	/* Algorithm parameters */
	double dq = 0.7;
	double pq = 0.05;
	double f0cut1 = 190;
	double f0cut2 = 320;
	double f0 = gsl_vector_get(fundf,index);
	if(f0 <= f0cut1) {
	    dq = 0.9;
	    pq = 0.05;
	} else if(f0 < f0cut2) {
	    dq = 0.55;
	    pq = 0.05;
	} else {
	    dq = 0.7;
	    pq = 0.05;
	}
	double d = 0.00001;
	int nramp = DEFAULT_NRAMP;

	/* Sanity check */
	if(dq + pq > 1.0)
	    dq = 1.0 - pq;

	/* Initialize */
	int i,j,t,t1 = 0,t2 = 0;

	/* Set weight according to GCIs */
	gsl_vector_set_all(weight,d);
	for(i=0;i<inds->size-1;i++) {
		t = gsl_vector_get(inds,i+1)-gsl_vector_get(inds,i);
		t1 = round(dq*t);
		t2 = round(pq*t);
		while(t1+t2 > t)
			t1 = t1-1;
		for(j=gsl_vector_get(inds,i)+t2;j<gsl_vector_get(inds,i)+t2+t1;j++)
			gsl_vector_set(weight,j,1);
		if(nramp > 0) {
			for(j=gsl_vector_get(inds,i)+t2;j<gsl_vector_get(inds,i)+t2+nramp;j++)
				gsl_vector_set(weight,j,(j-gsl_vector_get(inds,i)-t2+1)/(nramp+1));
			if(gsl_vector_get(inds,i)+t2+t1-nramp >= 0)
				for(j=gsl_vector_get(inds,i)+t2+t1-nramp;j<gsl_vector_get(inds,i)+t2+t1;j++)
					gsl_vector_set(weight,j,1.0-(j-gsl_vector_get(inds,i)-t2-t1+nramp+1)/(nramp+1));
		}
	}

	/* Set last weighting */
	i = i-1;
	int Nend = glottsig->size-(t2+gsl_vector_get(inds,i+1));
	if(t2+gsl_vector_get(inds,i+1) < glottsig->size) {
		if(t1+t2 < Nend) {
			for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+t1;j++)
				gsl_vector_set(weight,j,1);
			if(nramp > 0) {
				for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+nramp;j++)
					gsl_vector_set(weight,j,(j-gsl_vector_get(inds,i+1)-t2+1)/(nramp+1));
				for(j=gsl_vector_get(inds,i+1)+t2+t1-nramp;j<gsl_vector_get(inds,i+1)+t2+t1;j++)
					gsl_vector_set(weight,j,1.0-(j-gsl_vector_get(inds,i+1)-t2-t1+nramp+1)/(nramp+1));
			}
		} else {
			t1 = Nend-t2;
			for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+t1-1;j++)
				gsl_vector_set(weight,j,1);
			if(nramp > 0)
				for(j=gsl_vector_get(inds,i+1)+t2;j<gsl_vector_get(inds,i+1)+t2+nramp;j++)
					gsl_vector_set(weight,j,(j-gsl_vector_get(inds,i+1)-t2+1)/(nramp+1));
		}
	}

	// TEST
	/* Windowing */
	//for(i=0; i<weight->size; i++)
	//	gsl_vector_set(weight, i, gsl_vector_get(weight, i)*BLACKMAN(i,weight->size));

	/* Free memory */
	gsl_vector_free(inds);
}












/**
 * Function SXLP
 *
 * Calculate Stabilized eXtended Linear Prediction (SXLP) coefficients
 * (unstabilized can be evaluated by switching "stabilized" to zero)
 *
 * @param frame pointer to the samples
 * @param a pointer to coefficiets
 * @param p LPC degree
 * @param w weighting window length
 */
void SXLP(gsl_vector *frame, gsl_vector *a, int w, int stabilized) {

	int i,j,k,s,N = frame->size,p = a->size-1;
	double win = 0,sum = 0;
	gsl_vector *wframe = gsl_vector_alloc(frame->size);

	/* Windowing (choose window) */
	for(i=0; i<frame->size; i++) {
		if (WIN_TYPE == HANN_WIN) win = HANN(i,frame->size);
		else if (WIN_TYPE == BLACKMAN_WIN) win = BLACKMAN(i,frame->size);
		else if (WIN_TYPE == HAMMING_WIN) win = HAMMING(i,frame->size);
		else win = 1.0;
		gsl_vector_set(wframe, i, gsl_vector_get(frame, i)*win);
	}

	/* Initialize */
	double mem = (w-1.0)/w;
	gsl_vector *Z = gsl_vector_alloc(p+1);
	gsl_vector *Zp = gsl_vector_calloc(p+1);
	gsl_vector *ac = gsl_vector_calloc(p+1);
	gsl_matrix *D = gsl_matrix_calloc(N+p,p+1);
	for(i=0;i<N+p;i++)
		for(j=0;j<i+1;j++)
			if(i-j<p+1 && j<N)
				gsl_matrix_set(D,i,i-j,gsl_vector_get(wframe,j));

	/* Create matrix D */
	for(i=0;i<N+p;i++) {
		for(j=0;j<p+1;j++)
			gsl_vector_set(ac,j, mem*gsl_vector_get(ac,j) + (1.0-mem) * (fabs(gsl_matrix_get(D,i,0)) + fabs(gsl_matrix_get(D,i,j))));
		if(stabilized == 1) {
			for(j=1;j<p+1;j++)
				gsl_vector_set(Z,j,GSL_MAX(gsl_vector_get(ac,j),gsl_vector_get(Zp,j-1)));
			gsl_vector_set(Z,0,GSL_MAX(gsl_vector_get(ac,0),0));
		} else
			for(j=0;j<p+1;j++)
				gsl_vector_set(Z,j,gsl_vector_get(ac,j));
		for(j=0;j<p+1;j++)
			gsl_matrix_set(D,i,j,gsl_matrix_get(D,i,j)*gsl_vector_get(Z,j));
		if(stabilized == 1)
			for(j=0;j<p+1;j++)
				gsl_vector_set(Zp,j,gsl_vector_get(Z,j));
	}

	/* Autocorrelation matrix R (R = (DT*D)/N, size p*p) and vector b (size p) */
	gsl_matrix *R = gsl_matrix_calloc(p,p);
	gsl_vector *b = gsl_vector_calloc(p);
	for(i=0;i<D->size2;i++) {
		for(j=0;j<D->size2;j++) {
			if(i > 0 && j > 0) {
				for(k=0;k<D->size1;k++) {
					gsl_matrix_set(R,i-1,j-1,gsl_matrix_get(R,i-1,j-1) + gsl_matrix_get(D,k,i)*gsl_matrix_get(D,k,j));
				}
				gsl_matrix_set(R,i-1,j-1,gsl_matrix_get(R,i-1,j-1)/wframe->size);
			}
			if(i > 0 && j == 0) {
				for(k=0;k<D->size1;k++)
					gsl_vector_set(b,i-1,gsl_vector_get(b,i-1) + gsl_matrix_get(D,k,i)*gsl_matrix_get(D,k,j));
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

	/* Stabilize unstable filter by scaling the poles along the unit circle */
	if(stabilized == 0)
		Pole_stabilize(a);

	/* Free memory */
	gsl_matrix_free(D);
	gsl_vector_free(Z);
	gsl_vector_free(Zp);
	gsl_vector_free(ac);
	gsl_vector_free(wframe);
	gsl_matrix_free(R);
	gsl_vector_free(b);
	gsl_vector_free(a_temp);
	gsl_permutation_free(perm);
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

	/* Copy signal */
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
 * Function Filter
 *
 * FIR-Filter implementation.
 * NB! Filtering starts at i=p, so there is a pre-frame attached to the actual frame.
 *
 * @param frame pointer to the samples
 * @param result pointer to the filtered samples
 * @param a LPC coefficients
 */
void Filter(gsl_vector *frame, gsl_vector *result, gsl_vector *a) {

	double sum;
	int i,j,p = frame->size-result->size;

	/* Filter signal */
	for(i=p; i<frame->size; i++) {
		sum = 0;
		for(j=0; j<a->size; j++) {
			sum += gsl_vector_get(frame, i-j)*gsl_vector_get(a, j);
		}
		gsl_vector_set(result, i-p, sum);
	}
}





/**
 * Function RealRootScale
 *
 * Scale the real roots of the LPC-coefficients to near origo
 *
 * @param a pointer to the samples
 *
 */
void RealRootScale(gsl_vector *a) {

	int i,j,k;
	int n_c = a->size;
	int n_r = 2*(n_c-1);
    double coeffs[n_c];
    double roots[n_r];
    double real_roots[10] = {0};

    /* Copy coefficients to "coeffs" */
    gsl_vector_reverse(a);
    for(i=0; i<n_c; i++) {
    	coeffs[i] = gsl_vector_get(a, i);
    }

    /* Solve roots */
    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(n_c);
    gsl_poly_complex_solve(coeffs, n_c, w, roots);
    gsl_poly_complex_workspace_free(w);

    /* Find real roots */
    j = 0;
    for (i=0; i<n_r/2; i++) {
    	if(fabs(roots[2*i+1]) < ROOT_SCALE_ZERO_LIMIT) {
    		real_roots[j] = roots[2*i];
    		j++;
    	}
    }

    /* Scale roots */
    gsl_vector_reverse(a);
    gsl_vector *temp = gsl_vector_calloc(a->size);
    double y;
    for(i=0; i<10; i++) {
    	if(real_roots[i] != 0) {

    		/* Deconvolve */
    		y = 0;
    		for(j=0; j<n_c; j++) {
				gsl_vector_set(a, j, gsl_vector_get(a, j) + y*real_roots[i]);
				y = gsl_vector_get(a, j);
			}
			gsl_vector_set(a, a->size-1, 0);

			/* Scale the root */
			double scaled_factor[2] = {1, -1*ROOT_SCALE_FACTOR*real_roots[i]};

			/* Convolve with scaled root */
			for(j=0; j<a->size; j++) {
				y = 0;
				for(k=0; k<=GSL_MIN(j, 1); k++) {
					y += gsl_vector_get(a, j-k)*scaled_factor[k];
				}
				gsl_vector_set(temp, j, y);
			}

			/* Copy result to a */
			for(j=0; j<a->size; j++) {
				gsl_vector_set(a, j, gsl_vector_get(temp, j));
			}

    	} else {
    		break;
    	}
    }
    gsl_vector_free(temp);
}





/**
 * Function Convert_matrix_to_LSF
 *
 * Converts the LPC-coefficient matrix to Line Spectrum Frequencies (LSF)
 *
 * @param LSF pointer to the LSF matrix
 * @param a pointer to LPC-vector
 * @param index time index
 */
void Convert_matrix_to_LSF(gsl_matrix *LSF, gsl_vector *a, int index) {

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
			gsl_matrix_set(LSF, index, i, (i+1)*M_PI/(LSF->size2+1));
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
			gsl_vector_set(LSF, i, (i+1)*M_PI/(LSF->size+1));
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
 * Function MedFilt3_matrix
 *
 * 3-point median filtering for matrices.
 * Filter along the first dimension.
 *
 * @param frame pointer to matrix to be filtered
 *
 */
void MedFilt3_matrix(gsl_matrix *matrix) {

	int i,j;
	gsl_vector *temp = gsl_vector_alloc(matrix->size1);
	for(i=0;i<matrix->size2;i++) {
		for(j=0;j<matrix->size1;j++) {
			gsl_vector_set(temp,j,gsl_matrix_get(matrix,j,i));
		}
		MedFilt3(temp);
		for(j=0;j<matrix->size1;j++) {
			gsl_matrix_set(matrix,j,i,gsl_vector_get(temp,j));
		}
	}
	gsl_vector_free(temp);
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
 * Function Remove_mean
 *
 * Remove mean of given vector
 *
 * @param frame pointer to samples
 *
 */
void Remove_mean(gsl_vector *frame) {

	int i;
	double mean = 0;
	for(i=0; i<frame->size; i++) {
		mean += gsl_vector_get(frame, i);
	}
	mean = mean/(double)frame->size;
	for(i=0; i<frame->size; i++) {
		gsl_vector_set(frame, i, gsl_vector_get(frame, i)-mean);
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
 * Function lsf2poly
 *
 * Convert LSF to polynomial
 *
 * @param lsf_matrix
 * @param poly
 * @param index
 */
void lsf2poly(gsl_matrix *lsf_matrix, gsl_vector *poly, int index, int HMM) {

	int i,l = lsf_matrix->size2;
	gsl_vector *lsf_vector = gsl_vector_calloc(l);
	gsl_vector *fi_p = NULL, *fi_q = NULL;

	/* Get values to vector */
	for(i=0;i<l;i++)
		gsl_vector_set(lsf_vector,i,gsl_matrix_get(lsf_matrix,index,i));

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
	if((l-1)/2 > 0) {
		gsl_vector_free(fi_q);
	}
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
 * Function Convert_Hz2ERB
 *
 * Convert vector scale from Hz to ERB
 *
 * @param vector pointer to original HNR vector
 * @param vector_erb pointer to vector of ERB-scale HNR values
 *
 */
void Convert_Hz2ERB(gsl_vector *vector, gsl_vector *vector_erb, int FS) {

	int i,j,hnr_channels = vector_erb->size;
	gsl_vector *erb = gsl_vector_alloc(vector->size);
	gsl_vector *erb_sum = gsl_vector_calloc(hnr_channels);

	/* Evaluate ERB scale indices for vector */
	for(i=0;i<vector->size;i++)
		gsl_vector_set(erb,i,log10(0.00437*(i/(vector->size-1.0)*(FS/2.0))+1.0)/log10(0.00437*(FS/2.0)+1.0)*(hnr_channels-SMALL_NUMBER));

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
 * Function Harmonic_analysis
 *
 * Extract the magnitudes of N harmonics.
 * Extract the amount of aperiodicity according to smoothed upper and lower spectral envelopes.
 * Extract H1-H2
 *
 * @param frame pointer to voice source frame
 * @param harmonics pointer to harmonic matrix
 * @param hnr pointer to hnr matrix
 * @param h1h2 pointer to h1h2 vector
 * @param f0 fundamental frequency
 * @param FS sampling frequency
 * @param index time index
 *
 */
void Harmonic_analysis(gsl_vector *frame, gsl_matrix *harmonics, gsl_matrix *hnr, gsl_vector *h1h2, double f0, int FS, int index) {

	/* Variables */
	int i,j,guess_index = 0,h_values_size,harmonic_search_range;
	int hnr_channels = hnr->size2;
	double ind[MAX_HARMONICS] = {0};
	double guess_index_double = 0;
	double data[FFT_LENGTH_LONG] = {0};
	double h_values[MAX_HARMONICS] = {0};
	double n_values[MAX_HARMONICS] = {0};
	gsl_vector *fft = gsl_vector_calloc(FFT_LENGTH_LONG/2);
	gsl_vector *find;

	/* FFT (with windowing) */
	for (i=0; i<frame->size; i++)
		data[i] = gsl_vector_get(frame, i)*HANN(i,frame->size);
	gsl_fft_real_radix2_transform(data, 1, FFT_LENGTH_LONG);
	for(i=1; i<FFT_LENGTH_LONG/2; i++) {
		gsl_vector_set(fft, i, 20*log10(sqrt(pow(data[i], 2) + pow(data[FFT_LENGTH_LONG-i], 2))));
	}
	gsl_vector_set(fft, 0, 20*log10(fabs(data[0])));

	/* Set (possible) infinity values to zero */
	for(i=0;i<fft->size;i++) {
		if(!gsl_finite(gsl_vector_get(fft,i)))
			gsl_vector_set(fft,i,MIN_LOG_POWER);
	}

	/* Find the indices and magnitudes of the harmonic peaks */
	i = 0;
	while(1) {

		/* Define harmonics search range, decreasing to the higher frequencies */
		harmonic_search_range = GSL_MAX(HARMONIC_SEARCH_COEFF*f0/(FS/(double)FFT_LENGTH_LONG)*((fft->size-1-guess_index)/(double)(fft->size-1)),1.0);
		find = gsl_vector_alloc(harmonic_search_range);

		/* Estimate the index of the i_th harmonic
		 * Use an iterative estimation based on earlier values */
		if(i > 0) {
			guess_index_double = 0;
			for(j=0;j<i;j++)
				guess_index_double += ind[j]/(j+1.0)*(i+1.0);
			guess_index = (int)GSL_MAX(guess_index_double/j - (harmonic_search_range-1)/2.0,0);
		} else
			guess_index = (int)GSL_MAX(f0/(FS/(double)FFT_LENGTH_LONG) - (harmonic_search_range-1)/2.0,0);

		/* Stop search if the end (minus safe limit) of the fft vector or the maximum number of harmonics is reached */
		if(guess_index + rint(HNR_UNCERTAINTY_COEFF*FFT_LENGTH_LONG) > fft->size-1 || i > MAX_HARMONICS-1) {
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

	//if(index == 147) {
	//	printf("%lf\n",f0);
	//	VPrint1(fft);
	//	APrint1(h_values,h_values_size);
	//	APrint3(ind,h_values_size);
	//	APrint2(n_values,i);
	//	APrint4(ind_n,i);
	//	pause(1);
	//}

	/* Postfilter HNR and iterpolate vectors */
	gsl_vector *hnr_est = gsl_vector_alloc(h_values_size-1);
	gsl_vector *hnr_est_erb = gsl_vector_calloc(hnr_channels);
	for(i=0;i<hnr_est->size;i++)
		gsl_vector_set(hnr_est,i,n_values[i]-h_values[i]);
	MedFilt3(hnr_est);
	Convert_Hz2ERB(hnr_est, hnr_est_erb, FS);

	/* Set values to hnr matrix */
	for(i=0;i<hnr_channels;i++)
		gsl_matrix_set(hnr,index,i,gsl_vector_get(hnr_est_erb,i));

	/* Extract harmonics values, or actually their differences. First value is H1-H2. */
	double mag0 = gsl_vector_get(fft,ind[0]);
	for(i=0;i<harmonics->size2;i++)
		gsl_matrix_set(harmonics,index,i,mag0 - gsl_vector_get(fft,ind[i+1]));

	/* Set H1-H2 */
	gsl_vector_set(h1h2,index,gsl_matrix_get(harmonics,index,0));

	/* Free memory */
	gsl_vector_free(fft);
	gsl_vector_free(hnr_est);
	gsl_vector_free(hnr_est_erb);
}













/**
 * Function LSF_Postfilter
 *
 * Apply formant enhancement to LSFs
 *
 * @param lsf LSF-matrix
 * @param alpha postfilter coefficient alpha
 */
void LSF_Postfilter(gsl_matrix *lsf, PARAM *params) {

	if(params->formant_enh_coeff > 0) {
		int i,j;
		double d[lsf->size2-1];
		for(i=0;i<lsf->size1;i++) {
			for(j=0;j<lsf->size2-1;j++) {
				d[j] = params->formant_enh_coeff*(gsl_matrix_get(lsf,i,j+1) - gsl_matrix_get(lsf,i,j));
				if(j>0) {
					gsl_matrix_set(lsf,i,j, gsl_matrix_get(lsf,i,j-1) + d[j-1] + (pow(d[j-1],2)/(pow(d[j-1],2) + pow(d[j],2))) * ( gsl_matrix_get(lsf,i,j+1) - gsl_matrix_get(lsf,i,j-1) - d[j] - d[j-1] ) );
				}
			}
		}
	}
}









/**
 * Function Select_new_refined_values
 *
 * Select new parameter values for pulses according to smoothed and refined parameters.
 *
 * @param ...
 */
void Select_new_refined_values(gsl_vector *fundf, gsl_vector *gain, gsl_matrix *LSF, gsl_matrix *spectral_tilt,
		gsl_matrix *hnr_i, gsl_matrix *harmonics, gsl_vector *pgain, gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *phnr,
		gsl_matrix *pharm, gsl_vector *pulse_pos, gsl_vector *pulse_lengths, gsl_matrix *waveform,
		gsl_vector *h1h2, gsl_vector *ph1h2, gsl_vector *naq, gsl_vector *pnaq, PARAM *params) {

	/* Do not perform if parameters are not extracted */
	if(params->extract_pulselib_params == 0)
		return;

	/* Check that the number of pulses is not more than the maximum number of pulses */
	if(params->number_of_pulses > params->maxnumberofpulses)
		params->number_of_pulses = params->maxnumberofpulses;

	/* Select new parameter values for pulses according to smoothed and refined parameters
	 * Remeber to take into account that original parameter vectors are added with a few frames
	 * in the beginning and in the end (value empty_frames) */
	int i,j,empty_frames = rint((params->frame_length/(double)params->shift - 1)/2);
	for(i=0;i<params->number_of_pulses;i++) {
		gsl_vector_set(pgain,i,gsl_vector_get(gain,gsl_vector_get(pulse_pos,i)+empty_frames));
		gsl_vector_set(ph1h2,i,gsl_vector_get(h1h2,gsl_vector_get(pulse_pos,i)+empty_frames));
		gsl_vector_set(pnaq,i,gsl_vector_get(naq,gsl_vector_get(pulse_pos,i)+empty_frames));
		for(j=0;j<hnr_i->size2;j++)
			gsl_matrix_set(phnr,i,j,gsl_matrix_get(hnr_i,gsl_vector_get(pulse_pos,i)+empty_frames,j));
		for(j=0;j<harmonics->size2;j++)
			gsl_matrix_set(pharm,i,j,gsl_matrix_get(harmonics,gsl_vector_get(pulse_pos,i)+empty_frames,j));
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
	if(params->noise_reduction_analysis == 1) {
		int i;
		for(i=0;i<gain->size;i++)
			if(gsl_vector_get(gain,i) < params->noise_reduction_limit_db)
				gsl_vector_set(gain,i,gsl_vector_get(gain,i)-params->noise_reduction_db);
	}
}







/**
 * Function Fill_waveform_gaps
 *
 * Fill possible gaps in parameter waveform due to f0 postprocessing/in case waveform could not be analysed
 *
 * @param ...
 */
void Fill_waveform_gaps(gsl_vector *fundf, gsl_matrix *waveform) {

	int i,j,k,k1,k2;
	double sum;

	/* Copy original waveform */
	gsl_matrix *waveform_new = gsl_matrix_alloc(waveform->size1,waveform->size2);
	gsl_matrix_memcpy(waveform_new,waveform);

	/* Fix gaps in waveform */
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) > 0) {
			sum = 0;
			for(j=0;j<waveform->size2;j++)
				sum += fabs(gsl_matrix_get(waveform,i,j));
			if(sum == 0) {

				/* Gap found, find the next non-zero waveform */
				k1 = i+1;
				while(sum == 0) {
					if(k1 > waveform->size1-1)
						break;
					for(j=0;j<waveform->size2;j++)
						sum += fabs(gsl_matrix_get(waveform,k1,j));
					k1++;
				}
				sum = 0;
				k2 = i-1;
				while(sum == 0) {
					if(k2 < 0)
						break;
					for(j=0;j<waveform->size2;j++)
						sum += fabs(gsl_matrix_get(waveform,k2,j));
					k2--;
				}
				if(fabs(k1-i) < fabs(k2-i))
					k = k1-1;
				else
					k = k2+1;

				/* Fill the gap with the nearest found non-zero waveform */
				for(j=0;j<waveform->size2;j++)
					gsl_matrix_set(waveform_new,i,j,gsl_matrix_get(waveform,k,j));
			}
		} else {

			/* If f0 == 0, set zero */
			for(j=0;j<waveform->size2;j++)
				gsl_matrix_set(waveform,i,j,0);
		}
	}

	/* Copy fixed waveform matrix to the original one and free memory */
	gsl_matrix_memcpy(waveform,waveform_new);
	gsl_matrix_free(waveform_new);
}







/**
 * Function Fill_naq_gaps
 *
 * Fill possible gaps in parameter NAQ due to f0 postprocessing/in case waveform could not be analysed
 *
 * @param ...
 */
void Fill_naq_gaps(gsl_vector *fundf, gsl_vector *naq) {

	int i,k,k1,k2;

	/* Copy original vector */
	gsl_vector *naq_new = gsl_vector_alloc(naq->size);
	gsl_vector_memcpy(naq_new,naq);

	/* Fix gaps in naq */
	for(i=0;i<fundf->size;i++) {
		if(gsl_vector_get(fundf,i) > 0) {

			/* If gap is found, find the nearest non-zero value */
			if(gsl_vector_get(naq,i) == 0) {
				k1 = i+1;
				while(k1 < naq->size-1 && gsl_vector_get(naq,k1) == 0)
					k1++;
				k2 = i-1;
				while(k2 > 0 && gsl_vector_get(naq,k2) == 0)
					k2--;
				if(fabs(k1-i) < fabs(k2-i)) {
					if(k1 < naq->size-1)
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
				if(k < 0 || k > naq->size-1)
					break;

				/* Fill the gap with the nearest found non-zero value */
				gsl_vector_set(naq_new,i,gsl_vector_get(naq,k));

			}
		} else {

			/* If f0 == 0, set zero */
			gsl_vector_set(naq,i,0);
		}
	}

	/* Copy fixed vector to the original one and free memory */
	gsl_vector_memcpy(naq,naq_new);
	gsl_vector_free(naq_new);
}






/**
 * Function Construct_source
 *
 * Construct continuous source signal from the glottal inverse filtered signal
 *
 * @param ...
 */
void Construct_source(gsl_vector *source_signal, gsl_vector *glottal, gsl_vector *glottal_f0, int index, PARAM *params) {

	if(params->extract_source != 1)
		return;

	int i;
	int N = glottal->size;

	/* Get glottal samples to source signal */
	if(index == 0) {
		for(i=0; i<N/2; i++)
			gsl_vector_set(source_signal, i, gsl_vector_get(glottal,i));
		for(i=0; i<params->shift; i++)
			gsl_vector_set(source_signal, N/2+i, gsl_vector_get(glottal,N/2+i)*HANN(params->shift+i,2*params->shift));
	} else {
		for(i=0; i<2*params->shift; i++)
			gsl_vector_set(source_signal, (index-1)*params->shift+ceil(N/2)+i, gsl_vector_get(source_signal,(index-1)*params->shift+ceil(N/2)+i) + gsl_vector_get(glottal, ceil(N/2)-params->shift+i)*HANN(i,2*params->shift));
	}

	/* Get glottal samples to source signal from the longer frame (FO) */
	/*
	int Ng = glottal_f0->size;
	if(index == 0) {
		for(i=0; i<N/2; i++)
			gsl_vector_set(source_signal, i, gsl_vector_get(glottal_f0,i));
		for(i=0; i<shift; i++)
			gsl_vector_set(source_signal, N/2+i, gsl_vector_get(glottal_f0,N/2+i)*HANN(shift+i,2*shift));
	} else {
		for(i=0; i<2*shift; i++)
			gsl_vector_set(source_signal, (index-1)*shift+ceil(N/2)+i, gsl_vector_get(source_signal,(index-1)*shift+ceil(N/2)+i) + gsl_vector_get(glottal_f0, ceil(Ng/2)-shift+i)*HANN(i,2*shift));
	}
	*/
}




/**
 * Function Write_parameters_to_file
 *
 * Write speech parameters to file
 *
 * @param ...
 */
void Write_parameters_to_file(char *filename, gsl_matrix *LSF, gsl_matrix *LSF2, gsl_matrix *spectral_tilt,
		gsl_matrix *HNR, gsl_matrix* harmonics, gsl_matrix *waveform,gsl_vector *fundf, gsl_vector *gain,
		gsl_vector *h1h2, gsl_vector *naq, gsl_vector *source_signal, gsl_matrix *fftmatrix_vt,
		gsl_matrix *fftmatrix_src, PARAM *params) {

	/* Initialize */
	FILE *LSF_file = NULL, *LSF2_file = NULL, *gain_file = NULL, *fundf_file = NULL, *spectral_tilt_file = NULL, *params_file = NULL;
	FILE *HNR_file = NULL, *source_file = NULL, *harmonics_file = NULL, *waveform_file = NULL, *h1h2_file = NULL, *naq_file = NULL;
	FILE *FFTvt_file = NULL, *FFTsrc_file = NULL;
	int slen = strlen(filename)-4; // Remove ending ".wav"
	char orig[slen+1];
	char tmp[DEF_STRING_LEN];
	strncpy(orig, filename, slen);
	orig[slen] = '\0';

	/* Open files */
	if(params->extract_lsf == 1) {strcpy(tmp, orig);LSF_file = fopen(strcat(tmp, FILENAME_ENDING_LSF), "w"); if(LSF_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_gain == 1) {strcpy(tmp, orig);gain_file = fopen(strcat(tmp, FILENAME_ENDING_GAIN), "w"); if(gain_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_h1h2 == 1) {strcpy(tmp, orig);h1h2_file = fopen(strcat(tmp, FILENAME_ENDING_H1H2), "w"); if(h1h2_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_naq == 1) {strcpy(tmp, orig);naq_file = fopen(strcat(tmp, FILENAME_ENDING_NAQ), "w"); if(naq_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_f0 == 1) {strcpy(tmp, orig);fundf_file = fopen(strcat(tmp, FILENAME_ENDING_F0), "w"); if(fundf_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_tilt == 1) {strcpy(tmp, orig);spectral_tilt_file = fopen(strcat(tmp, FILENAME_ENDING_LSFSOURCE), "w"); if(spectral_tilt_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_info == 1) {strcpy(tmp, orig);params_file = fopen(strcat(tmp, FILE_ENDING_INFO), "w"); if(params_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_hnr == 1) {strcpy(tmp, orig);HNR_file = fopen(strcat(tmp, FILENAME_ENDING_HNR), "w"); if(HNR_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_harmonics == 1) {strcpy(tmp, orig);harmonics_file = fopen(strcat(tmp, FILENAME_ENDING_HARMONICS), "w");if(harmonics_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_waveform == 1) {strcpy(tmp, orig);waveform_file = fopen(strcat(tmp, FILENAME_ENDING_WAVEFORM), "w");if(waveform_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->extract_source == 1) {strcpy(tmp, orig);source_file = fopen(strcat(tmp, FILENAME_ENDING_SOURCESIGNAL), "w"); if(source_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->sep_vuv_spectrum == 1) {strcpy(tmp, orig);LSF2_file = fopen(strcat(tmp, FILENAME_ENDING_LSF2), "w"); if(LSF2_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->write_fftspectra_to_file == 1) {strcpy(tmp, orig);FFTvt_file = fopen(strcat(tmp, FILENAME_ENDING_FFT_VT), "w"); if(FFTvt_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}
	if(params->write_fftspectra_to_file == 1) {strcpy(tmp, orig);FFTsrc_file = fopen(strcat(tmp, FILENAME_ENDING_FFT_SRC), "w"); if(FFTsrc_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}}

	/* Write analysis parameters to file */
	if(params->extract_info == 1) {
		gsl_vector *parameters = gsl_vector_alloc(NPARAMS);
		gsl_vector_set(parameters, 0, params->frame_length_ms);
		gsl_vector_set(parameters, 1, params->shift_ms);
		gsl_vector_set(parameters, 2, (int)fundf->size); // This is inequal (but correct) to n_frames due to the "add_missing_frames" procedure!
		gsl_vector_set(parameters, 3, params->lpc_order_vt);
		gsl_vector_set(parameters, 4, params->lpc_order_gl);
		gsl_vector_set(parameters, 5, params->lambda_vt);
		gsl_vector_set(parameters, 6, params->lambda_gl);
		gsl_vector_set(parameters, 7, HNR->size2);
		gsl_vector_set(parameters, 8, harmonics->size2);
		gsl_vector_set(parameters, 9, params->number_of_pulses);
		gsl_vector_set(parameters, 10, params->pulsemaxlen_ms);
		gsl_vector_set(parameters, 11, params->rspulsemaxlen_ms);
		gsl_vector_set(parameters, 12, waveform->size2);
		gsl_vector_set(parameters, 13, params->FS);
		gsl_vector_set(parameters, 14, params->data_format);
		gsl_vector_fprintf(params_file, parameters, "%f");
		gsl_vector_free(parameters);
	}

	/* Write parameters to file */
	if(params->data_format == DATA_FORMAT_ID_ASCII) {
		if(params->extract_hnr == 1) gsl_matrix_fprintf(HNR_file, HNR, "%.5f");
		if(params->extract_lsf == 1) gsl_matrix_fprintf(LSF_file, LSF, "%.7f");
		if(params->extract_tilt == 1) gsl_matrix_fprintf(spectral_tilt_file, spectral_tilt, "%.7f");
		if(params->extract_gain == 1) gsl_vector_fprintf(gain_file, gain, "%.7f");
		if(params->extract_h1h2 == 1) gsl_vector_fprintf(h1h2_file, h1h2, "%.7f");
		if(params->extract_naq == 1) gsl_vector_fprintf(naq_file, naq, "%.7f");
		if(params->extract_f0 == 1) gsl_vector_fprintf(fundf_file, fundf, "%.7f");
		if(params->extract_harmonics == 1) gsl_matrix_fprintf(harmonics_file, harmonics, "%.7f");
		if(params->extract_waveform == 1) gsl_matrix_fprintf(waveform_file, waveform, "%.7f");
		if(params->extract_source == 1) gsl_vector_fprintf(source_file, source_signal, "%.7f");
		if(params->sep_vuv_spectrum == 1) gsl_matrix_fprintf(LSF2_file, LSF2, "%.7f");
		if(params->write_fftspectra_to_file == 1) gsl_matrix_fprintf(FFTvt_file, fftmatrix_vt, "%.7f");
		if(params->write_fftspectra_to_file == 1) gsl_matrix_fprintf(FFTsrc_file, fftmatrix_src, "%.7f");
	} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
		if(params->extract_hnr == 1) gsl_matrix_fwrite(HNR_file, HNR);
		if(params->extract_lsf == 1) gsl_matrix_fwrite(LSF_file, LSF);
		if(params->extract_tilt == 1) gsl_matrix_fwrite(spectral_tilt_file, spectral_tilt);
		if(params->extract_gain == 1) gsl_vector_fwrite(gain_file, gain);
		if(params->extract_h1h2 == 1) gsl_vector_fwrite(h1h2_file, h1h2);
		if(params->extract_naq == 1) gsl_vector_fwrite(naq_file, naq);
		if(params->extract_f0 == 1) gsl_vector_fwrite(fundf_file, fundf);
		if(params->extract_harmonics == 1) gsl_matrix_fwrite(harmonics_file, harmonics);
		if(params->extract_waveform == 1) gsl_matrix_fwrite(waveform_file, waveform);
		if(params->extract_source == 1) gsl_vector_fwrite(source_file, source_signal);
		if(params->sep_vuv_spectrum == 1) gsl_matrix_fwrite(LSF2_file, LSF2);
		if(params->write_fftspectra_to_file == 1) gsl_matrix_fwrite(FFTvt_file, fftmatrix_vt);
		if(params->write_fftspectra_to_file == 1) gsl_matrix_fwrite(FFTsrc_file, fftmatrix_src);
	}

	/* Close files */
	if(params->extract_lsf == 1) fclose(LSF_file);
	if(params->extract_gain == 1) fclose(gain_file);
	if(params->extract_h1h2 == 1) fclose(h1h2_file);
	if(params->extract_naq == 1) fclose(naq_file);
	if(params->extract_tilt == 1) fclose(spectral_tilt_file);
	if(params->extract_f0 == 1) fclose(fundf_file);
	if(params->extract_info == 1) fclose(params_file);
	if(params->extract_hnr == 1) fclose(HNR_file);
	if(params->extract_harmonics == 1) fclose(harmonics_file);
	if(params->extract_waveform == 1) fclose(waveform_file);
	if(params->extract_source == 1) fclose(source_file);
	if(params->sep_vuv_spectrum == 1) fclose(LSF2_file);
	if(params->write_fftspectra_to_file == 1) fclose(FFTvt_file);
	if(params->write_fftspectra_to_file == 1) fclose(FFTsrc_file);

	/* Save source signal to wav file */
	Save_source_signal_to_wavfile(filename, source_signal, params);

	/* Free memory */
	gsl_vector_free(source_signal);
	gsl_vector_free(gain);
	gsl_vector_free(h1h2);
	gsl_vector_free(naq);
	gsl_vector_free(fundf);
	gsl_matrix_free(spectral_tilt);
	gsl_matrix_free(LSF);
	gsl_matrix_free(LSF2);
	gsl_matrix_free(HNR);
	gsl_matrix_free(harmonics);
	gsl_matrix_free(waveform);
	gsl_matrix_free(fftmatrix_vt);
	gsl_matrix_free(fftmatrix_src);
}














/**
 * Function Select_unique_pulses
 *
 * Select only unique pulses from the pulses library, i.e. remove pulses that are extracted from different frames,
 * but are actually the same pulse instants. Among several same pulses, select the ones closest to global or local
 * average.
 *
 * @param ...
 */
void Select_unique_pulses(gsl_matrix **gpulses, gsl_matrix **gpulses_rs, gsl_vector **pulse_lengths, gsl_vector **pulse_pos, gsl_vector **pulse_inds,
		gsl_matrix **plsf, gsl_matrix **ptilt, gsl_matrix **phnr, gsl_matrix **pharm, gsl_matrix **pwaveform, gsl_vector **pgain,
		gsl_vector **ph1h2, gsl_vector **pnaq, PARAM *params) {

	/* Do not perform if there is no pulse library */
	if(params->extract_pulselib_params == 0 || params->number_of_pulses == 0)
		return;

	/* Initialize */
	int i,j;
	gsl_vector *best_pulse_inds = NULL;

	/* Do not perform if all pulses must be extracted */
	if(params->extract_only_unique_pulses == 1) {

		/* Initialize */
		int k,N = params->number_of_pulses;

		/* Evaluate global mean of the pulses */
		gsl_vector *global_mean = gsl_vector_calloc((*gpulses_rs)->size2);
		for(i=0;i<global_mean->size;i++)
			for(j=0;j<N;j++)
				gsl_vector_set(global_mean,i,gsl_vector_get(global_mean,i) + gsl_matrix_get((*gpulses_rs),j,i));
		for(i=0;i<global_mean->size;i++)
			gsl_vector_set(global_mean,i,gsl_vector_get(global_mean,i)/N);

		/* Find pulses with (approximately) the same starting sample index */
		gsl_vector *pulse_number = gsl_vector_alloc(N);
		gsl_vector *pulse_inds_tmp = gsl_vector_alloc((*pulse_inds)->size);
		gsl_vector_memcpy(pulse_inds_tmp,(*pulse_inds));
		double psit = PULSE_START_INDEX_THRESHOLD*(params->FS/16000.0);
		int new_flag;
		int ind = 0;
		for(i=0;i<N;i++) {
			new_flag = 0;
			if(gsl_vector_get(pulse_inds_tmp,i) > 0) {
				for(j=i;j<N;j++) {
					if(fabs(gsl_vector_get((*pulse_inds),i)-gsl_vector_get(pulse_inds_tmp,j)) < psit) {
						gsl_vector_set(pulse_number,j,ind);
						gsl_vector_set(pulse_inds_tmp,j,BIG_NEG_NUMBER);
						new_flag = 1;
					}
				}
				if(new_flag == 1)
					ind++;
			}
		}
		gsl_vector_free(pulse_inds_tmp);

		/* Select best estimate of the pulse from all the instances */
		int Nuniq = gsl_vector_max(pulse_number)+1;
		best_pulse_inds = gsl_vector_alloc(Nuniq);
		gsl_vector *local_mean = gsl_vector_alloc((*gpulses_rs)->size2);
		gsl_vector *pinds_tmp = gsl_vector_alloc(N);
		int npulses;
		for(i=0;i<Nuniq;i++) {

			/* Initialize */
			gsl_vector_set_zero(pinds_tmp);
			npulses = 0;

			/* Count and mark pulses belonging to the same time instance */
			for(j=0;j<N;j++) {
				if(gsl_vector_get(pulse_number,j) == i) {
					gsl_vector_set(pinds_tmp,npulses,j);
					npulses++;
				}
			}

			/* Evaluate best pulses:
			 *
			 * 1. Single instance, set directly */
			if(npulses == 1) {
				gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,0));

			/* 2. Two instances, select the one closest to global mean */
			} else if(npulses == 2) {
				double err1 = 0, err2 = 0;
				for(j=0;j<global_mean->size;j++) {
					err1 += fabs(gsl_vector_get(global_mean,j) - gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,0),j));
					err2 += fabs(gsl_vector_get(global_mean,j) - gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,1),j));
				}
				if(err1 < err2)
					gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,0));
				else
					gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,1));

			/* 3. More than two instances, select the one closest to local mean */
			} else {

				/* Eval local mean */
				gsl_vector_set_zero(local_mean);
				for(k=0;k<local_mean->size;k++)
					for(j=0;j<npulses;j++)
						gsl_vector_set(local_mean,k,gsl_vector_get(local_mean,k) + gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,j),k));
				for(k=0;k<local_mean->size;k++)
					gsl_vector_set(local_mean,k,gsl_vector_get(local_mean,k)/(double)npulses);

				/* Eval error */
				gsl_vector *error = gsl_vector_calloc(npulses);
				for(j=0;j<npulses;j++)
					for(k=0;k<local_mean->size;k++)
						gsl_vector_set(error,j,gsl_vector_get(error,j) + fabs(gsl_vector_get(local_mean,k) - gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,j),k)));

				/* Find minimum error */
				double min_error = BIG_POS_NUMBER;
				double min_index = 0;
				for(j=0;j<npulses;j++) {
					if(gsl_vector_get(error,j) < min_error) {
						min_error = gsl_vector_get(error,j);
						min_index = j;
					}
				}
				gsl_vector_free(error);

				/* Select pulse with minimum error */
				gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,min_index));
			}
		}
		gsl_vector_free(global_mean);
		gsl_vector_free(local_mean);
		gsl_vector_free(pinds_tmp);
		gsl_vector_free(pulse_number);

		/* Set number of pulses */
		params->number_of_pulses = Nuniq;
	}

	/* Resize pulse and parameter matrices */
	gsl_matrix *gpulses_new = gsl_matrix_alloc(params->number_of_pulses,(*gpulses)->size2);
	gsl_matrix *gpulses_rs_new = gsl_matrix_alloc(params->number_of_pulses,(*gpulses_rs)->size2);
	gsl_vector *pulse_pos_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *pulse_inds_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *pulse_lengths_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_matrix *plsf_new = gsl_matrix_alloc(params->number_of_pulses,(*plsf)->size2);
	gsl_matrix *ptilt_new = gsl_matrix_alloc(params->number_of_pulses,(*ptilt)->size2);
	gsl_matrix *pharm_new = gsl_matrix_alloc(params->number_of_pulses,(*pharm)->size2);
	gsl_matrix *phnr_new = gsl_matrix_alloc(params->number_of_pulses,(*phnr)->size2);
	gsl_matrix *pwaveform_new = gsl_matrix_alloc(params->number_of_pulses,(*pwaveform)->size2);
	gsl_vector *pgain_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *ph1h2_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *pnaq_new = gsl_vector_alloc(params->number_of_pulses);

	/* Set pulses and parameters */
	int index;
	for(i=0;i<params->number_of_pulses;i++) {

		/* Select index (depending on whether only unique or all pulses are saved  */
		if(params->extract_only_unique_pulses == 1)
			index = gsl_vector_get(best_pulse_inds,i);
		else
			index = i;

		/* Set pulses and parameters */
		for(j=0;j<(*gpulses)->size2;j++)
			gsl_matrix_set(gpulses_new,i,j,gsl_matrix_get((*gpulses),index,j));
		for(j=0;j<(*gpulses_rs)->size2;j++)
			gsl_matrix_set(gpulses_rs_new,i,j,gsl_matrix_get((*gpulses_rs),index,j));
		gsl_vector_set(pulse_pos_new,i,gsl_vector_get((*pulse_pos),index));
		gsl_vector_set(pulse_lengths_new,i,gsl_vector_get((*pulse_lengths),index));
		gsl_vector_set(pulse_inds_new,i,gsl_vector_get((*pulse_inds),index));
		for(j=0;j<(*plsf)->size2;j++)
			gsl_matrix_set(plsf_new,i,j,gsl_matrix_get((*plsf),index,j));
		for(j=0;j<(*ptilt)->size2;j++)
			gsl_matrix_set(ptilt_new,i,j,gsl_matrix_get((*ptilt),index,j));
		for(j=0;j<(*pharm)->size2;j++)
			gsl_matrix_set(pharm_new,i,j,gsl_matrix_get((*pharm),index,j));
		for(j=0;j<(*phnr)->size2;j++)
			gsl_matrix_set(phnr_new,i,j,gsl_matrix_get((*phnr),index,j));
		for(j=0;j<(*pwaveform)->size2;j++)
			gsl_matrix_set(pwaveform_new,i,j,gsl_matrix_get((*pwaveform),index,j));
		gsl_vector_set(pgain_new,i,gsl_vector_get((*pgain),index));
		gsl_vector_set(ph1h2_new,i,gsl_vector_get((*ph1h2),index));
		gsl_vector_set(pnaq_new,i,gsl_vector_get((*pnaq),index));
	}

	/* Free variables */
	if(params->extract_only_unique_pulses == 1)
		gsl_vector_free(best_pulse_inds);
	gsl_matrix_free(*gpulses);
	gsl_matrix_free(*gpulses_rs);
	gsl_vector_free(*pulse_inds);
	gsl_vector_free(*pulse_pos);
	gsl_vector_free(*pulse_lengths);
	gsl_matrix_free(*plsf);
	gsl_matrix_free(*ptilt);
	gsl_matrix_free(*pharm);
	gsl_matrix_free(*phnr);
	gsl_matrix_free(*pwaveform);
	gsl_vector_free(*pgain);
	gsl_vector_free(*ph1h2);
	gsl_vector_free(*pnaq);

	/* Set pointers */
	(*gpulses) = gpulses_new;
	(*gpulses_rs) = gpulses_rs_new;
	(*pulse_inds) = pulse_inds_new;
	(*pulse_pos) = pulse_pos_new;
	(*pulse_lengths) = pulse_lengths_new;
	(*plsf) = plsf_new;
	(*ptilt) = ptilt_new;
	(*pharm) = pharm_new;
	(*phnr) = phnr_new;
	(*pwaveform) = pwaveform_new;
	(*pgain) = pgain_new;
	(*ph1h2) = ph1h2_new;
	(*pnaq) = pnaq_new;
}





/**
 * Function Select_one_pulse_per_frame
 *
 * Select only one pulse per frame, i.e. remove pulses (except one) that are extracted from the same frame in order to remove bias
 * towards high-pitched pulses. Among several pulses in the frame, select the one closest to global or local average.
 *
 * @param ...
 */
void Select_one_pulse_per_frame(gsl_matrix **gpulses, gsl_matrix **gpulses_rs, gsl_vector **pulse_lengths, gsl_vector **pulse_pos, gsl_vector **pulse_inds,
		gsl_matrix **plsf, gsl_matrix **ptilt, gsl_matrix **phnr, gsl_matrix **pharm, gsl_matrix **pwaveform, gsl_vector **pgain,
		gsl_vector **ph1h2, gsl_vector **pnaq, PARAM *params) {

	/* Do not perform if there is no pulse library */
	if(params->extract_pulselib_params == 0 || params->number_of_pulses == 0)
		return;

	/* Initialize */
	int i,j;
	gsl_vector *best_pulse_inds = NULL;

	/* Do not perform if all pulses must be extracted */
	if(params->extract_one_pulse_per_frame == 1) {

		/* Initialize */
		int k,N = params->number_of_pulses;

		/* Evaluate global mean of the pulses */
		gsl_vector *global_mean = gsl_vector_calloc((*gpulses_rs)->size2);
		for(i=0;i<global_mean->size;i++)
			for(j=0;j<N;j++)
				gsl_vector_set(global_mean,i,gsl_vector_get(global_mean,i) + gsl_matrix_get((*gpulses_rs),j,i));
		for(i=0;i<global_mean->size;i++)
			gsl_vector_set(global_mean,i,gsl_vector_get(global_mean,i)/N);

		/* Find pulses from same frame */
		gsl_vector *pulse_number = gsl_vector_calloc(N);
		gsl_vector *pulse_pos_tmp = gsl_vector_alloc((*pulse_pos)->size);
		gsl_vector_memcpy(pulse_pos_tmp,(*pulse_pos));
		int new_flag;
		int ind = 0;
		for(i=0;i<N;i++) {
			new_flag = 0;
			for(j=i;j<N;j++) {
				if(fabs(gsl_vector_get((*pulse_pos),i) == gsl_vector_get(pulse_pos_tmp,j))) {
					gsl_vector_set(pulse_number,j,ind);
					gsl_vector_set(pulse_pos_tmp,j,BIG_NEG_NUMBER);
					new_flag = 1;
				}
			}
			if(new_flag == 1)
				ind++;
		}
		gsl_vector_free(pulse_pos_tmp);

		/* Select one pulse from each frame */
		int Nuniq = gsl_vector_max(pulse_number)+1;
		best_pulse_inds = gsl_vector_alloc(Nuniq);
		gsl_vector *local_mean = gsl_vector_alloc((*gpulses_rs)->size2);
		gsl_vector *pinds_tmp = gsl_vector_alloc(N);
		int npulses;
		for(i=0;i<Nuniq;i++) {

			/* Initialize */
			gsl_vector_set_zero(pinds_tmp);
			npulses = 0;

			/* Count and mark pulses belonging to the same frame */
			for(j=0;j<N;j++) {
				if(gsl_vector_get(pulse_number,j) == i) {
					gsl_vector_set(pinds_tmp,npulses,j);
					npulses++;
				}
			}

			/* Evaluate best pulses:
			 *
			 * 1. Single instance, set directly */
			if(npulses == 1) {
				gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,0));

			/* 2. Two instances, select the one closest to global mean */
			} else if(npulses == 2) {
				double err1 = 0, err2 = 0;
				for(j=0;j<global_mean->size;j++) {
					err1 += fabs(gsl_vector_get(global_mean,j) - gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,0),j));
					err2 += fabs(gsl_vector_get(global_mean,j) - gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,1),j));
				}
				if(err1 < err2)
					gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,0));
				else
					gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,1));

			/* 3. More than two instances, select the one closest to local mean */
			} else {

				/* Eval local mean */
				gsl_vector_set_zero(local_mean);
				for(k=0;k<local_mean->size;k++)
					for(j=0;j<npulses;j++)
						gsl_vector_set(local_mean,k,gsl_vector_get(local_mean,k) + gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,j),k));
				for(k=0;k<local_mean->size;k++)
					gsl_vector_set(local_mean,k,gsl_vector_get(local_mean,k)/(double)npulses);

				/* Eval error */
				gsl_vector *error = gsl_vector_calloc(npulses);
				for(j=0;j<npulses;j++)
					for(k=0;k<local_mean->size;k++)
						gsl_vector_set(error,j,gsl_vector_get(error,j) + fabs(gsl_vector_get(local_mean,k) - gsl_matrix_get((*gpulses_rs),gsl_vector_get(pinds_tmp,j),k)));

				/* Find minimum error */
				double min_error = BIG_POS_NUMBER;
				double min_index = 0;
				for(j=0;j<npulses;j++) {
					if(gsl_vector_get(error,j) < min_error) {
						min_error = gsl_vector_get(error,j);
						min_index = j;
					}
				}
				gsl_vector_free(error);

				/* Select pulse with minimum error */
				gsl_vector_set(best_pulse_inds,i,gsl_vector_get(pinds_tmp,min_index));
			}
		}
		gsl_vector_free(global_mean);
		gsl_vector_free(local_mean);
		gsl_vector_free(pinds_tmp);
		gsl_vector_free(pulse_number);

		/* Set number of pulses */
		params->number_of_pulses = Nuniq;
	}

	/* Resize pulse and parameter matrices */
	gsl_matrix *gpulses_new = gsl_matrix_alloc(params->number_of_pulses,(*gpulses)->size2);
	gsl_matrix *gpulses_rs_new = gsl_matrix_alloc(params->number_of_pulses,(*gpulses_rs)->size2);
	gsl_vector *pulse_pos_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *pulse_inds_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *pulse_lengths_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_matrix *plsf_new = gsl_matrix_alloc(params->number_of_pulses,(*plsf)->size2);
	gsl_matrix *ptilt_new = gsl_matrix_alloc(params->number_of_pulses,(*ptilt)->size2);
	gsl_matrix *pharm_new = gsl_matrix_alloc(params->number_of_pulses,(*pharm)->size2);
	gsl_matrix *phnr_new = gsl_matrix_alloc(params->number_of_pulses,(*phnr)->size2);
	gsl_matrix *pwaveform_new = gsl_matrix_alloc(params->number_of_pulses,(*pwaveform)->size2);
	gsl_vector *pgain_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *ph1h2_new = gsl_vector_alloc(params->number_of_pulses);
	gsl_vector *pnaq_new = gsl_vector_alloc(params->number_of_pulses);

	/* Set pulses and parameters */
	int index;
	for(i=0;i<params->number_of_pulses;i++) {

		/* Select index */
		if(params->extract_one_pulse_per_frame == 1)
			index = gsl_vector_get(best_pulse_inds,i);
		else
			index = i;

		/* Set pulses and parameters */
		for(j=0;j<(*gpulses)->size2;j++)
			gsl_matrix_set(gpulses_new,i,j,gsl_matrix_get((*gpulses),index,j));
		for(j=0;j<(*gpulses_rs)->size2;j++)
			gsl_matrix_set(gpulses_rs_new,i,j,gsl_matrix_get((*gpulses_rs),index,j));
		gsl_vector_set(pulse_pos_new,i,gsl_vector_get((*pulse_pos),index));
		gsl_vector_set(pulse_lengths_new,i,gsl_vector_get((*pulse_lengths),index));
		gsl_vector_set(pulse_inds_new,i,gsl_vector_get((*pulse_inds),index));
		for(j=0;j<(*plsf)->size2;j++)
			gsl_matrix_set(plsf_new,i,j,gsl_matrix_get((*plsf),index,j));
		for(j=0;j<(*ptilt)->size2;j++)
			gsl_matrix_set(ptilt_new,i,j,gsl_matrix_get((*ptilt),index,j));
		for(j=0;j<(*pharm)->size2;j++)
			gsl_matrix_set(pharm_new,i,j,gsl_matrix_get((*pharm),index,j));
		for(j=0;j<(*phnr)->size2;j++)
			gsl_matrix_set(phnr_new,i,j,gsl_matrix_get((*phnr),index,j));
		for(j=0;j<(*pwaveform)->size2;j++)
			gsl_matrix_set(pwaveform_new,i,j,gsl_matrix_get((*pwaveform),index,j));
		gsl_vector_set(pgain_new,i,gsl_vector_get((*pgain),index));
		gsl_vector_set(ph1h2_new,i,gsl_vector_get((*ph1h2),index));
		gsl_vector_set(pnaq_new,i,gsl_vector_get((*pnaq),index));
	}

	/* Free variables */
	if(params->extract_only_unique_pulses == 1)
		gsl_vector_free(best_pulse_inds);
	gsl_matrix_free(*gpulses);
	gsl_matrix_free(*gpulses_rs);
	gsl_vector_free(*pulse_inds);
	gsl_vector_free(*pulse_pos);
	gsl_vector_free(*pulse_lengths);
	gsl_matrix_free(*plsf);
	gsl_matrix_free(*ptilt);
	gsl_matrix_free(*pharm);
	gsl_matrix_free(*phnr);
	gsl_matrix_free(*pwaveform);
	gsl_vector_free(*pgain);
	gsl_vector_free(*ph1h2);
	gsl_vector_free(*pnaq);

	/* Set pointers */
	(*gpulses) = gpulses_new;
	(*gpulses_rs) = gpulses_rs_new;
	(*pulse_inds) = pulse_inds_new;
	(*pulse_pos) = pulse_pos_new;
	(*pulse_lengths) = pulse_lengths_new;
	(*plsf) = plsf_new;
	(*ptilt) = ptilt_new;
	(*pharm) = pharm_new;
	(*phnr) = phnr_new;
	(*pwaveform) = pwaveform_new;
	(*pgain) = pgain_new;
	(*ph1h2) = ph1h2_new;
	(*pnaq) = pnaq_new;
}





















/**
 * Function Write_pulselibrary_to_file
 *
 * Write pulse library parameters to file
 *
 * @param ...
 */
void Write_pulselibrary_to_file(char *filename, gsl_matrix *gpulses, gsl_matrix *gpulses_rs, gsl_vector *pulse_lengths, gsl_vector *pulse_pos, gsl_vector *pulse_inds,
		gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *phnr, gsl_matrix *pharm, gsl_matrix *pwaveform, gsl_vector *pgain, gsl_vector *ph1h2, gsl_vector *pnaq, PARAM *params) {

	/* Do not save data if not requested or pulses were not found */
	if(params->extract_pulselib_params == 1 && params->number_of_pulses > 0) {

		/* Initialize */
		int slen = strlen(filename)-4; // Remove ".wav" extension
		char orig[slen+1];
		char tmp[DEF_STRING_LEN];
		strncpy(orig, filename, slen);
		orig[slen] = '\0';

		/* Check that the number of pulses is not more than the maximum number of pulses */
		if(params->number_of_pulses > params->maxnumberofpulses)
			params->number_of_pulses = params->maxnumberofpulses;

		/* Save pulse data */
		FILE *gpulses_file,*gpulses_rs_file,*pulse_pos_file,*pulse_lengths_file,*pulse_inds_file;
		strcpy(tmp, orig);gpulses_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_PULSES), "w"); if(gpulses_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);gpulses_rs_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_RSPULSES), "w"); if(gpulses_rs_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pulse_pos_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_PULSEPOS), "w"); if(pulse_pos_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pulse_inds_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_PULSEINDS), "w"); if(pulse_inds_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pulse_lengths_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_PULSELENGTHS), "w"); if(pulse_lengths_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		if(params->data_format == DATA_FORMAT_ID_ASCII) {
			gsl_matrix_fprintf(gpulses_file, gpulses, "%.9f");
			gsl_matrix_fprintf(gpulses_rs_file, gpulses_rs, "%.9f");
			gsl_vector_fprintf(pulse_lengths_file, pulse_lengths, "%.0f");
			gsl_vector_fprintf(pulse_pos_file, pulse_pos, "%.0f");
			gsl_vector_fprintf(pulse_inds_file, pulse_inds, "%.0f");
		} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
			gsl_matrix_fwrite(gpulses_file, gpulses);
			gsl_matrix_fwrite(gpulses_rs_file, gpulses_rs);
			gsl_vector_fwrite(pulse_lengths_file, pulse_lengths);
			gsl_vector_fwrite(pulse_pos_file, pulse_pos);
			gsl_vector_fwrite(pulse_inds_file, pulse_inds);
		}
		fclose(gpulses_file);
		fclose(gpulses_rs_file);
		fclose(pulse_pos_file);
		fclose(pulse_inds_file);
		fclose(pulse_lengths_file);

		/* Save parameter data */
		FILE *plsf_file,*ptilt_file,*pharm_file,*phnr_file,*pgain_file,*pwaveform_file,*ph1h2_file,*pnaq_file;
		strcpy(tmp, orig);plsf_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_LSF), "w"); if(plsf_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);ptilt_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_TILT), "w"); if(ptilt_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pharm_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_HARM), "w"); if(pharm_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);phnr_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_HNR), "w"); if(phnr_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pwaveform_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_WAVEFORM), "w"); if(pwaveform_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pgain_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_GAIN), "w"); if(pgain_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);ph1h2_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_H1H2), "w"); if(ph1h2_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		strcpy(tmp, orig);pnaq_file = fopen(strcat(tmp, FILENAME_ENDING_PULSELIB_NAQ), "w"); if(pnaq_file==NULL) {printf("Error: Can't create file \"%s\".\n", tmp);}
		if(params->data_format == DATA_FORMAT_ID_ASCII) {
			gsl_matrix_fprintf(plsf_file, plsf, "%.7f");
			gsl_matrix_fprintf(ptilt_file, ptilt, "%.7f");
			gsl_matrix_fprintf(pharm_file, pharm, "%.7f");
			gsl_matrix_fprintf(phnr_file, phnr, "%.7f");
			gsl_matrix_fprintf(pwaveform_file, pwaveform, "%.7f");
			gsl_vector_fprintf(pgain_file, pgain, "%.7f");
			gsl_vector_fprintf(ph1h2_file, ph1h2, "%.7f");
			gsl_vector_fprintf(pnaq_file, pnaq, "%.7f");
		} else if(params->data_format == DATA_FORMAT_ID_BINARY) {
			gsl_matrix_fwrite(plsf_file, plsf);
			gsl_matrix_fwrite(ptilt_file, ptilt);
			gsl_matrix_fwrite(pharm_file, pharm);
			gsl_matrix_fwrite(phnr_file, phnr);
			gsl_matrix_fwrite(pwaveform_file, pwaveform);
			gsl_vector_fwrite(pgain_file, pgain);
			gsl_vector_fwrite(ph1h2_file, ph1h2);
			gsl_vector_fwrite(pnaq_file, pnaq);
		}
		fclose(plsf_file);
		fclose(ptilt_file);
		fclose(pharm_file);
		fclose(phnr_file);
		fclose(pwaveform_file);
		fclose(pgain_file);
		fclose(ph1h2_file);
		fclose(pnaq_file);
	}

	/* Free memory */
	gsl_matrix_free(gpulses);
	gsl_matrix_free(gpulses_rs);
	gsl_vector_free(pulse_inds);
	gsl_vector_free(pulse_pos);
	gsl_vector_free(pulse_lengths);
	gsl_matrix_free(plsf);
	gsl_matrix_free(ptilt);
	gsl_matrix_free(pharm);
	gsl_matrix_free(phnr);
	gsl_matrix_free(pwaveform);
	gsl_vector_free(pgain);
	gsl_vector_free(ph1h2);
	gsl_vector_free(pnaq);
}










/**
 * Function Free_variables
 *
 * Free analysis variables
 *
 * @param ...
 */
void Free_variables(gsl_vector *frame, gsl_vector *frame0, gsl_vector *signal, gsl_vector *glottal, gsl_vector *glottsig, gsl_vector *glottsig_f0,
		gsl_vector *uvgain, gsl_vector *f0_frame, gsl_vector *f0_frame0, gsl_vector *glottal_f0, gsl_matrix *bp_gain, gsl_matrix *fundf_candidates,
		gsl_matrix *fftmatrix_uv) {

	gsl_vector_free(frame);
	gsl_vector_free(frame0);
	gsl_vector_free(signal);
	gsl_vector_free(glottal);
	gsl_vector_free(glottsig);
	gsl_vector_free(glottsig_f0);
	gsl_vector_free(uvgain);
	gsl_vector_free(f0_frame);
	gsl_vector_free(f0_frame0);
	gsl_vector_free(glottal_f0);
	gsl_matrix_free(bp_gain);
	gsl_matrix_free(fundf_candidates);
	gsl_matrix_free(fftmatrix_uv);
}





/**
 * Function Save_source_signal_to_wavfile
 *
 * Save source signal to wav file
 *
 * @param filename
 * @param signal
 * @param params
 */
void Save_source_signal_to_wavfile(char *filename, gsl_vector *signal, PARAM *params) {

	if(params->extract_source == 0)
		return;

	/* Scale maximum to one if absmax is greater than one */
	int i;
	double absmax = GSL_MAX(gsl_vector_max(signal),-gsl_vector_min(signal));
	if(absmax > 1.0) {
		absmax = WAV_SCALE/absmax;
		for(i=0;i<signal->size;i++)
			gsl_vector_set(signal,i,gsl_vector_get(signal,i)*absmax);
	}

	/* Copy values to array */
	double *samples = (double *)calloc(signal->size, sizeof(double));
	for(i=0;i<signal->size;i++)
		samples[i] = gsl_vector_get(signal,i);

	/* Save synthesized speech to wav-file */
	SNDFILE *soundfile;
	SF_INFO sfinfo;
	sfinfo.samplerate = params->FS;
	sfinfo.channels = 1;
	sfinfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	char temp[DEF_STRING_LEN];
	strcpy(temp, filename);

	/* Open file with default ending or givend ending */
	int slen = strlen(filename)-4; // Remove ending ".wav"
	char tmp[DEF_STRING_LEN];
	strncpy(tmp, filename, slen);
	tmp[slen] = '\0';
	soundfile = sf_open(strcat(tmp,FILENAME_ENDING_SOURCESIGNAL_WAV), SFM_WRITE, &sfinfo);

	/* Check if success */
	if(soundfile==NULL) {
		printf("\n\nError creating file \"%s\": %s\n\n",temp,strerror(errno));
		return;
	}

	/* Write to file */
	sf_write_double(soundfile, samples, signal->size);

	/* Free memory */
	free(samples);
	sf_close(soundfile);
}







/**
 * Function Convert_F0_to_log
 *
 * Convert fundamental frequency vector to natural logarithm
 *
 * @param fundf
 * @param params
 */
void Convert_F0_to_log(gsl_vector *fundf, PARAM *params) {

	if(params->logf0 == 0)
		return;

	int i;
	for(i=0;i<fundf->size;i++)
		gsl_vector_set(fundf,i,log(gsl_vector_get(fundf,i)));
}






/**
 * Function EvalFFTSpectrum
 *
 * Calculate FFT of the signal and write spectrum to matrix
 *
 * @param signal pointer to signal vector
 * @param matrix matrix for the spectrum
 * @param index index of the parameters
 *
 */
void EvalFFTSpectrum(gsl_vector *signal, gsl_matrix *matrix, int index, PARAM *params) {

	/* Skip if not defined to be performed */
	if(params->write_fftspectra_to_file == 1) {

		/* Initialize */
		int i;
		double data[OUTPUT_FFT_LENGTH] = {0};

		/* FFT with windowing */
		for(i=0; i<signal->size; i++)
			data[i] = gsl_vector_get(signal, i)*HANN(i,signal->size);
		gsl_fft_real_radix2_transform(data, 1, OUTPUT_FFT_LENGTH);
		for(i=1; i<OUTPUT_FFT_LENGTH/2; i++) {
			gsl_matrix_set(matrix, index, i, sqrt(pow(data[i], 2) + pow(data[OUTPUT_FFT_LENGTH-i], 2)));
		}
		gsl_matrix_set(matrix, index, 0, data[0]);
	}
}



/**
 * Function WriteFFTSpectrumToFile
 *
 * Write spectrum matrix to file
 *
 * @param filename
 * @param fft matrix
 *
 */
void WriteFFTSpectrumToFile(char *filename, gsl_matrix *fft, PARAM *params) {

	/* Skip if not defined to be performed */
	if(params->write_fftspectra_to_file == 1) {

	}
}





































/*****************************************************************************/
/*                          FUNCTIONS NOT IN USE                             */
/*****************************************************************************/


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
    xi = x[0];
    while(i<length) {
    	gsl_vector_set(i_vector,i,gsl_spline_eval(spline, xi, acc));
    	xi += (len-1)/(double)(length-1);
    	i++;
    }

    /* Free memory */
    gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}







/**
 * Function Interpolate_poly
 *
 * Interpolate (polynomial) given vector to new vector of given length
 *
 * @param vector original vector
 * @param i_vector interpolated vector
 */
void Interpolate_poly(gsl_vector *vector, gsl_vector *i_vector) {

	int i,len = vector->size,length = i_vector->size;

	/* Read values to array */
	double x[len];
	double y[len];
	for(i=0; i<len; i++) {
		x[i] = i;
		y[i] = gsl_vector_get(vector,i);
	}

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_polynomial,len);
    gsl_spline_init(spline, x, y, len);
    double xi;
    i = 0;

    /* New implementation (27.3.2009, bug fix 8.2.2010) */
    xi = x[0];
    while(i<length) {
    	gsl_vector_set(i_vector,i,gsl_spline_eval(spline, xi, acc));
    	xi += (len-1)/(double)(length-1);
    	i++;
    }

    /* Free memory */
    gsl_spline_free(spline);
	gsl_interp_accel_free(acc);
}








/**
 * Function Pole_stabilize
 *
 * Check the validity of polynomial through mirroring poles outside the unit circle
 *
 * @param lpc polynomial vector
 */
void Pole_stabilize(gsl_vector *a) {

	int i,j,k;
	int n_c = a->size;
	int n_r = 2*(n_c-1);
	double r,yr,yc;
	double coeffs[n_c];
	double roots[n_r];
	gsl_vector *out_roots = gsl_vector_calloc(n_c);
	gsl_matrix *ac = gsl_matrix_calloc(2,a->size);
	gsl_matrix *temp = gsl_matrix_alloc(ac->size1,ac->size2);

	/* Copy coefficients to array coeffs and matrix ac */
	for(i=0; i<n_c; i++) {
		coeffs[i] = gsl_vector_get(a, a->size-i-1);
		gsl_matrix_set(ac,0,i,gsl_vector_get(a,i));
	}

	/* Solve roots */
	gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(n_c);
	gsl_poly_complex_solve(coeffs, n_c, w, roots);
	gsl_poly_complex_workspace_free(w);

	/* Find roots whose absolute value is greater than one */
	for(i=0; i<n_r/2; i++) {
		r = sqrt(powf(roots[2*i],2) + powf(roots[2*i+1],2));
		if(r > 1.0)
			gsl_vector_set(out_roots,i,1.0/powf(r,2));
	}

	/* Scale roots that are greater than one */
	for(i=0; i<n_c; i++) {
		if(gsl_vector_get(out_roots,i) > 0) {

			/* Deconvolve with root outside the unit circle (complex numbers) */
			yr = 0;yc = 0;
			for(j=0; j<n_c; j++) {
				gsl_matrix_set(ac, 0, j, gsl_matrix_get(ac, 0, j) + yr*roots[2*i] - yc*roots[2*i+1]);
				gsl_matrix_set(ac, 1, j, gsl_matrix_get(ac, 1, j) + yr*roots[2*i+1] + yc*roots[2*i]);
				yr = gsl_matrix_get(ac, 0, j);
				yc = gsl_matrix_get(ac, 1, j);
			}
			gsl_matrix_set(ac, 0, ac->size2-1, 0);
			gsl_matrix_set(ac, 1, ac->size2-1, 0);

			/* Scaled root */
			double scaled_root_r[2] = {1, -1*gsl_vector_get(out_roots,i)*roots[2*i]};
			double scaled_root_c[2] = {0, -1*gsl_vector_get(out_roots,i)*roots[2*i+1]};

			/* Convolve with scaled root (complex numbers) */
			gsl_matrix_set_all(temp,0);
			for(j=0; j<n_c; j++) {
				yr = 0;yc = 0;
				for(k=0; k<=GSL_MIN(j, 1); k++) {
					yr += gsl_matrix_get(ac, 0, j-k)*scaled_root_r[k] - gsl_matrix_get(ac, 1, j-k)*scaled_root_c[k];
					yc += gsl_matrix_get(ac, 0, j-k)*scaled_root_c[k] + gsl_matrix_get(ac, 1, j-k)*scaled_root_r[k];
				}
				gsl_matrix_set(temp, 0, j, yr);
				gsl_matrix_set(temp, 1, j, yc);
			}

			/* Copy result to ac */
			for(j=0; j<n_c; j++) {
				gsl_matrix_set(ac, 0, j, gsl_matrix_get(temp, 0, j));
				gsl_matrix_set(ac, 1, j, gsl_matrix_get(temp, 1, j));
			}
		}
	}

	/* Take only the real values and copy to original coefficient vector */
	for(j=0; j<n_c; j++)
		gsl_vector_set(a, j, gsl_matrix_get(ac, 0, j));

	/* Free memory */
	gsl_vector_free(out_roots);
	gsl_matrix_free(ac);
	gsl_matrix_free(temp);
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
 * Function IIR_filter
 *
 * IIR filter implementation (all-pole), coefficients in vector
 *
 * @param signal pointer to the samples
 * @param coeffs pointer to coefficients
 */
void IIR_filter(gsl_vector *signal, gsl_vector *coeffs) {

	int i,j,n = coeffs->size;
	double sum;

	/* Multiply filter coefficients by -1 */
	for(i=1;i<coeffs->size;i++) {
		gsl_vector_set(coeffs,i,-gsl_vector_get(coeffs,i));
	}

	/* Filter signal */
	for(i=0; i<signal->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, n-1); j++) {
			sum += gsl_vector_get(signal, i-j)*gsl_vector_get(coeffs,j);
		}
		gsl_vector_set(signal, i, sum);
	}
}



/**
 * Function FIR_filter_file
 *
 * FIR filter, coefficients from file
 *
 * @param signal pointer to the samples
 * @param filename name of the filter coefficient file
 * @param n length of the filter
 */
void FIR_filter_file(gsl_vector *signal, char *filename, int n) {

	int i,j;
	double sum;
	gsl_vector *temp = gsl_vector_alloc(signal->size);
	gsl_vector *coeffs = gsl_vector_alloc(n);
	FILE *COEFFS;

	/* Open file */
	COEFFS = fopen(filename, "r");
	if(COEFFS == NULL) {
		printf("\nError: Can't open file \"%s\".\n",filename);
	}

	/* Read coefficients */
	gsl_vector_fscanf(COEFFS,coeffs);
	fclose(COEFFS);

	/* Filter signal */
	for(i=0; i<signal->size; i++) {
		sum = 0;
		for(j=0; j<=GSL_MIN(i, n-1); j++) {
			sum += gsl_vector_get(signal, i-j)*gsl_vector_get(coeffs,j);
		}
		gsl_vector_set(temp, i, sum);
	}

	/* Copy "temp" samples to "signal" */
	for(i=0; i<signal->size; i++) {
		gsl_vector_set(signal, i, gsl_vector_get(temp, i));
	}

	/* Free memory */
	gsl_vector_free(temp);
	gsl_vector_free(coeffs);
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




// TEST
void Estimate_ERB_gain(gsl_vector *frame, gsl_matrix *erb_gain, int index, int FS) {

	/* Variables */
	int i,j,erb_channels = erb_gain->size2;
	double data[FFT_LENGTH_LONG] = {0};
	gsl_vector *fft = gsl_vector_calloc(FFT_LENGTH_LONG/2);
	gsl_vector *erb_inds = gsl_vector_alloc(FFT_LENGTH_LONG/2);
	gsl_vector *gain_erb = gsl_vector_calloc(erb_channels);
	gsl_vector *erb_index_sum = gsl_vector_calloc(erb_channels);

	/* FFT (with windowing) */
	for (i=0; i<frame->size; i++)
		data[i] = gsl_vector_get(frame, i)*HANN(i,frame->size);
	gsl_fft_real_radix2_transform(data, 1, FFT_LENGTH_LONG);
	for(i=1; i<FFT_LENGTH_LONG/2; i++) {
		gsl_vector_set(fft, i, 20*log10(sqrt(pow(data[i], 2) + pow(data[FFT_LENGTH_LONG-i], 2))));
	}
	gsl_vector_set(fft, 0, 20*log10(fabs(data[0])));

	/* Evaluate ERB scale indices */
	for(i=0;i<fft->size;i++)
		gsl_vector_set(erb_inds,i,log10(0.00437*(i/(fft->size-1.0)*(FS/2.0))+1.0)/log10(0.00437*(FS/2.0)+1.0)*(erb_channels-SMALL_NUMBER));

	/* Evaluate gain values according to ERB rate */
	for(i=0;i<fft->size;i++) {
		j = floor(gsl_vector_get(erb_inds,i));
		gsl_vector_set(gain_erb,j,gsl_vector_get(gain_erb,j)+gsl_vector_get(fft,i));
		gsl_vector_set(erb_index_sum,j,gsl_vector_get(erb_index_sum,j)+1);
	}

	/* Average values */
	for(i=0;i<erb_channels;i++)
		gsl_vector_set(gain_erb,i,gsl_vector_get(gain_erb,i)/gsl_vector_get(erb_index_sum,i));

	/* Set values to matrix */
	for(i=0;i<erb_channels;i++)
		gsl_matrix_set(erb_gain,index,i,gsl_vector_get(gain_erb,i));

	/* Free memory */
	gsl_vector_free(fft);
	gsl_vector_free(erb_inds);
	gsl_vector_free(erb_index_sum);
	gsl_vector_free(gain_erb);

}












/*********************************************************************/
/*                       TEST FUNCTIONS                              */
/*********************************************************************/

// TEST FUNCTION
// PRINT VECTOR TO FILE: p1.dat
void VPrint1(gsl_vector *vector) {
	FILE *f = fopen("p1.dat", "w");
	gsl_vector_fprintf(f, vector, "%.30f");
	fclose(f);
}

// TEST FUNCTION
// PRINT VECTOR TO FILE: p2.dat
void VPrint2(gsl_vector *vector) {
	FILE *f = fopen("p2.dat", "w");
	gsl_vector_fprintf(f, vector, "%.30f");
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
// PRINT MATRIX TO FILE: p3.dat
void MPrint1(gsl_matrix *matrix) {
	FILE *f = fopen("m1.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.30f");
	fclose(f);
}

// TEST FUNCTION
// PRINT MATRIX TO FILE: p4.dat
void MPrint2(gsl_matrix *matrix) {
	FILE *f = fopen("m2.dat", "w");
	gsl_matrix_fprintf(f, matrix, "%.30f");
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

