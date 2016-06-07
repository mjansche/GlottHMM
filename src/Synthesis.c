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
 * File Synthesis.c
 * Version: 1.1
 *
 */



/***********************************************/
/*                 INCLUDE                     */
/***********************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_vector.h>		/* GSL, Vector */
#include <gsl/gsl_matrix.h>		/* GSL, Matrix */
#include <gsl/gsl_fft_real.h>		/* GSL, FFT */
#include <libconfig.h>				/* Configuration file */
#include "SynthesisFunctions.h"





/*******************************************************************/
/*                          MAIN                                   */
/*******************************************************************/

int main(int argc, char *argv[]) {

	/* Start counter */
	double time1 = (double)clock();

	/* Check command line format */
	if(Check_command_line(argc) == EXIT_FAILURE)
		return EXIT_FAILURE;

	/* Read default configuration file and assign parameters */
	struct config_t *conf_def = Read_config(argv[2]);
	if(conf_def == NULL)
		return EXIT_FAILURE;
	PARAM params;
	if(Assign_config_parameters(argv[1],conf_def,&params,DEF_CONF) == EXIT_FAILURE) {
		config_destroy(conf_def);
		free(conf_def);
		return EXIT_FAILURE;
	}

	/* Read user configuration file and assign parameters */
	if(argc == 4) {
		struct config_t *conf_usr = Read_config(argv[3]);
		if(conf_usr == NULL)
			return EXIT_FAILURE;
		if(Assign_config_parameters(argv[1],conf_usr,&params,USR_CONF) == EXIT_FAILURE) {
			config_destroy(conf_usr);
			free(conf_usr);
			return EXIT_FAILURE;
		}
	}

	/* Print synthesis settings */
	Print_synthesis_settings_start(&params);

	/*******************************************************************/
	/*                 READ PULSE LIBRARY / PULSE                      */
	/*******************************************************************/


	/* Allocate space for pulse library and read data */
	gsl_matrix *pulses,*pulses_rs,*plsf,*ptilt,*pharm,*phnr,*pwaveform,*pca_pc,*pca_w_lib;
	gsl_vector *pgain,*pulse_lengths,*ph1h2,*pnaq,*pca_mean,*stoch_env,*stoch_sp;
	if(Read_pulse_library(&params,&pulses,&pulses_rs,&plsf,&ptilt,&pharm,&phnr,&pwaveform,&pca_pc,
			&pca_w_lib,&stoch_env,&stoch_sp,&pgain,&ph1h2,&pnaq,&pca_mean,&pulse_lengths) == EXIT_FAILURE)
		return EXIT_FAILURE;

	/* Read pulse file if pulse library is not used */
	gsl_vector *original_pulse = NULL;
	if(params.use_pulselib == 0) {
		original_pulse = ReadPulseFile(&params);
		if(original_pulse == NULL)
			return EXIT_FAILURE;
	}

	/* Read DNN pulse generation weights and allocate parameters */
	gsl_matrix **DNN_W = (gsl_matrix**)malloc(params.dnn_weight_dims->size/2*sizeof(gsl_matrix*));
	if(Read_DNN_weights(&params,DNN_W) == EXIT_FAILURE)
		return EXIT_FAILURE;
	gsl_vector **input_minmax = (gsl_vector**)malloc(sizeof(gsl_vector*));
	if(Read_input_minmax(&params,input_minmax) == EXIT_FAILURE)
			return EXIT_FAILURE;
	gsl_vector *dnnpulseindices = NULL;
	gsl_vector *dnnpulses = NULL;

	/* Define variables */
	gsl_vector *gain,*fundf,*excitation_voiced,*excitation_unvoiced;
	gsl_vector *gain_new,*pulse_clus_id,*h1h2,*naq,*resynthesis_pulse_index;
	gsl_matrix *LSF,*LSF2,*glflowsp,*glflowsp_new,*hnr,*hnr_new,*harmonics,*waveform,*pulse_clusters,*pca_w,*LSF_interp = NULL;
	int i;


	/****************************************************************************************/
	/*                          START SYNTHESING FILE(S)                                    */
	/****************************************************************************************/

	/* Start loop for synthesis list */
	for(i=0; i<params.synlistlen; i++) {

		/* Initialize params */
		if(Initialize_params(&params,i) == EXIT_FAILURE)
			return EXIT_FAILURE;

		/* Print synthesis settings */
		Print_synthesis_settings_middle(&params);

		/* Compatibility check */
		if(Compatibility_check(&params) == EXIT_FAILURE)
			return EXIT_FAILURE;


		/*****************************************************************************/
		/*                     LOAD AND MODIFY PARAMETERS                            */
		/*****************************************************************************/

		/* Allocate and read synthesis parameters, initialize pulse clustering if used */
		Allocate_params(&excitation_voiced,&excitation_unvoiced,&resynthesis_pulse_index,&gain_new,&glflowsp_new,&hnr_new,&params);
		if(Read_synthesis_parameters(&gain,&fundf,&LSF,&LSF2,&glflowsp,&hnr,&harmonics,&waveform,&h1h2,&naq,&pca_w,&params) == EXIT_FAILURE) continue;
		if(Pulse_clustering(&pulse_clus_id, &pulse_clusters, &params) == EXIT_FAILURE) continue;

		/* Miscellaneous operations */
		Convert_logF0_to_lin(fundf,&params);
		Merge_voiced_unvoiced_spectra(LSF,LSF2,fundf,&params);
		Integrate_LSFs(&LSF,&params);
		if(params.use_tilt == 1)
			Integrate_LSFs(&glflowsp,&params);
		LSF_fix_matrix(LSF);
		LSF_fix_matrix(glflowsp);
		Smooth_matrix(hnr,params.hnr_smooth_len);
		if(params.use_hmm == 0) Smooth_matrix(harmonics,params.harmonics_smooth_len);
		if(params.use_hmm == 0) MA(gain,params.gain_smooth_len);

		/* Modification for noise robust speech */
		Noise_robust_speech2(gain,harmonics,&params);

		/* Formant enhancement for LSFs */
		Postfilter(LSF,&params);

		/* Apply noise reduction */
		Noise_reduction(gain,&params);

		/* Normalize pulse library parameters according to synthesis parameters (normalize mean) */
		Normalize_pulse_library_var(plsf,ptilt,pharm,phnr,pwaveform,pgain,ph1h2,pnaq,LSF,glflowsp,
				harmonics,hnr,waveform,gain,h1h2,naq,fundf,&params);

		/* Adapt synthesis parameters according to pulse library parameters (normalize mean) */
		Adapt_synthesis_parameters_var(plsf,ptilt,pharm,phnr,pwaveform,pgain,ph1h2,pnaq,LSF,glflowsp,
				harmonics,hnr,waveform,gain,h1h2,naq,fundf,pulse_lengths,&params);

		/* Select vocal tract LSFs from pulse library */
		// TODO: Viterbi
		int winsize = 5;
		Select_LSFs_from_pulse_library(LSF,plsf,fundf,winsize,&params);

		/* Smooth and interpolate LSFs to signal length */
		LSF_interp = gsl_matrix_alloc(params.signal_length,params.lpc_order_vt);
		Smooth_interp_lsf(LSF_interp,LSF,params.signal_length,params.use_hmm,params.lsf_smooth_len);

		/* DNN pulse gen */
		if(params.use_dnn_pulsegen == 1) {
			dnnpulseindices = gsl_vector_alloc(params.signal_length);
			dnnpulses = gsl_vector_alloc(2*params.signal_length);
		}



		/*******************************************************************/
		/*                    CREATE EXCITATION                            */
		/*******************************************************************/


		printf("	- Creating excitation...\n");

		/* Create excitation for estimating the HNR of the voiced excitation */
		CreateExcitation(&params,excitation_voiced,excitation_unvoiced,fundf,gain,LSF,glflowsp,hnr,harmonics,waveform,h1h2,naq,
				original_pulse,pulses,pulses_rs,pulse_lengths,pgain,plsf,ptilt,phnr,pharm,pwaveform,ph1h2,pnaq,
				resynthesis_pulse_index,pulse_clus_id,pulse_clusters,gain,glflowsp_new,hnr_new,pca_mean,pca_pc,pca_w,pca_w_lib,
				stoch_env,stoch_sp,DNN_W,input_minmax,dnnpulseindices,dnnpulses);

		/* Compensate HNR etc */
		HNR_compensation(hnr,hnr_new,&params);
		Fill_pulse_indices(resynthesis_pulse_index);
		params.resynth = 1;

		/* Create excitation */
		CreateExcitation(&params,excitation_voiced,excitation_unvoiced,fundf,gain,LSF,glflowsp,hnr,harmonics,waveform,h1h2,naq,
				original_pulse,pulses,pulses_rs,pulse_lengths,pgain,plsf,ptilt,phnr,pharm,pwaveform,ph1h2,pnaq,
				resynthesis_pulse_index,pulse_clus_id,pulse_clusters,gain,glflowsp_new,hnr_new,pca_mean,pca_pc,pca_w,pca_w_lib,
				stoch_env,stoch_sp,DNN_W,input_minmax,dnnpulseindices,dnnpulses);
		Print_elapsed_time(&params);

		/* Compensate for the pre-emphasis during analysis (de-emphasis) */
		if(params.unvoiced_pre_emphasis == 1)
			Integrate(excitation_unvoiced, 0.97); // Slightly less than leaky factor


		/*******************************************************************/
		/*                    SPECTRAL MATCHING                            */
		/*******************************************************************/


		if((params.use_pulselib == 0 && params.use_dnn_pulsegen == 0) || (params.use_pulselib_pca == 1 && params.pca_spectral_matching == 1) ||
			 (params.use_dnn_pulsegen == 1 && params.use_dnn_specmatch == 1)) {
			printf("	- Spectral matching...\n");
			Spectral_match(excitation_voiced, glflowsp, glflowsp_new, &params);
			if(params.use_pulselib_pca == 0 && params.use_dnn_pulsegen == 0 && params.two_pitch_period_diff_pulse == 0)
				LipRadiation(excitation_voiced);
			Print_elapsed_time(&params);
		}



		/*******************************************************************/
		/*                    FILTER EXCITATION                            */
		/*******************************************************************/


		/* Combine excitations and filter */
		printf("	- Filtering excitation...\n");
		gsl_vector_add(excitation_voiced,excitation_unvoiced);
		Filter_excitation(excitation_voiced,LSF_interp,&params);
		Print_elapsed_time(&params);


		/*******************************************************************/
		/*              RESYNTHESIS FOR GAIN NORMALIZATION                 */
		/*******************************************************************/

		/* Evaluate new gain for re-synthesis */
		printf("	- Re-synthesis for gain normalization...\n");
		Evaluate_new_gain(excitation_voiced,gain_new,gain,fundf,&params);

		/* Re-synthesize with new gain */
		CreateExcitation(&params,excitation_voiced,excitation_unvoiced,fundf,gain_new,LSF,glflowsp,hnr,harmonics,waveform,h1h2,naq,
				original_pulse,pulses,pulses_rs,pulse_lengths,pgain,plsf,ptilt,phnr,pharm,pwaveform,ph1h2,pnaq,
				resynthesis_pulse_index,pulse_clus_id,pulse_clusters,gain,glflowsp_new,hnr_new,pca_mean,pca_pc,pca_w,pca_w_lib,
				stoch_env,stoch_sp,DNN_W,input_minmax,dnnpulseindices,dnnpulses);
		if(params.unvoiced_pre_emphasis == 1)
			Integrate(excitation_unvoiced, 0.97); // Slightly less than leaky factor

		/* Spectral matching */
		if((params.use_pulselib == 0 && params.use_dnn_pulsegen == 0) || (params.use_pulselib_pca == 1 && params.pca_spectral_matching == 1) ||
			 (params.use_dnn_pulsegen == 1 && params.use_dnn_specmatch == 1)) {
			Spectral_match(excitation_voiced, glflowsp, glflowsp_new, &params);
			if(params.use_pulselib_pca == 0 && params.use_dnn_pulsegen == 0 && params.two_pitch_period_diff_pulse == 0)
				LipRadiation(excitation_voiced);
		}

		/* Combine excitations and filter */
		gsl_vector_add(excitation_voiced,excitation_unvoiced);

		/* Save excitation to wav file */
		Save_excitation_to_wav(excitation_voiced,&params);

		/* Filter excitation */
		Filter_excitation(excitation_voiced,LSF_interp,&params);


		/*******************************************************************/
		/*              POSTPROCESSING AND SAVE TO FILE                    */
		/*******************************************************************/

		/* Filter out signal below F0 */
		Hp_filt_below_f0(excitation_voiced,fundf,&params);

		/* Compression of speech signal (only in noise robust speech) */
		Compression(excitation_voiced, params.compcoeff);

		/* Scale (if > 1) and save signal to file */
		Scale_signal(excitation_voiced,SCALE_IF_GREATER_THAN_ONE);
		Save_signal_to_file(excitation_voiced,&params,NULL);

		/* Report */
		Print_synthesis_settings_end(&params,time1);
	}


	/*******************************************************************/
	/*                     FREE MEMORY AND FINISH                      */
	/*******************************************************************/

	/* Free memory */
	Free_pulselib_variables(pulses,pulses_rs,pwaveform,pulse_lengths,plsf,ptilt,pharm,phnr,pgain,ph1h2,pnaq,pca_mean,pca_pc,pca_w_lib,stoch_env,stoch_sp,&params);
	Free_variables(original_pulse,excitation_voiced,excitation_unvoiced,fundf,gain,gain_new,LSF,LSF2,LSF_interp,glflowsp,glflowsp_new,
			hnr,hnr_new,harmonics,waveform,h1h2,naq,resynthesis_pulse_index,pulse_clus_id,pulse_clusters,pca_w,
			DNN_W,input_minmax,dnnpulseindices,dnnpulses,&params);

	/* Exit */
  	return EXIT_SUCCESS;
}

/***********/
/*   EOF   */
/***********/



