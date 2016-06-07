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
 * File Analysis.c
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
#include <sndfile.h> 				/* Read and write wav */
#include <gsl/gsl_vector.h>		/* GSL, Vector */
#include <gsl/gsl_matrix.h>		/* GSL, Matrix */
#include <gsl/gsl_fft_real.h>		/* GSL, FFT */
#include <libconfig.h>				/* Configuration file */
#include "AnalysisFunctions.h"







/*******************************************************************/
/*                          MAIN                                   */
/*******************************************************************/

int main(int argc, char *argv[]) {

	/* Check command line format */
	if(Check_command_line(argc) == EXIT_FAILURE)
		return EXIT_FAILURE;

	/* Read default configuration file and assign parameters */
	struct config_t *conf_def = Read_config(argv[2]);
	if(conf_def == NULL)
		return EXIT_FAILURE;
	PARAM params;
	if(Assign_config_parameters(conf_def,&params,DEF_CONF) == EXIT_FAILURE) {
		config_destroy(conf_def);
		free(conf_def);
		return EXIT_FAILURE;
	}

	/* Read user configuration file and assign parameters */
	if(argc == 4) {
		struct config_t *conf_usr = Read_config(argv[3]);
		if(conf_usr == NULL)
			return EXIT_FAILURE;
		if(Assign_config_parameters(conf_usr,&params,USR_CONF) == EXIT_FAILURE) {
			config_destroy(conf_usr);
			free(conf_usr);
			return EXIT_FAILURE;
		}
	}

	/* Check the validity of the parameters */
	if(Check_parameter_validity(&params) == EXIT_FAILURE)
		return EXIT_FAILURE;

	/* Read soundfile, define sampling frequency, invert signal if requested */
	gsl_vector *signal = Read_soundfile(argv[1], &params);
	if(signal == NULL)
		return EXIT_FAILURE;

	/* Read external F0 file if requested */
	gsl_vector *fundf;
	if(Read_external_f0(&fundf,&params) == EXIT_FAILURE) {
		gsl_vector_free(signal);
		return EXIT_FAILURE;
	}

	/* High-pass filtering */
	if(HighPassFilter(signal,&params) == EXIT_FAILURE) {
		gsl_vector_free(signal);
		return EXIT_FAILURE;
	}

	/*******************************************************************/
	/*                    ALLOCATE MEMORY                              */
	/*******************************************************************/

	/* Allocate memory for variables */
	int i,j,index;
	double f0 = 0;
	gsl_vector *frame,*frame0,*glottal,*gain,*uvgain,*f0_frame,*glottal_f0,*f0_frame0,*glottsig,*glottsig_f0;
	gsl_vector *source_signal,*h1h2,*naq,*ph1h2,*pnaq;
	gsl_matrix *fundf_candidates, *LSF, *LSF2, *bp_gain, *spectral_tilt, *HNR, *waveform, *harmonics,
			*fftmatrix_vt, *fftmatrix_src, *fftmatrix_uv;
	Allocate_variables(&frame,&frame0,&glottal,&gain,&uvgain,&f0_frame,&glottal_f0,&f0_frame0,&glottsig,
			&glottsig_f0,&source_signal,&fundf_candidates,&LSF,&LSF2,&bp_gain,&spectral_tilt,&HNR,
			&waveform,&harmonics,&h1h2,&naq,&fftmatrix_vt,&fftmatrix_src,&fftmatrix_uv,&params);


	/*******************************************************************/
	/*                   ALLOCATE PULSE LIBRARY                        */
	/*******************************************************************/

	/* Allocate memory for pulse library variables */
	gsl_vector *pulse_inds,*pulse_lengths,*pulse_pos,*pgain;
	gsl_matrix *gpulses,*gpulses_rs,*plsf,*ptilt,*pharm,*phnr,*pwaveform;
	Allocate_pulselib_variables(&gpulses,&gpulses_rs,&pulse_inds,&pulse_pos,&pulse_lengths,
			&plsf,&ptilt,&pharm,&phnr,&pwaveform,&pgain,&ph1h2,&pnaq,&params);


	/********************************************************************/
	/*                    EXTRACT PARAMETERS                            */
	/********************************************************************/

	/* Analysis: Start reading signal vector */
	for(index=0; index<params.n_frames; index++) {

		/* Print progress */
		Print_progress(index,params.n_frames,argv[1],params.FS);

		/* Get samples to frames */
		Get_samples_to_frames(signal,frame,frame0,f0_frame,f0_frame0,params.shift,index);

		/* Gain extraction */
		Gain(frame, gain, index, USE_WINDOWING);
		UnvoicedGain(frame, uvgain, index, USE_WINDOWING, params.unvoiced_frame_length);
		BandPassGain(frame, bp_gain, &params, index);

		/* Preliminary inverse filtering for GCI detection (e.g. in GCI-weighted SWLP) */
		if(params.lp_method == LP_METHOD_ID_WLP && params.lp_weighting == LP_WEIGHTING_ID_GCI) {
			InverseFiltering_long(f0_frame, f0_frame0, glottsig_f0, fundf, glottsig, index, &params);
			FundF(glottsig_f0, frame, fundf_candidates, fundf, bp_gain, index, &params);
			InverseFiltering_long(frame, frame0, glottsig, fundf, glottsig, index, &params);
		}

		/* Inverse filtering */
		InverseFiltering(frame, frame0, glottal, LSF, LSF2, spectral_tilt, fundf, glottsig, fftmatrix_vt,
				fftmatrix_src, fftmatrix_uv, index, &params);
		InverseFiltering_long(f0_frame, f0_frame0, glottal_f0, fundf, glottsig_f0, index, &params);

		/* Fundamental frequency estimation. Define a f0 value for the frame (even if unvoiced) */
		FundF(glottal_f0, frame, fundf_candidates, fundf, bp_gain, index, &params);
		f0 = Define_current_f0(fundf,f0,index);

		/* Estimate Harmonic-to-Noise Ratio (HNR) and magnitudes of the first N harmonics */
		Harmonic_analysis(glottal_f0, harmonics, HNR, h1h2, f0, params.FS, index);

		/* Extract pulses, pitch-synchronous spectrum estimation */
		Extract_pulses(f0_frame,glottal_f0,fundf,naq,gpulses,gpulses_rs,pulse_pos,pulse_inds,pulse_lengths,
				plsf,ptilt,pharm,phnr,pgain,ph1h2,pnaq,LSF,spectral_tilt,harmonics,h1h2,HNR,gain,pwaveform,waveform,index,&params);

		/* Construct source signal */
		Construct_source(source_signal, glottal, glottal_f0, index, &params);
	}


	/********************************************************************/
	/*                POSTPROCESSING OF PARAMETERS                      */
	/********************************************************************/

	/* Add missing frames to the beginning and the end of parameter vectors/matrices.
	 * This is due to the analysis scheme, where the analysis starts
	 * and ends with whole frames instead of zero-padding */
	Add_missing_frames(&fundf,&fundf_candidates,&gain,&uvgain,&HNR,&spectral_tilt,&LSF,&LSF2,&harmonics,&waveform,&naq,&h1h2,
			&fftmatrix_vt,&fftmatrix_src,&fftmatrix_uv,&params);

	/* Postprocess F0 */
	F0_postprocess(fundf,fundf_candidates,&params);

	/* Median filtering */
	MedFilt5_matrix(HNR);
	MedFilt5_matrix(harmonics);
	MedFilt5(h1h2);
	MedFilt5(naq);

	/* Replace unvoiced fetures */
	for(i=0; i<fundf->size; i++) {
		if(gsl_vector_get(fundf, i) == 0) {
			if(params.sep_vuv_spectrum == 0) {
				for(j=0; j<params.lpc_order_vt; j++)
					gsl_matrix_set(LSF, i, j, gsl_matrix_get(LSF2, i, j));
			}
			for(j=0; j<fftmatrix_vt->size2; j++)
				gsl_matrix_set(fftmatrix_vt, i, j, gsl_matrix_get(fftmatrix_uv, i, j));
			gsl_vector_set(gain,i,gsl_vector_get(uvgain,i));
		}
	}

	/* Noise reduction */
	Noise_reduction(gain,&params);

	/* Fill possible gaps in waveform and NAQ */
	Fill_waveform_gaps(fundf,waveform);
	Fill_naq_gaps(fundf,naq);

	/* Select new values to pulse parameters according to refined parameters */
	Select_new_refined_values(fundf,gain,LSF,spectral_tilt,HNR,harmonics,pgain,plsf,ptilt,
			phnr,pharm,pulse_pos,pulse_lengths,waveform,h1h2,ph1h2,naq,pnaq,&params);

	/* Select only unique pulses */
	Select_unique_pulses(&gpulses,&gpulses_rs,&pulse_lengths,&pulse_pos,&pulse_inds,&plsf,
			&ptilt,&phnr,&pharm,&pwaveform,&pgain,&ph1h2,&pnaq,&params);

	/* Select only one pulse per frame */
	Select_one_pulse_per_frame(&gpulses,&gpulses_rs,&pulse_lengths,&pulse_pos,&pulse_inds,&plsf,
				&ptilt,&phnr,&pharm,&pwaveform,&pgain,&ph1h2,&pnaq,&params);


	/********************************************************************/
	/*                      FORMANT ENHANCEMENT                         */
	/********************************************************************/

	/* Formant enhancement by LSFs*/
	if(params.formant_enh_method == FORMANT_ENH_ID_LSF)
		LSF_Postfilter(LSF, &params);

	/* Formant Enhancement by re-estimating LPC and modifying the autocorrelation */
	if(params.formant_enh_method == FORMANT_ENH_ID_LPC)
		LPC_Postfilter(LSF, &params);

	/* Differentiate LSFs if requested */
	if(params.differential_lsf == 1) {
		Differentiate_LSFs(&LSF);
		Differentiate_LSFs(&LSF2);
		Differentiate_LSFs(&spectral_tilt);
	}

	/* Convert F0 to logarithmic scale */
	Convert_F0_to_log(fundf,&params);


	/********************************************************************/
	/*                    WRITE PARAMETERS TO FILE                      */
	/********************************************************************/

	/* Open files for writing parameters, free memory */
	Write_parameters_to_file(argv[1],LSF,LSF2,spectral_tilt,HNR,harmonics,waveform,fundf,gain,h1h2,naq,source_signal,
			fftmatrix_vt,fftmatrix_src,&params);
	Free_variables(frame,frame0,signal,glottal,glottsig,glottsig_f0,uvgain,f0_frame,f0_frame0,glottal_f0,bp_gain,
			fundf_candidates,fftmatrix_uv);


	/**************************************************************/
	/*            WRITE  PULSE LIBRARY TO FILE                    */
	/**************************************************************/

	/* Select unique pulses and write pulse library data to file, free memory */
	Write_pulselibrary_to_file(argv[1],gpulses,gpulses_rs,pulse_lengths,pulse_pos,pulse_inds,plsf,
			ptilt,phnr,pharm,pwaveform,pgain,ph1h2,pnaq,&params);

	/* Finish */
	printf("\nFinished analysis.\n\n");
	return EXIT_SUCCESS;
}

/***********/
/*   EOF   */
/***********/

