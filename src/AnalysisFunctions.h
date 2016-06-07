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
 * File AnalysisFunctions.h
 * Version: 1.1
 *
 */


#ifndef ANALYSISFUNCTIONS_H_
#define ANALYSISFUNCTIONS_H_



/* Parameters */
#define USE_DEF_CONF 1
#define LEAK 0.99
#define USE_PRE_EMPH 1
#define GAIN_FREQUENCY_BANDS 5
#define FFT_LENGTH 4096
#define FFT_LENGTH_LONG 8192
#define OUTPUT_FFT_LENGTH 512
#define WIN_TYPE 2
#define HARMONIC_SEARCH_COEFF 0.5
#define DEFAULT_F0 100.0
#define PROGRESS_UPDATE_INTERVAL 40
#define ROOT_SCALING 0
#define ROOT_SCALE_FACTOR 0.5
#define ROOT_SCALE_MIN_DEGREE 5
#define ROOT_SCALE_ZERO_LIMIT 0.01
#define F0_INTERP_SAMPLES 7
#define NUMBER_OF_F0_CANDIDATES 2
#define F0_FILL_RANGE 6
#define MAX_HARMONICS 1000
#define HNR_UNCERTAINTY_COEFF 0.01
#define POWER_SPECTRUM_FRAME_LEN 4096
#define HARMONIC_SEARCH_RANGE 15
#define LIP_RADIATION 0.99
#define F0_CUM_ERROR_RANGE 30
#define F0_CUM_ERROR_LIM 7.0
#define MAX_NUMBER_OF_PULSES_IN_FRAME 50
#define LSF_EPSILON 0.01
#define PULSE_START_INDEX_THRESHOLD 10.0

/* Constants */
#define VERSION "1.1"
#define DATE "20/01/2015"
#define DEF_CONF 1
#define USR_CONF 2
#define NPARAMS 15
#define HANN_WIN 1
#define BLACKMAN_WIN 2
#define HAMMING_WIN 3
#define USE_WINDOWING 1
#define NO_WINDOWING 0
#define NOTMAX -5.0
#define E_REF 0.00001
#define BIG_NEG_NUMBER -1000
#define BIG_POS_NUMBER 1000
#define SMALL_NUMBER 0.00001
#define WAV_SCALE 0.999
#define MIN_LOG_POWER -60
#define DEF_STRING_LEN 300
#define LP_METHOD_LPC "LPC"
#define LP_METHOD_WLP "WLP"
#define LP_METHOD_XLP "XLP"
#define LP_METHOD_ID_LPC 1
#define LP_METHOD_ID_WLP 2
#define LP_METHOD_ID_XLP 3
#define LP_WEIGHTING_STE "STE"
#define LP_WEIGHTING_GCI "GCI"
#define LP_WEIGHTING_ID_STE 1
#define LP_WEIGHTING_ID_GCI 2
#define FORMANT_ENH_NONE "NONE"
#define FORMANT_ENH_LSF "LSF"
#define FORMANT_ENH_LPC "LPC"
#define FORMANT_ENH_ID_NONE 0
#define FORMANT_ENH_ID_LSF 1
#define FORMANT_ENH_ID_LPC 2
#define DATA_FORMAT_ASCII "ASCII"
#define DATA_FORMAT_BINARY "BINARY"
#define DATA_FORMAT_ID_ASCII 1
#define DATA_FORMAT_ID_BINARY 2
#define DEFAULT_NRAMP 14

/* File endings */
#define FILE_ENDING_INFO ".infofile"
#define FILENAME_ENDING_LSF ".LSF"
#define FILENAME_ENDING_LSF2 ".LSF2"
#define FILENAME_ENDING_LSFSOURCE ".LSFsource"
#define FILENAME_ENDING_GAIN ".Gain"
#define FILENAME_ENDING_F0 ".F0"
#define FILENAME_ENDING_HNR ".HNR"
#define FILENAME_ENDING_HARMONICS ".Harmonics"
#define FILENAME_ENDING_WAVEFORM ".Waveform"
#define FILENAME_ENDING_SOURCESIGNAL ".SourceSignal"
#define FILENAME_ENDING_SOURCESIGNAL_WAV ".SourceSignal.wav"
#define FILENAME_ENDING_H1H2 ".H1H2"
#define FILENAME_ENDING_NAQ ".NAQ"
#define FILENAME_ENDING_FFT_VT ".FFT_VT"
#define FILENAME_ENDING_FFT_SRC ".FFT_SRC"
#define FILENAME_ENDING_PULSELIB_PULSES ".PULSELIB.pulses"
#define FILENAME_ENDING_PULSELIB_RSPULSES ".PULSELIB.rspulses"
#define FILENAME_ENDING_PULSELIB_PULSEPOS ".PULSELIB.pulsepos"
#define FILENAME_ENDING_PULSELIB_PULSEINDS ".PULSELIB.pulseinds"
#define FILENAME_ENDING_PULSELIB_PULSELENGTHS ".PULSELIB.pulselengths"
#define FILENAME_ENDING_PULSELIB_LSF ".PULSELIB.lsf"
#define FILENAME_ENDING_PULSELIB_TILT ".PULSELIB.lsfsource"
#define FILENAME_ENDING_PULSELIB_HARM ".PULSELIB.harmonics"
#define FILENAME_ENDING_PULSELIB_HNR ".PULSELIB.hnr"
#define FILENAME_ENDING_PULSELIB_GAIN ".PULSELIB.gain"
#define FILENAME_ENDING_PULSELIB_WAVEFORM ".PULSELIB.waveform"
#define FILENAME_ENDING_PULSELIB_H1H2 ".PULSELIB.h1h2"
#define FILENAME_ENDING_PULSELIB_NAQ ".PULSELIB.naq"

/* Paths to configuration file */
#define SAMPLING_FREQUENCY "SAMPLING_FREQUENCY"
#define FRAME_LENGTH "FRAME_LENGTH"
#define UNVOICED_FRAME_LENGTH "UNVOICED_FRAME_LENGTH"
#define FRAME_SHIFT "FRAME_SHIFT"
#define LPC_ORDER "LPC_ORDER"
#define LPC_ORDER_SOURCE "LPC_ORDER_SOURCE"
#define DIFFERENTIAL_LSF "DIFFERENTIAL_LSF"
#define WARPING_VT "WARPING_VT"
#define WARPING_GL "WARPING_GL"
#define HNR_CHANNELS "HNR_CHANNELS"
#define F0_FRAME_LENGTH "F0_FRAME_LENGTH"
#define NUMBER_OF_HARMONICS "NUMBER_OF_HARMONICS"
#define SEPARATE_VOICED_UNVOICED_SPECTRUM "SEPARATE_VU_SPECTRUM"
#define DATA_FORMAT "DATA_FORMAT"
#define USE_IAIF "USE_IAIF"
#define LPC_ORDER_GL_IAIF "LPC_ORDER_GL_IAIF"
#define USE_MOD_IAIF "USE_MOD_IAIF"
#define LP_METHOD "LP_METHOD"
#define LP_WEIGHTING "LP_WEIGHTING"
#define LP_STABILIZED "LP_STABILIZED"
#define FORMANT_PRE_ENH_METHOD "FORMANT_PRE_ENH_METHOD"
#define FORMANT_PRE_ENH_COEFF "FORMANT_PRE_ENH_COEFF"
#define FORMANT_PRE_ENH_LPC_DELTA "FORMANT_PRE_ENH_LPC_DELTA"
#define F0_MIN "F0_MIN"
#define F0_MAX "F0_MAX"
#define VOICING_THRESHOLD "VOICING_THRESHOLD"
#define ZCR_THRESHOLD "ZCR_THRESHOLD"
#define USE_F0_POSTPROCESSING "USE_F0_POSTPROCESSING"
#define F0_CHECK_RANGE "F0_CHECK_RANGE"
#define RELATIVE_F0_THRESHOLD "RELATIVE_F0_THRESHOLD"
#define MAX_NUMBER_OF_PULSES "MAX_NUMBER_OF_PULSES"
#define PITCH_SYNCHRONOUS_ANALYSIS "PITCH_SYNCHRONOUS_ANALYSIS"
#define PULSEMAXLEN "PULSEMAXLEN"
#define RESAMPLED_PULSELEN "RESAMPLED_PULSELEN"
#define WAVEFORM_SAMPLES "WAVEFORM_SAMPLES"
#define MAX_PULSE_LEN_DIFF "MAX_PULSE_LEN_DIFF"
#define LPC_ORDER_G "LPC_ORDER_G"
#define INVERT_SIGNAL "INVERT_SIGNAL"
#define HPFILTER_FILENAME "HPFILTER_FILENAME"
#define HP_FILTERING "HP_FILTERING"
#define EXTRACT_F0 "EXTRACT_F0"
#define EXTRACT_GAIN "EXTRACT_GAIN"
#define EXTRACT_LSF "EXTRACT_LSF"
#define EXTRACT_LSFSOURCE "EXTRACT_LSFSOURCE"
#define EXTRACT_HNR "EXTRACT_HNR"
#define EXTRACT_HARMONICS "EXTRACT_HARMONICS"
#define EXTRACT_H1H2 "EXTRACT_H1H2"
#define EXTRACT_NAQ "EXTRACT_NAQ"
#define EXTRACT_WAVEFORM "EXTRACT_WAVEFORM"
#define EXTRACT_INFOFILE "EXTRACT_INFOFILE"
#define EXTRACT_PULSELIB "EXTRACT_PULSELIB"
#define EXTRACT_SOURCE "EXTRACT_SOURCE"
#define USE_EXTERNAL_F0 "USE_EXTERNAL_F0"
#define EXTERNAL_F0_FILENAME "EXTERNAL_F0_FILENAME"
#define NOISE_REDUCTION_ANALYSIS "NOISE_REDUCTION_ANALYSIS"
#define NOISE_REDUCTION_LIMIT_DB "NOISE_REDUCTION_LIMIT_DB"
#define NOISE_REDUCTION_DB "NOISE_REDUCTION_DB"
#define LOG_F0 "LOG_F0"
#define EXTRACT_ONLY_UNIQUE_PULSES "EXTRACT_ONLY_UNIQUE_PULSES"
#define EXTRACT_ONE_PULSE_PER_FRAME "EXTRACT_ONE_PULSE_PER_FRAME"
#define EXTRACT_FFT_SPECTRA "EXTRACT_FFT_SPECTRA"
#define UNVOICED_PRE_EMPHASIS "UNVOICED_PRE_EMPHASIS"

/* Macros */
#define HANN(x,n) (0.5*(1.0-cos(2.0*M_PI*(x)/((n)-1.0))))
#define BLACKMAN(x,n) (0.42-0.5*cos(2.0*M_PI*(x)/((n)-1))+0.08*cos(4.0*M_PI*(x)/((n)-1.0)))
#define HAMMING(x,n) (0.53836 - 0.46164*(cos(2.0*M_PI*(x)/((n)-1.0))))
#define RAND() (2.0*(rand()/(double)RAND_MAX)-1.0)
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

/* Structures */
typedef struct {
	 int FS;
	 int frame_length;
	 int f0_frame_length;
	 int unvoiced_frame_length;
	 int pulsemaxlen;
	 int rspulsemaxlen;
	 int shift;
	 int n_frames;
	 int signal_length;
	 int lpc_order_vt;
	 int lpc_order_gl;
	 int differential_lsf;
	 int invert_signal;
	 int extract_pulselib_params;
	 int pitch_synchronous_analysis;
	 int f0_check_range;
	 int f0_postprocessing;
	 int maxnumberofpulses;
	 int number_of_pulses;
	 int waveform_samples;
	 int lp_method;
	 int lp_weighting;
	 int lp_stabilized;
	 int use_gif;
	 int use_iaif;
	 int lpc_order_gl_iaif;
	 int use_mod_iaif;
	 int hnr_channels;
	 int number_of_harmonics;
	 int formant_enh_method;
	 int sep_vuv_spectrum;
	 int extract_source;
	 int hp_filtering;
	 int extract_f0;
	 int extract_gain;
	 int extract_lsf;
	 int extract_tilt;
	 int extract_hnr;
	 int extract_harmonics;
	 int extract_h1h2;
	 int extract_naq;
	 int extract_waveform;
	 int extract_info;
	 int data_format;
	 int use_external_f0;
	 int noise_reduction_analysis;
	 int logf0;
	 int extract_only_unique_pulses;
	 int extract_one_pulse_per_frame;
	 int write_fftspectra_to_file;
	 int unvoiced_pre_emphasis;
	 double noise_reduction_limit_db;
	 double noise_reduction_db;
	 double frame_length_ms;
	 double unvoiced_frame_length_ms;
	 double shift_ms;
	 double pulsemaxlen_ms;
	 double rspulsemaxlen_ms;
	 double f0_frame_length_ms;
	 double lambda_vt;
	 double lambda_gl;
	 double voicing_threshold;
	 double ZCR;
	 double fmin;
	 double fmax;
	 double max_pulse_len_diff;
	 double relative_f0_threshold;
	 double formant_enh_lpc_delta;
	 double formant_enh_coeff;
	 char *hpfilter_filename;
	 char *external_f0_filename;
} PARAM;



/* Analysis functions */
int Check_command_line(int argc);
int Assign_config_parameters(struct config_t *conf, PARAM *params, int conf_type);
int Check_parameter_validity(PARAM *params);
gsl_vector *Read_soundfile(char *filename, PARAM *params);
int Read_external_f0(gsl_vector **fundf, PARAM *params);
int HighPassFilter(gsl_vector *signal, PARAM *params);
void Extract_pulses(gsl_vector *speech_frame, gsl_vector *frame_orig, gsl_vector *fundf, gsl_vector *naq, gsl_matrix *gpulses,
		gsl_matrix *gpulses_rs, gsl_vector *pulse_pos, gsl_vector *pulse_inds, gsl_vector *pulse_lengths, gsl_matrix *plsf,
		gsl_matrix *ptilt, gsl_matrix *pharm, gsl_matrix *phnr, gsl_vector *pgain, gsl_vector *ph1h2, gsl_vector *pnaq,
		gsl_matrix *lsf, gsl_matrix *spectral_tilt,	gsl_matrix *harmonics, gsl_vector *h1h2, gsl_matrix *hnr_i,gsl_vector *gain,
		gsl_matrix *pwaveform, gsl_matrix *waveform, int index, PARAM *params);
void Print_progress(int index, int n_frames, char *filename, int FS);
void Get_samples_to_frames(gsl_vector *signal, gsl_vector *frame, gsl_vector *frame0, gsl_vector *f0_frame, gsl_vector *f0_frame0, int shift, int index);
double Define_current_f0(gsl_vector *fundf, double f0, int index);
gsl_vector *Find_GCIs(gsl_vector *frame_orig, gsl_vector *fundf, int FS, int index);
void Gain(gsl_vector *frame, gsl_vector *gain, int index, int windowing);
void UnvoicedGain(gsl_vector *frame, gsl_vector *gain, int index, int windowing, int unvoiced_frame_length);
void BandPassGain(gsl_vector *frame, gsl_matrix *bp_gain, PARAM *params, int index);
gsl_vector *Conv(gsl_vector *conv1, gsl_vector *conv2);
void Interpolate(gsl_vector *vector, gsl_vector *i_vector);
void lsf2poly(gsl_matrix *lsf_matrix, gsl_vector *poly, int index, int HMM);
void lsf_vector2poly(gsl_vector *lsf_vector, gsl_vector *poly);
void InverseFiltering(gsl_vector *frame, gsl_vector *frame0, gsl_vector *glottal, gsl_matrix *LSF, gsl_matrix *LSF2,
		gsl_matrix *spectral_tilt, gsl_vector *fundf, gsl_vector *glottsig, gsl_matrix *fftmatrix_vt,
		gsl_matrix *fftmatrix_src, gsl_matrix *fftmatrix_uv, int index, PARAM* params);
void InverseFiltering_long(gsl_vector *frame,gsl_vector *frame0,gsl_vector *glottal, gsl_vector *fundf, gsl_vector *glottsig, int index, PARAM *params);
int EvalFileLength(const char *name);
void Allocate_variables(gsl_vector **frame, gsl_vector **frame0, gsl_vector **glottal, gsl_vector **gain, gsl_vector **uvgain,
		gsl_vector **f0_frame, gsl_vector **glottal_f0, gsl_vector **f0_frame0, gsl_vector **glottsig, gsl_vector **glottsig_f0, gsl_vector **source_signal,
		gsl_matrix **fundf_candidates, gsl_matrix **LSF, gsl_matrix **LSF2, gsl_matrix **bp_gain, gsl_matrix **spectral_tilt, gsl_matrix **HNR,
		gsl_matrix **waveform, gsl_matrix **harmonics, gsl_vector **h1h2, gsl_vector **naq, gsl_matrix **fftmatrix_vt, gsl_matrix **fftmatrix_src,  gsl_matrix **fftmatrix_uv, PARAM *params);
void Allocate_pulselib_variables(gsl_matrix **gpulses, gsl_matrix **gpulses_rs, gsl_vector **pulse_inds, gsl_vector **pulse_pos, gsl_vector **pulse_lengths, gsl_matrix **plsf,
		gsl_matrix **ptilt, gsl_matrix **pharm, gsl_matrix **phnr, gsl_matrix **pwaveform, gsl_vector **pgain, gsl_vector **ph1h2, gsl_vector **pnaq, PARAM *params);
void Integrator(gsl_vector *frame, double leak);
void LPC(gsl_vector *frame, gsl_vector *a);
void WLPC(gsl_vector *frame, gsl_vector *a, double lambda);
void SWLP(gsl_vector *frame, gsl_vector *a, int M, int lag, int lp_weighting, gsl_vector *fundf, int FS, int index, gsl_vector *glottsig, int stabilized);
void SXLP(gsl_vector *frame, gsl_vector *a, int w, int stabilized);
void Filter(gsl_vector *frame, gsl_vector *glottal, gsl_vector *a);
void RealRootScale(gsl_vector *a);
void Convert_matrix_to_LSF(gsl_matrix *LSF, gsl_vector *a, int index);
void Convert_vector_to_LSF(gsl_vector *a, gsl_vector *LSF);
void Chebyshev(gsl_vector *T);
void MedFilt3(gsl_vector *frame);
void MedFilt3_matrix(gsl_matrix *matrix);
void MedFilt5_matrix(gsl_matrix *matrix);
void MedFilt5(gsl_vector *frame);
void Remove_mean(gsl_vector *frame);
void FundF(gsl_vector *frame, gsl_vector *signal, gsl_matrix *fundf_candidates,	gsl_vector *fundf, gsl_matrix *bp_gain, int index, PARAM *params);
double Parabolic_interpolation(gsl_vector *r, int maxind, PARAM *params);
void F0_postprocess(gsl_vector *fundf, gsl_matrix *fundf_candidates, PARAM *params);
void Fundf_postprocessing(gsl_vector *fundf, gsl_vector *fundf_orig, gsl_matrix *fundf_candidates, PARAM *params);
void Fill_f0_gaps(gsl_vector *fundf, PARAM *params);
struct config_t *Read_config(char *filename);
void AllPassDelay(gsl_vector *signal, double lambda);
void WFilter(gsl_vector *signal, gsl_vector *result, gsl_vector *A, gsl_vector *B, double lambda);
double *NFArray(int size);
void alphas2sigmas(double *alp, double *sigm, double lambda, int dim);
void Harmonic_analysis(gsl_vector *frame, gsl_matrix *harmonics, gsl_matrix *hnr, gsl_vector *h1h2, double f0, int FS, int index);
void Convert_Hz2ERB(gsl_vector *vector, gsl_vector *vector_erb, int FS);
void Differentiate(gsl_vector *signal, double leak);
void Differentiate_noleak(gsl_vector *signal);
void Differentiate_LSFs(gsl_matrix **m1);
gsl_matrix *Construct_C(int size);
void Add_missing_frames(gsl_vector **fundf,gsl_matrix **fundf_candidates,gsl_vector **gain,gsl_vector **uvgain,gsl_matrix **hnr,
		gsl_matrix **spectral_tilt,gsl_matrix **LSF,gsl_matrix **LSF2,gsl_matrix **harmonics,gsl_matrix **waveform,
		gsl_vector **h1h2, gsl_vector **naq, gsl_matrix **fftmatrix_vt, gsl_matrix **fftmatrix_src, gsl_matrix **fftmatrix_uv, PARAM *params);
void LSF_Postfilter(gsl_matrix *lsf, PARAM *params);
void LPC_Postfilter(gsl_matrix *LSF, PARAM *params);
void Pre_emphasis(gsl_vector *signal, double coef);
void ModPowerSpectrum(gsl_vector *s, gsl_vector *formants, int n, double gamma, double delta);
void GetFormants(gsl_vector *s, gsl_vector *formants, int *n);
void Select_new_refined_values(gsl_vector *fundf, gsl_vector *gain, gsl_matrix *LSF, gsl_matrix *spectral_tilt,
		gsl_matrix *hnr_i, gsl_matrix *harmonics, gsl_vector *pgain, gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *phnr,
		gsl_matrix *pharm, gsl_vector *pulse_pos, gsl_vector *pulse_lengths, gsl_matrix *waveform,
		gsl_vector *h1h2, gsl_vector *ph1h2, gsl_vector *naq, gsl_vector *pnaq, PARAM *params);
void Fill_waveform_gaps(gsl_vector *fundf, gsl_matrix *waveform);
void Fill_naq_gaps(gsl_vector *fundf, gsl_vector *naq);
void Eval_STE_weight(gsl_vector *frame, gsl_vector *ste, int M, int lag);
void Eval_GCI_weight1(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index);
void Eval_GCI_weight2(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index);
void Eval_GCI_weight3(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index);
void Construct_source(gsl_vector *source_signal, gsl_vector *glottal, gsl_vector *glottal_f0, int index, PARAM *params);
void Write_parameters_to_file(char *filename, gsl_matrix *LSF, gsl_matrix *LSF2, gsl_matrix *spectral_tilt, gsl_matrix *HNR, gsl_matrix* harmonics, gsl_matrix *waveform,
		gsl_vector *fundf, gsl_vector *gain, gsl_vector *h1h2, gsl_vector *naq, gsl_vector *source_signal, gsl_matrix *fftmatrix_vt, gsl_matrix *fftmatrix_src, PARAM *params);
void Write_pulselibrary_to_file(char *filename, gsl_matrix *gpulses, gsl_matrix *gpulses_rs, gsl_vector *pulse_lengths, gsl_vector *pulse_pos, gsl_vector *pulse_inds,
		gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *phnr, gsl_matrix *pharm, gsl_matrix *pwaveform, gsl_vector *pgain,
		gsl_vector *ph1h2, gsl_vector *pnaq, PARAM *params);
void Free_variables(gsl_vector *frame, gsl_vector *frame0, gsl_vector *signal, gsl_vector *glottal, gsl_vector *glottsig, gsl_vector *glottsig_f0,
		gsl_vector *uvgain, gsl_vector *f0_frame, gsl_vector *f0_frame0, gsl_vector *glottal_f0, gsl_matrix *bp_gain, gsl_matrix *fundf_candidates,
		gsl_matrix *fftmatrix_uv);
void Save_source_signal_to_wavfile(char *filename, gsl_vector *signal, PARAM *params);
void Noise_reduction(gsl_vector *gain, PARAM *params);
void Convert_F0_to_log(gsl_vector *fundf, PARAM *params);
void Select_unique_pulses(gsl_matrix **gpulses, gsl_matrix **gpulses_rs, gsl_vector **pulse_lengths, gsl_vector **pulse_pos, gsl_vector **pulse_inds,
		gsl_matrix **plsf, gsl_matrix **ptilt, gsl_matrix **phnr, gsl_matrix **pharm, gsl_matrix **pwaveform, gsl_vector **pgain,
		gsl_vector **ph1h2, gsl_vector **pnaq, PARAM *params);
void Select_one_pulse_per_frame(gsl_matrix **gpulses, gsl_matrix **gpulses_rs, gsl_vector **pulse_lengths, gsl_vector **pulse_pos, gsl_vector **pulse_inds,
		gsl_matrix **plsf, gsl_matrix **ptilt, gsl_matrix **phnr, gsl_matrix **pharm, gsl_matrix **pwaveform, gsl_vector **pgain,
		gsl_vector **ph1h2, gsl_vector **pnaq, PARAM *params);
void LSF_stabilize(gsl_vector *a);
void Pole_stabilize(gsl_vector *a);
void EvalFFTSpectrum(gsl_vector *signal, gsl_matrix *matrix, int index, PARAM *params);
void WriteFFTSpectrumToFile(char *filename, gsl_matrix *fft, PARAM *params);

/* Currently not in use */
void Interpolate_poly(gsl_vector *vector, gsl_vector *i_vector);
void Interpolate_lin(gsl_vector *vector, gsl_vector *i_vector);
void LSF_fix_vector(gsl_vector *lsf);
void IIR_filter(gsl_vector *signal, gsl_vector *coeffs);
void FIR_filter_file(gsl_vector *signal, char *filename, int n);
void Estimate_ERB_gain(gsl_vector *frame, gsl_matrix *erb_gain, int index, int FS);

/* Test functions */
void pause(int print);
void VPrint1(gsl_vector *vector);
void VPrint2(gsl_vector *vector);
void VPrint3(gsl_vector *vector);
void VPrint4(gsl_vector *vector);
void VPrint5(gsl_vector *vector);
void VPrint6(gsl_vector *vector);
void VPrint7(gsl_vector *vector);
void VPrint8(gsl_vector *vector);
void VPrint9(gsl_vector *vector);
void VPrint10(gsl_vector *vector);
void APrint1(double *array, int size);
void APrint2(double *array, int size);
void APrint3(double *array, int size);
void APrint4(double *array, int size);
void MPrint1(gsl_matrix *matrix);
void MPrint2(gsl_matrix *matrix);
void MPrint3(gsl_matrix *matrix);
void MPrint4(gsl_matrix *matrix);
int Find_matrix_NaN_Inf(gsl_matrix *m);
int Find_vector_NaN_Inf(gsl_vector *v);

#endif

