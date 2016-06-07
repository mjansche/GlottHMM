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
 * File SynthesisFunctions.h
 * Version: 1.1
 *
 */


#ifndef SYNTHESISFUNCTIONS_H_
#define SYNTHESISFUNCTIONS_H_

/* Parameters */
#define USE_DEF_CONF 1
#define LIP_RADIATION 0.99
#define USE_PRE_EMPH 1
#define FREQUENCY_BANDS 5
#define MIN_FFT_LENGTH 4096
#define WIN_TYPE 2
#define LEAK 0.99
#define WAV_SCALE 0.999
#define MAX_HARMONICS 300
#define HARMONIC_SEARCH_COEFF 0.5
#define HNR_UNCERTAINTY_COEFF 0.01
#define MIN_LOG_POWER (-60)
#define POWER_SPECTRUM_WIN 20
#define POWER_SPECTRUM_FRAME_LEN 4096
#define LSF_EPSILON 0.01

/* Constants */
#define VERSION "1.1"
#define DATE "20/01/2015"
#define DEF_CONF 1
#define USR_CONF 2
#define HANN_WIN 1
#define BLACKMAN_WIN 2
#define HAMMING_WIN 3
#define E_REF 0.00001
#define MIN_P_TILT 3
#define NPARAMS 15
#define BIG_NEG_NUMBER (-1000)
#define BIG_POS_NUMBER 1000
#define BIGGER_POS_NUMBER 1000000
#define SMALL_NUMBER 0.00001
#define USE_WINDOWING 1
#define NO_WINDOWING 0
#define SCALE_IF_GREATER_THAN_ONE 0
#define FORCE_MAX_TO_ONE 1
#define NO_SMOOTHING 1
#define DEF_STRING_LEN 300
#define POSTFILTER_NONE "NONE"
#define POSTFILTER_LSF "LSF"
#define POSTFILTER_LPC "LPC"
#define POSTFILTER_ID_NONE 0
#define POSTFILTER_ID_LSF 1
#define POSTFILTER_ID_LPC 2
#define DATA_FORMAT_ASCII "ASCII"
#define DATA_FORMAT_BINARY "BINARY"
#define DATA_FORMAT_ID_ASCII 1
#define DATA_FORMAT_ID_BINARY 2
#define STOCH_SP_LPC_ORDER 31

/* File endings */
#define FILENAME_ENDING_SYNTHESIS ".syn.wav"
#define FILENAME_ENDING_EXCITATION ".exc.wav"
#define FILENAME_ENDING_INFO ".infofile"
#define FILENAME_ENDING_LSF ".LSF"
#define FILENAME_ENDING_LSF2 ".LSF2"
#define FILENAME_ENDING_LSFSOURCE ".LSFsource"
#define FILENAME_ENDING_GAIN ".Gain"
#define FILENAME_ENDING_F0 ".F0"
#define FILENAME_ENDING_HNR ".HNR"
#define FILENAME_ENDING_HARMONICS ".Harmonics"
#define FILENAME_ENDING_WAVEFORM ".Waveform"
#define FILENAME_ENDING_H1H2 ".H1H2"
#define FILENAME_ENDING_NAQ ".NAQ"
#define FILENAME_ENDING_PULSELIB_PULSES ".pulses"
#define FILENAME_ENDING_PULSELIB_RSPULSES ".rspulses"
#define FILENAME_ENDING_PULSELIB_PULSELENGTHS ".pulselengths"
#define FILENAME_ENDING_PULSELIB_LSF ".lsf"
#define FILENAME_ENDING_PULSELIB_TILT ".lsfsource"
#define FILENAME_ENDING_PULSELIB_HARM ".harmonics"
#define FILENAME_ENDING_PULSELIB_HNR ".hnr"
#define FILENAME_ENDING_PULSELIB_WAVEFORM ".waveform"
#define FILENAME_ENDING_PULSELIB_GAIN ".gain"
#define FILENAME_ENDING_PULSELIB_H1H2 ".h1h2"
#define FILENAME_ENDING_PULSELIB_NAQ ".naq"
#define FILENAME_ENDING_PULSELIB_STOCH_ENV ".stoch_env"
#define FILENAME_ENDING_PULSELIB_STOCH_SP ".stoch_sp"
#define FILENAME_ENDING_PULSECLUSTER ".pulsecl"

/* Paths to configuration file */
#define SAMPLING_FREQUENCY "SAMPLING_FREQUENCY"
#define FRAME_LENGTH "FRAME_LENGTH"
#define FRAME_SHIFT "FRAME_SHIFT"
#define LPC_ORDER "LPC_ORDER"
#define LPC_ORDER_SOURCE "LPC_ORDER_SOURCE"
#define WARPING_VT "WARPING_VT"
#define WARPING_GL "WARPING_GL"
#define HNR_CHANNELS "HNR_CHANNELS"
#define F0_FRAME_LENGTH "F0_FRAME_LENGTH"
#define NUMBER_OF_HARMONICS "NUMBER_OF_HARMONICS"
#define SEPARATE_VOICED_UNVOICED_SPECTRUM "SEPARATE_VU_SPECTRUM"
#define DIFFERENTIAL_LSF "DIFFERENTIAL_LSF"
#define DATA_FORMAT "DATA_FORMAT"
#define GLOTTAL_PULSE_NAME "GLOTTAL_PULSE_NAME"
#define PULSE_LIBRARY_NAME "PULSE_LIBRARY_NAME"
#define USE_PULSE_LIBRARY "USE_PULSE_LIBRARY"
#define SYNTHESIZE_MULTIPLE_FILES "SYNTHESIZE_MULTIPLE_FILES"
#define SYNTHESIS_LIST "SYNTHESIS_LIST"
#define NUMBER_OF_PULSE_CANDIDATES "NUMBER_OF_PULSE_CANDIDATES"
#define USE_PULSE_INTERPOLATION "USE_PULSE_INTERPOLATION"
#define CONCATENATION_COST "CONCATENATION_COST"
#define TARGET_COST "TARGET_COST"
#define ADD_NOISE_PULSELIB "ADD_NOISE_PULSELIB"
#define PULSE_ERROR_BIAS "PULSE_ERROR_BIAS"
#define USE_PULSE_CLUSTERING "USE_PULSE_CLUSTERING"
#define MAX_PULSES_IN_CLUSTER "MAX_PULSES_IN_CLUSTER"
#define PARAMETER_WEIGHTS "PARAMETER_WEIGHTS"
#define NOISE_GAIN_VOICED "NOISE_GAIN_VOICED"
#define USE_HMM "USE_HMM"
#define GAIN_UNVOICED "GAIN_UNVOICED"
#define NORM_GAIN_SMOOTH_V_LEN "NORM_GAIN_SMOOTH_V_LEN"
#define NORM_GAIN_SMOOTH_UV_LEN "NORM_GAIN_SMOOTH_UV_LEN"
#define GAIN_UNVOICED_FRAME_LENGTH "GAIN_UNVOICED_FRAME_LENGTH"
#define GAIN_VOICED_FRAME_LENGTH "GAIN_VOICED_FRAME_LENGTH"
#define FILTER_UPDATE_INTERVAL_VT "FILTER_UPDATE_INTERVAL_VT"
#define FILTER_UPDATE_INTERVAL_GL "FILTER_UPDATE_INTERVAL_GL"
#define LSF_SMOOTH_LEN "LSF_SMOOTH_LEN"
#define HARMONICS_SMOOTH_LEN "HARMONICS_SMOOTH_LEN"
#define GAIN_SMOOTH_LEN "GAIN_SMOOTH_LEN"
#define HNR_SMOOTH_LEN "HNR_SMOOTH_LEN"
#define GLFLOWSP_SMOOTH_LEN "LSFSOURCE_SMOOTH_LEN"
#define PITCH "PITCH"
#define SPEED "SPEED"
#define NOISE_LOW_FREQ_LIMIT "NOISE_LOW_FREQ_LIMIT"
#define USE_HARMONIC_MODIFICATION "USE_HARMONIC_MODIFICATION"
#define HP_FILTER_F0 "HP_FILTER_F0"
#define NOISE_ROBUST_SPEECH "NOISE_ROBUST_SPEECH"
#define POSTFILTER_METHOD "POSTFILTER_METHOD"
#define POSTFILTER_COEFFICIENT "POSTFILTER_COEFFICIENT"
#define WAVEFORM_SAMPLES "WAVEFORM_SAMPLES"
#define USE_TILT "USE_LSFSOURCE"
#define USE_HNR "USE_HNR"
#define USE_HARMONICS "USE_HARMONICS"
#define USE_H1H2 "USE_H1H2"
#define USE_NAQ "USE_NAQ"
#define USE_WAVEFORM "USE_WAVEFORM"
#define WRITE_EXCITATION_TO_WAV "WRITE_EXCITATION_TO_WAV"
#define JITTER "JITTER"
#define NOISE_REDUCTION_SYNTHESIS "NOISE_REDUCTION_SYNTHESIS"
#define NOISE_REDUCTION_LIMIT_DB "NOISE_REDUCTION_LIMIT_DB"
#define NOISE_REDUCTION_DB "NOISE_REDUCTION_DB"
#define NORMALIZE_PULSELIB "NORMALIZE_PULSELIB"
#define ADAPT_TO_PULSELIB "ADAPT_TO_PULSELIB"
#define ADAPT_COEFF "ADAPT_COEFF"
#define USE_PULSELIB_LSF "USE_PULSELIB_LSF"
#define USE_PULSE_PCA "USE_PULSE_PCA"
#define PCA_SPECTRAL_MATCHING "PCA_SPECTRAL_MATCHING"
#define AVERAGE_N_ADJACENT_PULSES "AVERAGE_N_ADJACENT_PULSES"
#define LOG_F0 "LOG_F0"
#define HNR_COMPENSATION "HNR_COMPENSATION"
#define DNN_INPUT_NORMALIZED "DNN_INPUT_NORMALIZED"
#define DNN_NUMBER_OF_STACKED_FRAMES "DNN_NUMBER_OF_STACKED_FRAMES"
#define UNVOICED_PRE_EMPHASIS "UNVOICED_PRE_EMPHASIS"
#define TWO_PITCH_PERIOD_DIFF_PULSE "TWO_PITCH_PERIOD_DIFF_PULSE"

/* Pulse library PCA */
#define USE_PULSELIB_PCA "USE_PULSELIB_PCA"
#define PCA_ORDER "PCA_ORDER"
#define PCA_ORDER_SYNTHESIS "PCA_ORDER_SYNTHESIS"
#define PCA_PULSE_LENGTH "PCA_PULSE_LENGTH"
#define FILENAME_ENDING_PULSELIB_PCA_PC ".pca_pc"
#define FILENAME_ENDING_PULSELIB_PCA_MEAN ".pca_mean"
#define FILENAME_ENDING_PULSELIB_PCA_W ".pca_w"

/* DNN pulse generation */
#define USE_DNN_PULSEGEN "USE_DNN_PULSEGEN"
#define USE_DNN_PULSELIB_SEL "USE_DNN_PULSELIB_SEL"
#define USE_DNN_SPECMATCH "USE_DNN_SPECMATCH"
#define DNN_WEIGHT_DIMS "DNN_WEIGHT_DIMS"
#define DNN_WEIGHT_PATH "DNN_WEIGHT_PATH"
#define FILENAME_DNN_W "dnn.w"
#define FILENAME_INPUT_MINMAX "input.minmax"

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
	 int gain_voiced_frame_length;
	 int gain_unvoiced_frame_length;
	 int shift;
	 int pulsemaxlen;
	 int rspulsemaxlen;
	 int number_of_pulses;
	 int signal_length;
	 int n_frames;
	 int lpc_order_vt;
	 int lpc_order_gl;
	 int use_hmm;
	 int use_pulselib;
	 int filter_update_interval_vt;
	 int filter_update_interval_gl;
	 int glflowsp_smooth_len;
	 int hnr_channels;
	 int number_of_harmonics;
	 int norm_gain_smooth_v_len;
	 int norm_gain_smooth_uv_len;
	 int lsf_smooth_len;
	 int gain_smooth_len;
	 int hnr_smooth_len;
	 int max_pulses_in_cluster;
	 int use_harmonic_modification;
	 int hpfiltf0;
	 int postfilter_modulation;
	 int harmonics_smooth_len;
	 int noise_robust_speech;
	 int postfilter_method;
	 int sep_vuv_spectrum;
	 int n_pulsecandidates;
	 int resynth;
	 int pulse_clustering;
	 int pulse_interpolation;
	 int differential_lsf;
	 int add_noise_pulselib;
	 int multisyn;
	 int synlistlen;
	 int synfilenumber;
	 int hnr_reestimated;
	 int waveform_samples;
	 int use_tilt;
	 int use_hnr;
	 int use_harmonics;
	 int use_h1h2;
	 int use_naq;
	 int use_waveform;
	 int data_format;
	 int write_excitation_to_wav;
	 int noise_reduction_synthesis;
	 int normalize_pulselib;
	 int adapt_to_pulselib;
	 int use_pulselib_lsf;
	 int use_pulselib_pca;
	 int pca_order;
	 int pca_order_synthesis;
	 int pca_spectral_matching;
	 int pca_pulse_length;
	 int use_pulse_pca;
	 int average_n_adjacent_pulses;
	 int logf0;
	 int use_dnn_pulsegen;
	 int use_dnn_pulselib_sel;
	 int use_dnn_specmatch;
	 int hnr_compensation;
	 int dnn_input_normalized;
	 int dnn_number_of_stacked_frames;
	 int unvoiced_pre_emphasis;
	 int two_pitch_period_diff_pulse;
	 double adapt_coeff;
	 double noise_reduction_limit_db;
	 double noise_reduction_db;
	 double jitter;
	 double filter_update_interval_vt_ms;
	 double filter_update_interval_gl_ms;
	 double time_temp;
	 double frame_length_ms;
	 double gain_voiced_frame_length_ms;
	 double gain_unvoiced_frame_length_ms;
	 double shift_ms;
	 double pulsemaxlen_ms;
	 double rspulsemaxlen_ms;
	 double f0_frame_length_ms;
	 double pitch;
	 double speed;
	 double postfilter_alpha;
	 double absmax;
	 double rho;
	 double lambda_vt;
	 double lambda_gl;
	 double gain_unvoiced;
	 double noise_gain_voiced;
	 double noise_low_freq_limit;
	 double compcoeff;
	 double concatenation_cost;
	 double target_cost;
	 double pulse_error_bias;
	 double pulse_tilt_decrease_coeff;
	 char *pulse_filename;
	 char *pulselibrary_filename;
	 char *synlist_filename;
	 char *dnnpath;
	 char **synlist;
	 gsl_vector *paramweights;
	 gsl_vector *dnn_weight_dims;
} PARAM;


int Check_command_line(int argc);
struct config_t *Read_config(char *filename);
int Assign_config_parameters(const char *filename, struct config_t *conf, PARAM *params, int conf_type);
int ReadSynthesisList(PARAM *params);
int Read_DNN_weights(PARAM *params, gsl_matrix **DNN_W);
int Read_input_minmax(PARAM *params, gsl_vector **input_minmax);
int Read_pulse_library(PARAM *params,gsl_matrix **pulses,gsl_matrix **pulses_rs,gsl_matrix **plsf,gsl_matrix **ptilt,gsl_matrix **pharm,
		gsl_matrix **phnr,gsl_matrix **pwaveform, gsl_matrix **pca_pc,gsl_matrix **pca_w_lib,gsl_vector **stoch_env,gsl_vector **stoch_sp,gsl_vector **pgain, gsl_vector **ph1h2,
		gsl_vector **pnaq, gsl_vector **pca_mean, gsl_vector **pulse_lengths);
gsl_vector *ReadPulseFile(PARAM *params);
int Initialize_params(PARAM *params,int synfilenumber);
int Compatibility_check(PARAM *params);
void Allocate_params(gsl_vector **excitation_voiced, gsl_vector **excitation_unvoiced, gsl_vector **resynthesis_pulse_index,
		gsl_vector **gain_new, gsl_matrix **glflowsp_new, gsl_matrix **hnr_new, PARAM *params);
int Read_synthesis_parameters(gsl_vector **gain, gsl_vector **fundf, gsl_matrix **LSF, gsl_matrix **LSF2, gsl_matrix **glflowsp,
		gsl_matrix **hnr, gsl_matrix **harmonics, gsl_matrix **waveform, gsl_vector **h1h2, gsl_vector **naq, gsl_matrix **pca_w, PARAM *params);
int Pulse_clustering(gsl_vector **pulse_clus_id, gsl_matrix **pulse_clusters, PARAM *params);
void Merge_voiced_unvoiced_spectra(gsl_matrix *LSF, gsl_matrix *LSF2, gsl_vector *fundf, PARAM *params);
void Integrate_LSFs(gsl_matrix **LSF, PARAM *params);
void Noise_robust_speech1(PARAM *params);
void Noise_robust_speech2(gsl_vector *gain, gsl_matrix *harmonics, PARAM *params);
void Postfilter(gsl_matrix *LSF, PARAM *params);
void Highpassfilter_fft(gsl_vector *signal);
void Highpassfilter_fir(gsl_vector *signal, gsl_vector *coeffs);
void Interpolate_matrix(gsl_matrix *matrix, gsl_matrix *imatrix);


void CreateExcitation(PARAM *params,gsl_vector *excitation_voiced,gsl_vector *excitation_unvoiced,gsl_vector *fundf,gsl_vector *gain,
	      gsl_matrix *lsf,gsl_matrix *glflowsp,gsl_matrix *hnr,gsl_matrix *harmonics,gsl_matrix *waveform,gsl_vector *h1h2,
	      gsl_vector *naq,gsl_vector *original_pulse,gsl_matrix *pulses,gsl_matrix *pulses_rs,gsl_vector *pulse_lengths,
	      gsl_vector *pgain,gsl_matrix *plsf,gsl_matrix *ptilt,gsl_matrix *phnr,gsl_matrix *pharm,gsl_matrix *pwaveform,
	      gsl_vector *ph1h2, gsl_vector *pnaq,gsl_vector *resynthesis_pulse_index,gsl_vector *pulse_clus_id,
	      gsl_matrix *pulse_clusters,gsl_vector *oldgain,gsl_matrix *glflowsp_new,gsl_matrix *hnr_new,
	      gsl_vector *pca_mean, gsl_matrix *pca_pc, gsl_matrix *pca_w, gsl_matrix *pca_w_lib, gsl_vector *stoch_env, gsl_vector *stoch_sp,
	      gsl_matrix **DNN_W, gsl_vector **input_minmax, gsl_vector *dnnpulseindices, gsl_vector *dnnpulses);

void Generate_DNN_pulse(gsl_vector *pulse, gsl_vector *fundf, gsl_vector *gain, gsl_vector *naq, gsl_vector *h1h2,
	gsl_matrix *hnr, gsl_matrix *glflowsp, gsl_matrix *lsf, gsl_matrix **DNN_W, gsl_vector **input_minmax, int index, PARAM *params);

void Generate_DNN_pulse_PCA(gsl_vector *pulse, gsl_vector *fundf, gsl_vector *gain, gsl_vector *naq, gsl_vector *h1h2,
	gsl_matrix *hnr, gsl_matrix *glflowsp, gsl_matrix *lsf, gsl_matrix **DNN_W, gsl_vector **input_minmax, gsl_matrix *pca_pc, int index, PARAM *params);

void Spectral_match(gsl_vector *signal, gsl_matrix *flow, gsl_matrix *flow_new, PARAM *params);

void Filter_excitation(gsl_vector *excitation, gsl_matrix *LSF, PARAM *params);

void Print_elapsed_time(PARAM *params);
void Print_synthesis_settings_start(PARAM *params);
void Print_synthesis_settings_middle(PARAM *params);
void Print_synthesis_settings_end(PARAM *params, double time1);

/* Gain modification */
void Evaluate_new_gain(gsl_vector *signal, gsl_vector *gain_new, gsl_vector *gain, gsl_vector *fundf, PARAM *params);
void Gain_eval(gsl_vector *signal, gsl_vector *gain, gsl_vector *fundf, PARAM *params);
void uvGain(gsl_vector *frame, gsl_vector *gain, int index, int windowing, int unvoiced_frame_length);
void Modify_gain(gsl_vector *gain_new, gsl_vector *gain);

void Scale_signal(gsl_vector *signal, int mode);
int Save_signal_to_file(gsl_vector *signal, PARAM *params, char *alternative_filename_ending);

void Reestimate_hnr(gsl_vector *excitation_voiced, gsl_vector *excitation_unvoiced, gsl_matrix *hnr, gsl_vector *fundf, PARAM *params);
void HNR_eval_vu(gsl_vector *signal_v,gsl_vector *signal_uv,gsl_matrix *hnr,gsl_vector *fundf,PARAM *params);
void Gain_normalization(gsl_vector *signal,gsl_vector *gain,int frame_length,int shift,double speed,double pitch,double gain_threshold);
void Evaluate_target_error(gsl_matrix *lsf, gsl_matrix *glflowsp, gsl_matrix *harmonics, gsl_matrix *hnr_i, gsl_matrix *waveform, gsl_vector *h1h2, gsl_vector *naq, gsl_vector *gain, gsl_vector *fundf,
		gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *pharm, gsl_matrix *phnr, gsl_matrix *pwaveform, gsl_vector *ph1h2, gsl_vector *pnaq, gsl_vector *pgain, gsl_matrix *pca_w, gsl_matrix *pca_w_lib, gsl_vector *pulse_lengths,
		int frame_index, gsl_vector *inds, gsl_vector *eparam, gsl_vector *pulse_clus_id, gsl_matrix *pulse_clusters, PARAM *params);
int Find_best_pulse_combination(gsl_vector *eparam, gsl_vector *eparam_next, gsl_vector *inds, gsl_vector *inds_next, gsl_vector *prevpulse, gsl_matrix *pulses_rs, int n_pulsecandidates, double target_cost, double concatenation_cost, int useprevpulse, int dist_to_unvoiced, double perror_bias);
void Evaluate_concatenation_error_viterbi(gsl_vector *target, gsl_vector* target_next, gsl_vector *inds, gsl_vector *inds_next,
					 gsl_matrix *pulses_rs, gsl_matrix *v_scores, gsl_matrix *v_best, int index, int dist_to_unvoiced, double additional_cost, PARAM *params);

void Phase_manipulation(gsl_vector *pulse, gsl_matrix *hnr, gsl_matrix *harmonics, int indec, PARAM *params);
void LSF_fix_vector(gsl_vector *lsf);
void LSF_fix_matrix(gsl_matrix *lsf);
void FIR_filter(gsl_vector *signal, gsl_vector *coeffs);
void IIR_filter(gsl_vector *signal, gsl_vector *coeffs);
void FIR_filter_file(gsl_vector *signal, char *filename, int n);
void Smooth_matrix(gsl_matrix *matrix,int len);
void LipRadiation(gsl_vector *excitation_voiced);
void Interpolate(gsl_vector *vector, gsl_vector *i_vector);
void Interpolate_lin(gsl_vector *vector, gsl_vector *i_vector);
void Interpolate_fract(gsl_vector *vector, gsl_vector *i_vector, double flen);
void Chebyshev(gsl_vector *T);
void Convert_to_LSF(gsl_matrix *LSF, gsl_vector *a, int index);
void Smooth(gsl_vector *vector, int length);
void MA(gsl_vector *vector, int length);
void lsf2poly(gsl_matrix *lsf_matrix, gsl_vector *poly, int index, int HMM);
void LSF_Postfilter(gsl_matrix *lsf, double alpha);
void LSF_Postfilter_mod(gsl_matrix *lsf, double alpha, gsl_vector *modulation);
void LPC_Postfilter(gsl_matrix *LSF, double gamma, int frame_length);
void Spectral_tilt(gsl_vector *excitation_voiced,gsl_vector *h1h2);
gsl_vector *Conv(gsl_vector *conv1, gsl_vector *conv2);
void Smooth_interp_lsf(gsl_matrix *LSF_i, gsl_matrix *LSF, int signal_len, int use_hmm, int lsf_smooth_len);
void Modify_h1h2(gsl_vector *h1h2);
gsl_vector *ReadFileDouble(char *name);
gsl_vector *ReadFileFloat(char *name);
void alphas2sigmas(double *alp, double *sigm, double lam, int dim);
double *NFArray(int size);
gsl_vector *WLPC(gsl_vector *frame, int p, double lambda);
void AllPassDelay(gsl_vector *signal, double lambda);
void Filter(gsl_vector *frame, gsl_vector *glottal, gsl_vector *a, int LPC_p);
void MA_voiced(gsl_vector *vector, gsl_vector *fundf, int len);
void MA_unvoiced(gsl_vector *vector, gsl_vector *fundf, int len);
void WFilter(gsl_vector *signal, gsl_vector *A, gsl_vector *B, double lambda);
void Upper_lower_envelope(gsl_vector *frame, gsl_matrix *hnr, double f0, int index, PARAM *params);
void MedFilt3(gsl_vector *frame);
void Convert_ERB2Hz(gsl_vector *vector_erb, gsl_vector *vector, PARAM *params);
void Convert_Hz2ERB(gsl_vector *vector, gsl_vector *vector_erb ,PARAM *params);
gsl_vector *Create_pulse_train(gsl_vector *pulse, gsl_vector *original_pulse, gsl_matrix *hnr, gsl_vector *fundf, gsl_matrix *harmonics, int frame_index, int sample_index, PARAM *params);
gsl_vector *Create_pulse_train_diff(gsl_vector *pulse, gsl_vector *original_pulse, gsl_matrix *hnr, gsl_vector *fundf, gsl_matrix *harmonics, int frame_index, int sample_index, PARAM *params);
void Analyse_pulse_train_spectrum(gsl_vector *pulse_train, gsl_matrix *glflowsp_new, int N, int sample_index, int frame_index, PARAM *params);
void MedFilt5_matrix(gsl_matrix *matrix);
void MedFilt5(gsl_vector *frame);
void HNR_compensation(gsl_matrix *hnr, gsl_matrix *hnr_new, PARAM *params);
void FillHNRValues(gsl_matrix *hnr, int frame_index_old, int frame_index);
void Differentiate(gsl_vector *signal, double leak);
void Differentiate_noleak(gsl_vector *signal);
void FillUnvoicedSyntheticSourceSpectrum(gsl_matrix *glflowsp, gsl_matrix *glflowsp_new, int frame_index, int frame_index_old);
void FillUnvoicedSyntheticSourceSpectrum_FLAT(gsl_matrix *glflowsp, gsl_matrix *glflowsp_new, int frame_index, int frame_index_old);
int EvalFileLength(const char *name, PARAM *params);
gsl_matrix *Construct_C(int size);
double Mean(gsl_vector *vector);
double NonZeroMean(gsl_vector *vector);
void Unwrap(gsl_vector *phase);
void Unwrap2(gsl_vector *phase, double tol);
void Wrap(gsl_vector *phase);
void Hp_filt_below_f0(gsl_vector *signal, gsl_vector *fundf, PARAM *params);
void Modify_harmonics(gsl_matrix *harmonics, double scale);
void Compression(gsl_vector *signal, double k);
void ModPowerSpectrum(gsl_vector *s, gsl_vector *formants, int n, double gamma);
void GetFormants(gsl_vector *s, gsl_vector *formants, int *n);
void Eval_STE_weight(gsl_vector *frame, gsl_vector *ste, int M, int lag);
void Eval_GCI_weight(gsl_vector *glottsig, gsl_vector *weight, gsl_vector *fundf, int FS, int index);
void LSF_stabilize(gsl_vector *a);
void Convert_vector_to_LSF(gsl_vector *a, gsl_vector *LSF);
void lsf_vector2poly(gsl_vector *lsf_vector, gsl_vector *poly);
gsl_vector *Find_GCI(gsl_vector *frame_orig, gsl_vector *fundf, int FS, int index);
void SWLP(gsl_vector *frame, gsl_vector *a, int M, int lag, int weighting, gsl_vector *fundf, int FS, int index, gsl_vector *glottsig, int stabilized);
void LSF_MOD(gsl_matrix *lsf);
void Evaluate_matrix_mean_std(gsl_vector *std, gsl_matrix *data);
void Evaluate_vector_std(gsl_vector *std, gsl_vector *data);
void Integrate_matrix(gsl_matrix *matrix);
void Integrate(gsl_vector *vector, double leak);
void Free_pulselib_variables(gsl_matrix *pulses, gsl_matrix *pulses_rs, gsl_matrix *pwaveform, gsl_vector *pulse_lengths,
		gsl_matrix *plsf, gsl_matrix *ptilt, gsl_matrix *pharm, gsl_matrix *phnr, gsl_vector *pgain, gsl_vector *ph1h2,
		gsl_vector *pnaq, gsl_vector *pca_mean, gsl_matrix *pca_pc, gsl_matrix *pca_w_lib, gsl_vector *stoch_env, gsl_vector *stoch_sp, PARAM *params);
void Free_variables(gsl_vector *original_pulse, gsl_vector *excitation_voiced, gsl_vector *excitation_unvoiced, gsl_vector *fundf, gsl_vector *gain,
		gsl_vector *gain_new, gsl_matrix *LSF, gsl_matrix *LSF2, gsl_matrix *LSF_interp, gsl_matrix *glflowsp, gsl_matrix *glflowsp_new,
		gsl_matrix *hnr, gsl_matrix *hnr_new, gsl_matrix *harmonics, gsl_matrix *waveform, gsl_vector *h1h2, gsl_vector *naq,
		gsl_vector *resynthesis_pulse_index, gsl_vector *pulse_clus_id, gsl_matrix *pulse_clusters, gsl_matrix *pca_w, 
		gsl_matrix **DNN_W, gsl_vector **input_minmax, gsl_vector *dnnpulseindices, gsl_vector *dnnpulses, PARAM *params);
void Save_excitation_to_wav(gsl_vector *excitation, PARAM *params);
void Noise_reduction(gsl_vector *gain, PARAM *params);
void Normalize_pulse_library_var(gsl_matrix *plsf,gsl_matrix *ptilt,gsl_matrix *pharm,gsl_matrix *phnr,gsl_matrix *pwaveform,
		gsl_vector *pgain,gsl_vector *ph1h2,gsl_vector *pnaq,gsl_matrix *lsf,gsl_matrix *tilt,gsl_matrix *harm,gsl_matrix *hnr,
		gsl_matrix *waveform,gsl_vector *gain,gsl_vector *h1h2,gsl_vector *naq, gsl_vector *fundf,PARAM *params);
void Adapt_synthesis_parameters_var(gsl_matrix *plsf,gsl_matrix *ptilt,gsl_matrix *pharm,gsl_matrix *phnr,gsl_matrix *pwaveform,
		gsl_vector *pgain,gsl_vector *ph1h2,gsl_vector *pnaq,gsl_matrix *lsf,gsl_matrix *tilt,gsl_matrix *harm,gsl_matrix *hnr,
		gsl_matrix *waveform,gsl_vector *gain,gsl_vector *h1h2,gsl_vector *naq,gsl_vector *fundf,gsl_vector *pulse_lengths,PARAM *params);
void Select_LSFs_from_pulse_library(gsl_matrix *lsf, gsl_matrix *plsf, gsl_vector *fundf, int win, PARAM *params);
void Average_pulses(gsl_vector *pulse, gsl_vector *fundf, gsl_vector *resynthesis_pulse_index, gsl_matrix *pulses, gsl_vector *pulse_lengths, int frame_index, PARAM *params);
void Shift_right_and_add(gsl_vector *vector, double value);
void Convert_logF0_to_lin(gsl_vector *fundf, PARAM *params);
gsl_vector *Truncate_pulse(gsl_vector *pulse, gsl_vector *fundf, int sample_index, int frame_index, PARAM *params);

void Fill_pulse_indices(gsl_vector *ind);

/* Not in use */
void Create_highband_excitation(gsl_vector *excitation_highband, gsl_matrix *erbgain, PARAM *params);

/* Test functions */
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
void VPrint11(gsl_vector *vector);
void VPrint12(gsl_vector *vector);
void MPrint1(gsl_matrix *matrix);
void MPrint2(gsl_matrix *matrix);
void MPrint3(gsl_matrix *matrix);
void MPrint4(gsl_matrix *matrix);
void MPrint5(gsl_matrix *matrix);
void MPrint6(gsl_matrix *matrix);
void APrint1(double *array, int size);
void APrint2(double *array, int size);
void APrint3(double *array, int size);
void APrint4(double *array, int size);
void TimePrint(double start_time);
void pause(int print);
int Find_matrix_NaN_Inf(gsl_matrix *m);
int Find_vector_NaN_Inf(gsl_vector *v);

#endif

