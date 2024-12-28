#include "audio.h"

#define PI 3.1415926

#define UINT24_MAX 16777215
#define INT24_MAX 8388607
#define INT24_MIN -8388608

bool is_little_endian_system;

union sample {
	float 		f32;
	double 		f64;

	int8_t 		s8;
	int16_t 	s16;
	int32_t 	s32;

	uint8_t 	u8;
	uint16_t 	u16;
	uint32_t 	u32;

	uint8_t 	i64_arr [8];
};


union uint32bit {
		uint8_t bytes[4];
		uint32_t uint32;
	};

union uint16bit {
	uint8_t bytes[2];
	uint16_t uint16;
};

/*-----------------------------------------------------------------------------*/
/* AUXILIARY FUNCTIONS */

int time_to_sample_number ( Audio* audio, double time ) {
	return (int)(time*audio->sample_rate);
}

bool is_little_endian () {
	union endian_test {
		uint16_t var1;
		uint8_t var2 [2];
	};

	union endian_test et;

	et.var1 = 0xffee;

	if ( et.var2[0] == 0xee ) {
		return true;
	}
	else {
		return false;
	}
}

void switch_endianness ( union sample* s, unsigned int byte_quantity ) {
	union sample s_temp = *s;

	for (unsigned int i=0; i<byte_quantity; i++) {
		s->i64_arr[i] = s_temp.i64_arr[byte_quantity-i-1];
	}
}

void serialize_transcoded_audio ( char* byte_sequence, AudioTranscoded* audio ) {
	unsigned int it = 0;

	for (int i=0; i<audio->sample_quantity_per_channel; i++) {
		for (int j=0; j<audio->channels; j++) {
			if (audio->format == S8 || audio->format == U8) {
				byte_sequence[it] = (char)audio->raw_audio[j][i].i64_arr[0];
				it++;
			}
			else if (audio->format == S16LE || audio->format == S16BE || audio->format == U16LE || audio->format == U16BE) {
				for (int k=0; k<2; k++) {
					byte_sequence[it] = (char)audio->raw_audio[j][i].i64_arr[k];
					it++;
				}
			}
			else if (audio->format == S24LE || audio->format == S24BE || audio->format == U24LE || audio->format == U24BE) {
				for (int k=0; k<3; k++) {
					byte_sequence[it] = (char)audio->raw_audio[j][i].i64_arr[k];
					it++;
				}
			}
			else if (audio->format == S32LE || audio->format == S32BE || audio->format == U32LE || audio->format == U32BE || audio->format == F32LE || audio->format == F32BE) {
				for (int k=0; k<4; k++) {
					byte_sequence[it] = (char)audio->raw_audio[j][i].i64_arr[k];
					it++;
				}
			}
			else if (audio->format == F64LE || audio->format == F64BE) {
				for (int k=0; k<8; k++) {
					byte_sequence[it] = (char)audio->raw_audio[j][i].i64_arr[k];
					it++;
				}
			}
		}
	}
}

bool audios_are_compatible ( Audio* audio1, Audio* audio2 ) {
	if ( audio1->sample_rate != audio2->sample_rate ) {
		return false;
	}
	if ( audio1->channels != audio2->channels ) {
		return false;
	}
	return true;
}

/*-----------------------------------------------------------------------------*/


void init_audio ( Audio* audio, unsigned int sample_rate, double duration, unsigned int channels ) {
	if ( duration < 0.0 ) {
		return;
	}

	audio->sample_quantity_per_channel = (unsigned int)(sample_rate * duration);
	audio->sample_quantity_total = (unsigned int)(sample_rate * duration * channels);
	audio->sample_rate = sample_rate;
	audio->duration = duration;
	audio->channels = channels;
	audio->current_place = 0.0;
	
	audio->raw_audio = (double**) malloc( sizeof(double*) * channels );
	for (int i=0; i<channels; i++) {
		audio->raw_audio[i] = (double*) malloc( sizeof(double) * audio->sample_quantity_per_channel );
	}
	
	is_little_endian_system = is_little_endian();

	srand(time(NULL));
}

void init_transcoded_audio ( AudioTranscoded* audio, unsigned int sample_rate, unsigned int format, double duration, unsigned int channels ) {
	if ( duration < 0.0 ) {
		return;
	}
	if ( format != F32LE && format != F32BE && format != F64LE && format != F64BE && format != S8 && format != S16LE && format != S16BE && format != S32LE && format != S32BE && format != U8 && format != U16LE && format != U16BE && format != U32LE && format != U32BE ) {
		return;
	}

	if ( format == S8 || format == U8 ) {
		audio->bytes_per_sample = 1;
	}
	else if ( format == S16LE || format == S16BE || format == U16LE || format == U16BE ) {
		audio->bytes_per_sample = 2;
	}
	else if ( format == S24LE || format == S24BE || format == U24LE || format == U24BE ) {
		audio->bytes_per_sample = 3;
	}
	else if ( format == S32LE || format == S32BE || format == U32LE || format == U32BE || format == F32LE || format == F32BE ) {
		audio->bytes_per_sample = 4;
	}
	else if ( format == F64LE || format == F64BE ) {
		audio->bytes_per_sample = 8;
	}

	audio->sample_quantity_per_channel = (unsigned int)(sample_rate * duration);
	audio->sample_quantity_total = (unsigned int)(sample_rate * duration * channels);
	audio->sample_rate = sample_rate;
	audio->format = format;
	audio->duration = duration;
	audio->channels = channels;
	audio->current_place = 0.0;
	
	audio->raw_audio = (union sample**) malloc( sizeof(union sample*) * channels );
	for (int i=0; i<channels; i++) {
		audio->raw_audio[i] = (union sample*) malloc( sizeof(union sample) * audio->sample_quantity_per_channel );
	}
	
	is_little_endian_system = is_little_endian();
}

void free_audio ( Audio* audio ) {
	for (int i=0; i<audio->channels; i++) {
		free( audio->raw_audio[i] );
	}
	free( audio->raw_audio );
}

void free_transcoded_audio ( AudioTranscoded* audio ) {
	for (int i=0; i<audio->channels; i++) {
		free( audio->raw_audio[i] );
	}
	free( audio->raw_audio );
}


/*-----------------------------------------------------------------------------*/
/* AUDIO GENERATORS */

void white_noise ( Audio* audio, double max_volume, double start, double duration ) {
	if ( max_volume < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = (((double)rand() / (double)RAND_MAX) - 0.5) * 2.0 * max_volume * DBL_MAX;
		}
	}
	
	audio->current_place = start+duration;
}

void sine_wave ( Audio* audio, double frequency, double max_volume, double start, double duration ) {
	if ( max_volume < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	double n;

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	double it_step = 1.0 / audio->sample_rate;
	double it = 0.0;

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			n = (double)(sin( 2.0*PI * (frequency/audio->sample_rate) * i ) * max_volume * DBL_MAX);
			audio->raw_audio[j][i] = n;
			it += it_step;
		}
	}
	
	audio->current_place = start+duration;
}

void cosine_wave ( Audio* audio, double frequency, double max_volume, double start, double duration ) {
	if ( max_volume < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	double n;

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	double it_step = 1.0 / audio->sample_rate;
	double it = 0.0;

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			n = (double)(cos( 2.0*PI * (frequency/audio->sample_rate) * i ) * max_volume * DBL_MAX);
			audio->raw_audio[j][i] = n;
			it += it_step;
		}
	}

	audio->current_place = start+duration;
}

void square_wave ( Audio* audio, double frequency, double max_volume, double start, double duration ) {
	if ( max_volume < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	double n;

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;;
	}

	unsigned int lambda = audio->sample_rate / frequency;
	
	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		if ( i%lambda < lambda/2 ) {
			n = (double)(max_volume * DBL_MAX);
		}
		else {
			n = (double)(-max_volume * DBL_MAX);
		}
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = n;
		}
	}
	
	audio->current_place = start+duration;
}

void sawtooth_wave ( Audio* audio, double frequency, double max_volume, double start, double duration ) {
	if ( max_volume < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	double level;
	double n;

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int lambda = audio->sample_rate / frequency;
	
	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		level = ((2.0*max_volume)/lambda)*(i%(lambda/2));

		if ( i%lambda >= lambda/2 ) {
			level = -(max_volume-level);
		}
		
		n = (double)(level*DBL_MAX);
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = n;
		}
	}
	
	audio->current_place = start+duration;
}


void triangular_wave ( Audio* audio, double frequency, double max_volume, double start, double duration ) {
	if ( max_volume < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	double level;
	double n;

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int lambda = audio->sample_rate / frequency;

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		if ( i%lambda < lambda/4 ) {
			level = max_volume * ((double)(i%(lambda/4))/(double)(lambda/4));
		}
		else if ( i%lambda >= lambda/4 && i%lambda < lambda/2 ) {
			level = -max_volume * ((double)(i%(lambda/4))/(double)(lambda/4)) + max_volume;
		}
		else if ( i%lambda >= lambda/2 && i%lambda < (lambda*3)/4 ) {
			level = -max_volume * ((double)(i%(lambda/4))/(double)(lambda/4));
		}
		else {
			level = max_volume * ((double)(i%(lambda/4))/(double)(lambda/4)) - max_volume;
		}

		n = (double)(level*DBL_MAX);

		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = n;
		}
	}

	audio->current_place = start+duration;
}


void constant_level ( Audio* audio, double level, double start, double duration ) {
	if ( level < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = (double)(level*DBL_MAX);
		}
	}
	
	audio->current_place = start+duration;
}

void silence ( Audio* audio, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	constant_level(audio, 0.0, start, start+duration);
	
	audio->current_place = start+duration;
}



/*-----------------------------------------------------------------------------*/
/* AUDIO PROCESSORS */

void am_modulate ( Audio* audio, Audio* carrier ) {
	if ( !audios_are_compatible(audio, carrier) ) {
		return;
	}

	for (int i=0; i<audio->sample_quantity_per_channel; i++) {
		audio->raw_audio[0][i] = (double)(audio->raw_audio[0][i] * ((double)carrier->raw_audio[0][i]/DBL_MAX));
	}
}

void amplify ( Audio* audio, double amp_factor, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = (double)(audio->raw_audio[j][i]*amp_factor);
		}
	}
}

void shift_vertically ( Audio* audio, double shift, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = audio->raw_audio[j][i] + (double)(shift*DBL_MAX);
		}
	}
}

void shift_horizontally ( Audio* audio, double shift, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample;
	unsigned int end_sample;
	
	int copy_start;
	int copy_end;
	int copy_offset;
	
	if (shift < 0.0) {
		shift = -shift;
		copy_start = time_to_sample_number(audio, start-shift);
		copy_end = time_to_sample_number(audio, start+duration-shift);
		copy_offset = time_to_sample_number(audio, shift);
		
		start_sample = time_to_sample_number(audio, start+duration-shift);
		end_sample = time_to_sample_number(audio, start+duration);
		
		for ( int i=copy_start; i<copy_end; i++ ) {
			if (i < 0) {
				continue;
			}
			for (int j=0; j<audio->channels; j++) {
				audio->raw_audio[j][i] = audio->raw_audio[j][i+copy_offset];
			}
		}
		for ( int i=start_sample; i<end_sample; i++ ) {
			for ( int j=0; j<audio->channels; j++ ) {
				audio->raw_audio[j][i] = 0.0;
			}
		}
	}
	else {
		copy_start = time_to_sample_number(audio, start+shift+duration);
		copy_end = time_to_sample_number(audio, start+shift);
		copy_offset = time_to_sample_number(audio, shift);
		
		start_sample = time_to_sample_number(audio, start);
		end_sample = time_to_sample_number(audio, start+shift);
	
		for ( unsigned int i=copy_start; i>=copy_end; i-- ) {
			if (i >= audio->sample_quantity_per_channel) {
				continue;
			}
			for (int j=0; j<audio->channels; j++) {
				audio->raw_audio[j][i] = audio->raw_audio[j][i-copy_offset];
			}
		}
		for ( unsigned int i=start_sample; i<end_sample; i++ ) {
			for (int j=0; j<audio->channels; j++) {
				audio->raw_audio[j][i] = 0.0;
			}
		}
	}
}

void limit_volume ( Audio* audio, double max_volume, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	amplify(audio, max_volume, start, start+duration);
}

void resample ( Audio* dest_audio, Audio* src_audio, unsigned int sample_rate ) {
	init_audio(dest_audio, sample_rate, src_audio->duration, src_audio->channels);
	
	double resample_step = (double)src_audio->sample_rate / (double)sample_rate;
	double resample_it = 0.0;

	for (int i=0; i<dest_audio->sample_quantity_per_channel; i++) {
		for (int j=0; j<dest_audio->channels; j++) {
			dest_audio->raw_audio[j][i] = src_audio->raw_audio[j][(unsigned int)resample_it];
		}
		resample_it += resample_step;
	}
}

void insert_audio ( Audio* audio, Audio* paste_audio, double start_audio, double start_paste, double duration_paste ) {
	if ( !audios_are_compatible(audio, paste_audio) ) {
		return;
	}
	if ( start_audio < 0.0 || start_paste < 0.0 || duration_paste < 0.0 ) {
		return;
	}

	unsigned int start_audio_sample = time_to_sample_number(audio, start_audio);
	unsigned int start_paste_sample = time_to_sample_number(paste_audio, start_paste);
	unsigned int end_paste_sample = time_to_sample_number(paste_audio, start_paste+duration_paste);

	unsigned int paste_it = start_paste_sample;

	for (unsigned int i=start_audio_sample; i<start_audio_sample+(end_paste_sample-start_paste_sample); i++) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = paste_audio->raw_audio[j][paste_it++];
		}
	}
	
	audio->current_place = start_audio+duration_paste;
}

void append_audio ( Audio* dest_audio, Audio* src_audio, double start_src_audio, double duration_src_audio ) {
	if ( !audios_are_compatible(dest_audio, src_audio) ) {
		return;
	}
	if ( start_src_audio < 0.0 || duration_src_audio < 0.0 ) {
		return;
	}

	unsigned int start = (unsigned int)(start_src_audio * src_audio->sample_rate);
	unsigned int end = start_src_audio + (unsigned int)(duration_src_audio*src_audio->sample_rate);

	for (int i=start; i<end; i++) {
		for (int j=0; j<dest_audio->channels; j++) {
			dest_audio->raw_audio[j][i+start] = src_audio->raw_audio[j][i];
		}
	}

	dest_audio->current_place = dest_audio->current_place + duration_src_audio;
}


void append_sound (Audio* audio, Sound* sound, double frequency, void(*generator_function)(Audio*, double, double, double, double)) {
	if ( frequency < 0.0 ) {
		return;
	}

	double sound_duration = sound->delay + sound->attack + sound->hold + sound->decay;

	Audio sound_audio;
	init_audio(&sound_audio, audio->sample_rate, sound_duration, audio->channels);

	generator_function(&sound_audio, frequency, sound->max_level, sound->delay, sound_duration-sound->delay);
	fade_in(&sound_audio, sound->attack, 1.0, sound->delay);
	fade_out(&sound_audio, sound->decay, MIN_LEVEL, sound->delay+sound->attack+sound->hold);

	append_audio(audio, &sound_audio, BEGIN, sound_audio.duration);
}


void fade_in ( Audio* audio, double duration, double to_level, double start ) {
	if ( to_level < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int duration_in_samples = (unsigned int)(duration * audio->sample_rate);
	unsigned int start_sample = time_to_sample_number(audio, start);

	double amp_factor = to_level / duration_in_samples;
	double fade_it = 0.0;

	for ( unsigned int i=start_sample; i<start_sample+duration_in_samples; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = (double)(audio->raw_audio[j][i] * fade_it);
		}
		fade_it += amp_factor;
	}
	
	audio->current_place = start+duration;
}

void fade_out( Audio* audio, double duration, double to_level, double start ) {
	if ( to_level < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int duration_in_samples = (unsigned int)(duration * audio->sample_rate);
	unsigned int start_sample = time_to_sample_number(audio, start);

	double amp_factor = (1.0-to_level) / duration_in_samples;
	double fade_it = 1.0;

	for ( unsigned int i=start_sample; i<start_sample+duration_in_samples; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = (double)(audio->raw_audio[j][i] * fade_it);
		}
		fade_it -= amp_factor;
	}
	
	audio->current_place = start+duration;
}


void make_positive ( Audio* audio, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			if (audio->raw_audio[j][i] < 0) {
				audio->raw_audio[j][i] = -audio->raw_audio[j][i];
			}
		}
	}
}

void make_negative ( Audio* audio, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			if (audio->raw_audio[j][i] > 0) {
				audio->raw_audio[j][i] = -audio->raw_audio[j][i];
			}
		}
	}
}

void invert ( Audio* audio, double start, double duration ) {
	if ( duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] = -audio->raw_audio[j][i];
		}
	}
}


void compress ( Audio* audio, double max_level, double start, double duration ) {
	if ( max_level < 0.0 || duration < 0.0 || start+duration >= audio->duration ) {
		return;
	}

	if (start == END) {
		start = audio->duration;
	}
	if (start == CURRENT) {
		start = audio->current_place;
	}

	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		for (int j=0; j<audio->channels; j++) {
			if ( abs(audio->raw_audio[j][i]) > (double)(max_level*DBL_MAX) ) {
				if (audio->raw_audio[j][i] < 0.0) {
					audio->raw_audio[j][i] = (double)(-max_level*DBL_MAX);
				}
				else {
					audio->raw_audio[j][i] = (double)(max_level*DBL_MAX);
				}
			}
		}
	}
}


void add_audio ( Audio* audio, Audio* audio_to_add ) {
	if ( !audios_are_compatible(audio, audio_to_add) ) {
		return;
	}

	for (unsigned int i=0; i<audio->sample_quantity_per_channel; i++) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] += audio_to_add->raw_audio[j][i];
		}
	}
}

void subtract_audio ( Audio* audio, Audio* audio_to_subtract ) {
	if ( !audios_are_compatible(audio, audio_to_subtract) ) {
		return;
	}

	for (unsigned int i=0; i<audio->sample_quantity_per_channel; i++) {
		for (int j=0; j<audio->channels; j++) {
			audio->raw_audio[j][i] -= audio_to_subtract->raw_audio[j][i];
		}
	}
}


/*-----------------------------------------------------------------------------*/
/* INPUT / OUTPUT */

void write_audio_to_file ( Audio* audio, unsigned int format, unsigned int sample_rate, char* file_path ) {
	unsigned int bytes_per_sample;

	if ( format == S8 || format == U8 ) {
		bytes_per_sample = 1;
	}
	else if ( format == S16LE || format == S16BE || format == U16LE || format == U16BE ) {
		bytes_per_sample = 2;
	}
	else if ( format == S24LE || format == S24BE || format == U24LE || format == U24BE ) {
		bytes_per_sample = 3;
	}
	else if ( format == S32LE || format == S32BE || format == U32LE || format == U32BE || format == F32LE || format == F32BE ) {
		bytes_per_sample = 4;
	}
	else if ( format == F64LE || format == F64BE ) {
		bytes_per_sample = 8;
	}
	else {
		return;
	}

	unsigned int bytes_quantity = bytes_per_sample * audio->sample_quantity_total;
	char* byte_sequence = (char*) malloc(bytes_quantity);

	FILE* fptr;
	fptr = fopen(file_path, "wb");

	Audio audio_resampled;
	resample(&audio_resampled, audio, sample_rate);

	AudioTranscoded audio_transcoded;
	transcode_raw_audio(&audio_transcoded, &audio_resampled, format);


	serialize_transcoded_audio(byte_sequence, &audio_transcoded);

	for (int i=0; i<bytes_quantity; i++) {
		fputc(byte_sequence[i], fptr);
	}

	free(byte_sequence);
	
	fclose(fptr);
}

void create_wav_file (Audio* audio, unsigned int format, char* file_path) {
	AudioTranscoded audio_transcoded;
	transcode_raw_audio( &audio_transcoded, audio, format );


	char wav_header [44];

	wav_header[0] = 'R';
	wav_header[1] = 'I';
	wav_header[2] = 'F';
	wav_header[3] = 'F';

	union uint32bit total_file_size;

	total_file_size.uint32 = audio_transcoded.sample_quantity_total * audio_transcoded.bytes_per_sample + 44;

	wav_header[4] = total_file_size.bytes[0];
	wav_header[5] = total_file_size.bytes[1];
	wav_header[6] = total_file_size.bytes[2];
	wav_header[7] = total_file_size.bytes[3];

	wav_header[8] = 'W';
	wav_header[9] = 'A';
	wav_header[10] = 'V';
	wav_header[11] = 'E';

	wav_header[12] = 'f';
	wav_header[13] = 'm';
	wav_header[14] = 't';
	wav_header[15] = ' ';

	union uint32bit length_of_format_data;
	length_of_format_data.uint32 = 16;

	wav_header[16] = length_of_format_data.bytes[0];
	wav_header[17] = length_of_format_data.bytes[1];
	wav_header[18] = length_of_format_data.bytes[2];
	wav_header[19] = length_of_format_data.bytes[3];

	union uint16bit format_type;
	format_type.uint16 = 1;

	wav_header[20] = format_type.bytes[0];
	wav_header[21] = format_type.bytes[1];

	union uint16bit number_of_channels;
	number_of_channels.uint16 = audio_transcoded.channels;

	wav_header[22] = number_of_channels.bytes[0];
	wav_header[23] = number_of_channels.bytes[1];

	union uint32bit sample_rate;
	sample_rate.uint32 = audio_transcoded.sample_rate;

	wav_header[24] = sample_rate.bytes[0];
	wav_header[25] = sample_rate.bytes[1];
	wav_header[26] = sample_rate.bytes[2];
	wav_header[27] = sample_rate.bytes[3];

	union uint32bit byte28to31;

	byte28to31.uint32 = audio_transcoded.sample_rate * audio_transcoded.bytes_per_sample * audio_transcoded.channels;

	wav_header[28] = byte28to31.bytes[0];
	wav_header[29] = byte28to31.bytes[1];
	wav_header[30] = byte28to31.bytes[2];
	wav_header[31] = byte28to31.bytes[3];

	union uint16bit byte32to33;

	byte32to33.uint16 = (uint16_t)(((double)audio_transcoded.sample_rate * (double)audio_transcoded.channels) / 8.1);

	wav_header[32] = byte32to33.bytes[0];
	wav_header[33] = byte32to33.bytes[1];

	union uint16bit bits_per_sample;

	bits_per_sample.uint16 = audio_transcoded.bytes_per_sample * 8;
	wav_header[34] = bits_per_sample.bytes[0];
	wav_header[35] = bits_per_sample.bytes[1];


	wav_header[36] = 'd';
	wav_header[37] = 'a';
	wav_header[38] = 't';
	wav_header[39] = 'a';

	union uint32bit data_section_size;
	data_section_size.uint32 = audio_transcoded.sample_quantity_total * audio_transcoded.bytes_per_sample;

	wav_header[40] = data_section_size.bytes[0];
	wav_header[41] = data_section_size.bytes[1];
	wav_header[42] = data_section_size.bytes[2];
	wav_header[43] = data_section_size.bytes[3];

	char* byte_sequence = (char*) malloc(data_section_size.uint32 + 44);

	for(int i=0; i<44; i++) {
		byte_sequence[i] = (char) wav_header[i];
	}

	FILE* fptr;
	fptr = fopen(file_path, "wb");

	serialize_transcoded_audio(&byte_sequence[44], &audio_transcoded);

	for (int i=0; i<data_section_size.uint32+44; i++) {
		fputc(byte_sequence[i], fptr);
	}

	free(byte_sequence);

	fclose(fptr);
}

void print_binary_audio ( Audio* audio ) {
	for (int i=0; i<audio->sample_quantity_per_channel; i++) {
		for (int j=0; j<audio->channels; j++) {
			putchar( (char)audio->raw_audio[j][i] );
		}
	}
}


void load_audio_from_raw_file ( Audio* audio, char* file_path, unsigned int sample_rate, unsigned int format, unsigned int channels ) {
	unsigned int bytes_per_sample = 0;
	
	if (format == S8 || format == U8) {
		bytes_per_sample = 1;
	}
	else if (format == S16LE || format == S16BE || format == U16LE || format == U16BE) {
		bytes_per_sample = 2;
	}
	else if ( format == S24LE || format == S24BE || format == U24LE || format == U24BE ) {
		bytes_per_sample = 3;
	}
	else if (format == S32LE || format == S32BE || format == U32LE || format == U32BE || format == F32LE || format == F32BE) {
		bytes_per_sample = 4;
	}
	else if (format == F64LE || format == F64BE) {
		bytes_per_sample = 8;
	}
	else {
		return;
	}


	FILE* audio_file = fopen(file_path, "r");

	fseek(audio_file, 0, SEEK_END);
	long file_size = ftell(audio_file);

	rewind(audio_file);
	

	double length = (double)file_size / (double)sample_rate / (double)channels / (double)bytes_per_sample;
	
	init_audio(audio, sample_rate, length, channels);

	char c = fgetc(audio_file);
	union sample s;
	
	int it = 0;

	while ( c != EOF ) {
		
		switch (format) {
			case F32LE:
				for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 4);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/FLT_MAX)*s.f32);
				break;
			case F32BE:
				for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 4);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/FLT_MAX)*s.f32);
				break;
			case F64LE:
				for (int i=0; i<8; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 8);
				}
				audio->raw_audio[it%channels][it/channels] = s.f64;
				break;
			case F64BE:
				for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 8);
				}
				audio->raw_audio[it%channels][it/channels] = s.f64;
				break;

			case S8:
				s.i64_arr[0] = (uint8_t)c;
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT8_MAX)*s.s8);
				break;
			case S16LE:
				for (int i=0; i<2; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 2);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT16_MAX)*s.s16);
				break;
			case S16BE:
				for (int i=0; i<2; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 2);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT16_MAX)*s.s16);

				break;
			case S24LE:
				for (int i=0; i<3; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 3);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT24_MAX)*s.s32);
				break;
			case S24BE:
				for (int i=0; i<3; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 3);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT24_MAX)*s.s32);

				break;
			case S32LE:
				for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 4);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT32_MAX)*s.s32);
				break;
			case S32BE:
			for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 4);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/INT32_MAX)*s.s32);
				break;


			// TODO: MIGHT NOT WORK !!!
			case U8:
				s.i64_arr[0] = (uint8_t)c;
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT8_MAX)*s.u8);
				break;
			case U16LE:
				for (int i=0; i<2; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 2);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT16_MAX)*s.u16);
				break;
			case U16BE:
				for (int i=0; i<2; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 2);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT16_MAX)*s.u16);
				break;
			case U24LE:
				for (int i=0; i<3; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 3);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT24_MAX)*s.u32);
				break;
			case U24BE:
				for (int i=0; i<3; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 3);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT24_MAX)*s.u32);
				break;
			case U32LE:
				for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 4);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT32_MAX)*s.u32);
				break;
			case U32BE:
				for (int i=0; i<4; i++) {
					s.i64_arr[i] = (uint8_t)c;
					c = fgetc(audio_file);
				}
				if ( !is_little_endian_system ) {
					switch_endianness(&s, 4);
				}
				audio->raw_audio[it%channels][it/channels] = (double)((DBL_MAX/UINT32_MAX)*s.u32);
				break;
			}
			
			it++;
	}

	fclose(audio_file);
}

void load_audio_from_wav_file ( Audio* audio, char* file_path ) {
	FILE* wav_file = fopen(file_path, "r");

	fseek(wav_file, 0L, SEEK_END);
	unsigned int file_size = ftell(wav_file);
	rewind(wav_file);

	char* file = malloc(file_size);

	for (unsigned int i=0; i<file_size; i++) {
		file[i] = fgetc(wav_file);
	}

	union uint32bit sample_rate;
	sample_rate.bytes[0] = file[24];
	sample_rate.bytes[1] = file[25];
	sample_rate.bytes[2] = file[26];
	sample_rate.bytes[3] = file[27];

	union uint32bit data_section_size;
	data_section_size.bytes[0] = file[40];
	data_section_size.bytes[1] = file[41];
	data_section_size.bytes[2] = file[42];
	data_section_size.bytes[3] = file[43];

	union uint16bit number_of_channels;
	number_of_channels.bytes[0] = file[22];
	number_of_channels.bytes[1] = file[23];

	union uint16bit bytes_per_sample;
	bytes_per_sample.bytes[0] = file[34];
	bytes_per_sample.bytes[1] = file[35];
	bytes_per_sample.uint16 /= 8;


	double duration = (double)data_section_size.uint32 / (double)number_of_channels.uint16 / (double)bytes_per_sample.uint16 / (double)sample_rate.uint32;

	init_audio(audio, sample_rate.uint32, duration, number_of_channels.uint16);

	union sample current_sample;
	double current_sample_double;

	unsigned int sample_it = 0;


	for (unsigned int i=0; i<data_section_size.uint32; i+=bytes_per_sample.uint16) {
		for (int j=0; j<bytes_per_sample.uint16; j++) {
			current_sample.i64_arr[j] = file[i+j+44];
		}

		switch (bytes_per_sample.uint16) {
			case 1:
				current_sample_double = (double)current_sample.s8 * (DBL_MAX/INT8_MAX);
				break;
			case 2:
				current_sample_double = (double)current_sample.s16 * (DBL_MAX/INT16_MAX);
				break;
			case 3:
				current_sample_double = (double)current_sample.s32 * (DBL_MAX/INT24_MAX);
				break;
			case 4:
				current_sample_double = (double)current_sample.s32 * (DBL_MAX/INT32_MAX);
				break;
		}

		audio->raw_audio[sample_it%audio->channels][sample_it/audio->channels] = current_sample_double;
		sample_it++;

		for (int j=0; j<8; j++) {
			current_sample.i64_arr[j] = 0;
		}
	}
};

/*-----------------------------------------------------------------------------*/
/* TRANSCODER */
void transcode_raw_audio ( AudioTranscoded* dest_audio, Audio* src_audio, unsigned int format ) {
	init_transcoded_audio(dest_audio, src_audio->sample_rate, format, src_audio->duration, src_audio->channels);

	for (int i=0; i<src_audio->sample_quantity_per_channel; i++) {
		for (int j=0; j<src_audio->channels; j++) {

			switch (format) {
				case F32LE:
					dest_audio->raw_audio[j][i].f32 = (float)(src_audio->raw_audio[j][i] / DBL_MAX);
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 4);
					}
					break;
				case F32BE:
					dest_audio->raw_audio[j][i].f32 = (float)(src_audio->raw_audio[j][i] / DBL_MAX);
					if ( is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 4);
					}
					break;
				case F64LE:
					dest_audio->raw_audio[j][i].f64 = (double)src_audio->raw_audio[j][i] / DBL_MAX;
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 8);
					}
					break;
				case F64BE:
					dest_audio->raw_audio[j][i].f64 = (double)src_audio->raw_audio[j][i] / DBL_MAX;
					if ( is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 8);
					}
					break;

				case S8:
					dest_audio->raw_audio[j][i].s8 = (int8_t)(src_audio->raw_audio[j][i] * (INT8_MAX/DBL_MAX));
					break;
				case S16LE:
					dest_audio->raw_audio[j][i].s16 = (int16_t)(src_audio->raw_audio[j][i] * (INT16_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 2);
					}
					break;
				case S16BE:
					dest_audio->raw_audio[j][i].s16 = (int16_t)(src_audio->raw_audio[j][i] * (INT16_MAX/DBL_MAX));
					if ( is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 2);
					}
					break;
				case S24LE:
					dest_audio->raw_audio[j][i].s32 = (int32_t)(src_audio->raw_audio[j][i] * (INT24_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 3);
					}
					break;
				case S24BE:
					dest_audio->raw_audio[j][i].s32 = (int32_t)(src_audio->raw_audio[j][i] * (INT24_MAX/DBL_MAX));
					if ( is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 3);
					}
					break;
				case S32LE:
					dest_audio->raw_audio[j][i].s32 = (int32_t)(src_audio->raw_audio[j][i] * (INT32_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 4);
					}
					break;
				case S32BE:
					dest_audio->raw_audio[j][i].s32 = (int32_t)(src_audio->raw_audio[j][i] * (INT32_MAX/DBL_MAX));
					if ( is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 4);
					}
					break;

				// TODO: MIGHT NOT WORK !!!
				case U8:
					dest_audio->raw_audio[j][i].u8 = (uint8_t)(src_audio->raw_audio[j][i] * ((UINT8_MAX)/DBL_MAX));
					break;
				case U16LE:
					dest_audio->raw_audio[j][i].u16 = (uint16_t)(src_audio->raw_audio[j][i] * (UINT16_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 2);
					}
					break;
				case U16BE:
					dest_audio->raw_audio[j][i].u16 = (uint16_t)(src_audio->raw_audio[j][i] * (UINT16_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 2);
					}
					break;
				case U24LE:
					dest_audio->raw_audio[j][i].u32 = (uint32_t)(src_audio->raw_audio[j][i] * (UINT24_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 3);
					}
					break;
				case U24BE:
					dest_audio->raw_audio[j][i].u32 = (uint32_t)(src_audio->raw_audio[j][i] * (UINT24_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 3);
					}
					break;
				case U32LE:
					dest_audio->raw_audio[j][i].u32 = (uint32_t)(src_audio->raw_audio[j][i] * (UINT32_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 4);
					}
					break;
				case U32BE:
					dest_audio->raw_audio[j][i].u32 = (uint32_t)(src_audio->raw_audio[j][i] * (UINT32_MAX/DBL_MAX));
					if ( !is_little_endian_system ) {
						switch_endianness(&dest_audio->raw_audio[j][i], 4);
					}
					break;
			}
		}
	}
}


/*-----------------------------------------------------------------------------*/
/* INFORMATION */

double get_highest_level ( Audio* audio, double start, double duration ) {
	double n = abs(audio->raw_audio[0][0]);
	
	unsigned int start_sample = time_to_sample_number(audio, start);
	unsigned int end_sample = time_to_sample_number(audio, start+duration);

	for ( unsigned int i=start_sample; i<end_sample; i++ ) {
		if ( abs(audio->raw_audio[0][i]) > n ) {
			n = abs(audio->raw_audio[0][i]);
		}
	}

	return (double)n/(double)DBL_MAX;
}
