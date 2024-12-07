#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

/*-----------------------------------------------------------------------------*/
/* MACROS */

enum audio_positions {
	// Begin of the audio
	BEGIN =		0,

	// End of the audio
	END =	   -1,

	// Current position
	CURRENT =  -2
};

/*
Binary audio formats
LE: Little Endian
BE: Big Endian
*/
enum audio_formats {
	// Floating point
	F32LE,
	F32BE,
	F64LE,
	F64BE,

	// Signed integer
	S8,
	S16LE,
	S16BE,
	S24LE,
	S24BE,
	S32LE,
	S32BE,

	// Unsigned integer
	U8,
	U16LE,
	U16BE,
	U24LE,
	U24BE,
	U32LE,
	U32BE
};

enum volume_levels {
	MIN_LEVEL,
	MAX_LEVEL
};


/*-----------------------------------------------------------------------------*/
/* DATA STRUCTURE */

typedef struct {
	double**		raw_audio;
	unsigned int 	sample_rate;
	double			current_place;
	double 			duration;
	unsigned int 	channels;
	unsigned int 	sample_quantity_per_channel;
	unsigned int 	sample_quantity_total;
} Audio;

typedef struct {
	union sample** 	raw_audio;
	unsigned int 	sample_rate;
	unsigned int	current_place;
	double 			duration;
	unsigned int 	channels;
	unsigned int 	sample_quantity_per_channel;
	unsigned int 	sample_quantity_total;
	unsigned int	bytes_per_sample;
	unsigned int 	format;
} AudioTranscoded;


typedef struct {
	double max_level;
	double delay;
	double attack;
	double hold;
	double decay;
	double release;
} Sound;


/*-----------------------------------------------------------------------------*/


/*
Initialize the audio struct (Memory allocation and parameter settings)
 - audio:	 	Pointer to an uninitialized Audio struct
 - sample_rate: Sample rate
 - (format): 	Binary audio format (enum audio_formats)
 - duration: 	Duration of the audio in seconds
 - channels: 	Number of channels
*/
void init_audio ( Audio* audio, unsigned int sample_rate, double duration, unsigned int channels );
void init_transcoded_audio ( AudioTranscoded* audio, unsigned int sample_rate, unsigned int format, double duration, unsigned int channels );

/*
Freeing the memory reserved by the audio structs
- audio: Pointer to the Audio struct to be destroyed
*/
void free_audio ( Audio* audio );
void free_transcoded_audio( AudioTranscoded* audio );

/*-----------------------------------------------------------------------------*/
/* AUDIO GENERATORS */

/*
White noise generator
- audio: 		Pointer to an initialized Audio struct
- max_volume: 	Maximum volume (0.0 -> minimum volume, 1.0 -> maximum volume)
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void white_noise ( Audio* audio, double max_volume, double start, double duration );

/*
Sine wave generator
- audio:		Pointer to an initialized Audio struct
- frequency: 	Frequency of the sine wave in Hz
- max_volume:	Maximum volume (0.0 -> minimum volume, 1.0 -> maximum volume)
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void sine_wave ( Audio* audio, double frequency, double max_volume, double start, double duration );

/*
Sine wave generator
- audio:		Pointer to an initialized Audio struct
- frequency: 	Frequency of the sine wave in Hz
- max_volume:	Maximum volume (0.0 -> minimum volume, 1.0 -> maximum volume)
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void cosine_wave ( Audio* audio, double frequency, double max_volume, double start, double duration );

/*
Square wave generator
- audio: 		Pointer to an initialized Audio struct
- frequency: 	Frequency of the square wave in Hz
- max_volume:	Maximum volume (0.0 -> minimum volume, 1.0 -> maximum volume)
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void square_wave ( Audio* audio, double frequency, double max_volume, double start, double duration );

/*
Sawtooth wave generator
- audio: 		Pointer to an initialized Audio struct
- frequency: 	Frequency of the sawtooth wave in Hz
- max_volume: 	Maximum volume (0.0 -> minimum volume, 1.0 -> maximum volume)
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void sawtooth_wave ( Audio* audio, double frequency, double max_volume, double start, double duration );

/*
Triangular wave generator
- audio: 		Pointer to an initialized Audio struct
- frequency: 	Frequency of the triangular wave in Hz
- max_volume: 	Maximum volume (0.0 -> minimum volume, 1.0 -> maximum volume)
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void triangular_wave ( Audio* audio, double frequency, double max_volume, double start, double duration );

/*
Constant level generator (Generates a constant audio level)
- audio: 	Pointer to an initialized Audio struct
- level: 	Audio level
- start: 	Point in time at which the audio should start
- duration: Duration of the section
*/
void constant_level ( Audio* audio, double level, double start, double duration );

/*
Create silence
- audio: 	Pointer to an initialized Audio struct
- start: 	Point in time at which the audio should start
- duration: Duration of the section
*/
void silence ( Audio* audio, double start, double duration );



/*-----------------------------------------------------------------------------*/
/* AUDIO PROCESSORS */

/*
Amplitude modulator
- audio: 	Audio to be modulated by carrier
- carrier:  Carrier to modulate the audio
*/
void am_modulate ( Audio* audio, Audio* carrier );

/*
Amplifier
- audio: 		Pointer to an initialized Audio struct
- amp_factor:	Amplification factor
- start: 		Point in time at which the amplification should start
- duration: 	Duration of the section
*/
void amplify ( Audio* audio, double amp_factor, double start, double duration );

/*
Vertical shifter (Shift the audio wave up or down by a shift value)
- audio: 	Pointer to an initialized Audio struct
- shift: 	Shift value in seconds by which the audio wave should be shifted up or down
- start: 	Point in time at which the shift should start
- duration: Duration of the section
*/
void shift_vertically ( Audio* audio, double shift, double start, double duration );

/*
Vertical shifter (Shift the audio wave up or down by a shift value)
- audio: 	Pointer to an initialized Audio struct
- shift: 	Shift value in seconds by which the audio wave should be shifted up or down
- start: 	Point in time at which the shift should start
- duration: Duration of the section
*/
void shift_horizontally ( Audio* audio, double shift, double start, double duration );

/*
Volume limiter (Limits the audio volume to a specified volume)
- audio: 		Pointer to an initialized Audio struct
- max_volume: 	Maximum volume level
- start: 		Point in time at which the audio should start
- duration: 	Duration of the section
*/
void limit_volume ( Audio* audio, double max_volume, double start, double duration );

/*
Resampler (Change the audio sample rate)
- dest_audio: 	Pointer to an uninitialized Audio struct to store the resampled audio
- src_audio:	Pointer to an initialized Audio struct to be resampled
- sample_rate: 	New sample rate of the audio
*/
void resample ( Audio* dest_audio, Audio* src_audio, unsigned int sample_rate );

/*
Insert audio into another audio
- audio: 			Audio to insert paste_audio into
- paste_audio:		Audio to be inserted into audio
- start_audio: 		Starting point of audio to insert from paste_audio
- start_paste: 		Starting point of paste_audio to insert into audio
- duration_paste: 	Duration of the section to be inserted
*/
void insert_audio ( Audio* audio, Audio* paste_audio, double start_audio, double start_paste, double duration_paste );


/*
Append audio
- dest_audio: 			Pointer to the destination Audio struct to which the audio should be appended
- src_audio: 			Pointer to the source Audio struct from which the audio should be appended to dest_audio
- start_src_audio: 		Starting point of the audio excerpt from src_audio
- duration_src_audio: 	Duration of the audio excerpt from src_audio
*/
void append_audio ( Audio* dest_audio, Audio* src_audio, double start_src_audio, double duration_src_audio );

/*
Append audio
- audio: 				Pointer to the destination Audio struct to which the audio should be appended
- sound: 				Pointer to the source Audio struct from which the audio should be appended to dest_audio
- frequency: 			Starting point of the audio excerpt from src_audio
- generator_function: 	Duration of the audio excerpt from src_audio
*/
void append_sound (Audio* audio, Sound* sound, double frequency, void(*generator_function)(Audio*, double, double, double, double));

/*
Audio fade in
- audio: 	Pointer to an initialized Audio struct
- duration: Duration of the fade in
- to_level: Volume level to which the audio should be faded in
- start: 	Point in time at which the fade in should start
*/
void fade_in ( Audio* audio, double duration, double to_level, double start );

/*
Audio fade out
- audio: 	Pointer to an initialized Audio struct
- duration: Duration of the fade out
- to_level: Volume level to which the audio should be faded out
- start: 	Point in time at which the fade out should start
*/
void fade_out( Audio* audio, double duration, double to_level, double start );

/*
Make all audio values positive
- audio: 	Pointer to an initialized Audio struct
- start:	Point in time at which the audio should start
- duration: Duration of the section
*/
void make_positive ( Audio* audio, double start, double duration );


/*
Make all audio values negative
- audio: 	Pointer to an initialized Audio struct
- start: 	Point in time at which the audio should start
- duration: Duration of the section
*/
void make_negative ( Audio* audio, double start, double duration );

/*
Audio inverter
- audio: 	Pointer to an initialized Audio struct
- start:	Point in time at which the audio should start
- duration: Duration of the section
*/
void invert ( Audio* audio, double start, double duration );

/*
Audio compressor (Set all samples louder than max_level to max_level)
- audio: 	 Pointer to an initialized Audio struct
- max_level: Maximum audio level after which the volume of louder samples are set to max_level
- start: 	 Point in time at which the audio should start
- duration:  Duration of the section
*/
void compress ( Audio* audio, double max_level, double start, double duration );

/*
Add two audios together
- audio: 		Pointer to an initialized Audio struct into which to store the audio sum
- audio_to_add: Pointer to an initialized Audio struct (Audio to add to audio)
*/
void add_audio ( Audio* audio, Audio* audio_to_add );

/*
Subtract audio from other audio
- audio: 				Pointer to an initialized Audio struct into which to store the audio difference
- audio_to_subtract: 	Pointer to an initialized Audio struct (Audio to subtract from audio)
*/
void subtract_audio ( Audio* audio, Audio* audio_to_subtract );


/*-----------------------------------------------------------------------------*/
/* INPUT / OUTPUT */

/*
File writer (Write the raw audio data to a file)
- audio: 		Pointer to an initialized Audio struct to be written to a file
- sample_rate:	Sample rate to transcode the audio to before writing it to a file
- format:		Binary audio format for writing the audio to a file (enum audio_formats)
- file_path: 	File path of the which into which the raw audio data should be written
*/
void write_audio_to_file ( Audio* audio, unsigned int format, unsigned int sample_rate, char* file_path );

/*
Output the raw audio data to stdout
- audio: Pointer to an initialized Audio struct to be outputted onto stdout
*/
void print_binary_audio ( Audio* audio );

/*
File reader (Read raw audio from file into Audio struct)
- audio: 		Pointer to an initialized Audio struct into which the audio data from the file should be loaded
- file_path: 	File path of the file from which the raw audio data should be loaded
- sample_rate: 	Sample rate per channel of the raw audio
- format:		Binary audio format for writing the audio to a file (enum audio_formats)
- channels: 	Number of audio channels
*/
void load_audio_from_file ( Audio* audio, char* file_path, unsigned int sample_rate, unsigned int format, unsigned int channels );


/*-----------------------------------------------------------------------------*/
/* TRANSCODER */

/*
Transcode the audio to another binary format
- dest_audio:	Pointer to an uninitialized AudioTranscoded struct to write the transcoded audio into
- src_audio:	Pointer to an initialized Audio struct to be transcoded
- format:		Binary audio format for writing the audio to a file (enum audio_formats)
*/
void transcode_raw_audio ( AudioTranscoded* dest_audio, Audio* src_audio, unsigned int format );

/*-----------------------------------------------------------------------------*/
/* INFORMATION */

/*
Get the maximum audio level (Distance from 0.0)
- audio: Pointer to an initialized Audio struct
- start: Point in time at which the audio should start
- duration: Duration of the section
*/
double get_highest_level ( Audio* audio, double start, double duration );
