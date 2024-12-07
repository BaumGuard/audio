#include "audio.h"

/*
Amplitude Shift Keying (ASK)
- audio:            Pointer to an uninitialized Audio struct
- number_of_bits:   Number of bits to be encoded using ASK
- frequency:        Audio frequency used for encoding
- baud_rate:        Baud rate of the encoded signal
- level_0:          Volume level of the sine wave of an encoded 0
- level_1:          Volume level of the sine wave on an encoded 1
*/
void ask (Audio* audio, unsigned int number_of_bits, double frequency, double baud_rate, double level_0, double level_1);

/*
Frequency Shift Keying (FSK)
- audio:            Pointer to an uninitialized Audio struct
- number_of_bits:   Number of bits to be encoded using FSK
- baud_rate:        Baud rate of the encoded signal
- frequency_0:      Frequency of the sine wave of an encoded 0
- frequency_1:      Frequency of the sine wave on an encoded 1
*/
void fsk (Audio* audio, unsigned int number_of_bits, double baud_rate, double frequency_0, double frequency_1);

/*
Quadrature Phase Shift Keying (QPSK)
- audio:            Pointer to an uninitialized Audio struct
- number_of_bits:   Number of bits to be encoded using ASK
- baud_rate:        Baud rate of the encoded signal
- frequency:        Frequency of the sine wave and cosine wave used for encoding
*/
void qpsk (Audio* audio, unsigned int number_of_bits, double baud_rate, double frequency);

/*
Offset Quadrature Phase Shift Keying (QPSK)
- audio:            Pointer to an uninitialized Audio struct
- number_of_bits:   Number of bits to be encoded using ASK
- baud_rate:        Baud rate of the encoded signal
- frequency:        Frequency of the sine wave and cosine wave used for encoding
*/
void oqpsk (Audio* audio, unsigned int number_of_bits, double baud_rate, double frequency);
