#include "digital_modulations.h"

void ask (Audio* audio, unsigned int number_of_bits, double frequency, double baud_rate, double level_0, double level_1) {
    double duration = (double)number_of_bits * (1.0/baud_rate);

    init_audio(audio, 48000, duration, 1);

    int bit;

    for (int i=0; i<number_of_bits; i++) {
        bit = rand() % 2;

        if (bit == 0) {
            sine_wave(audio, frequency, level_0, CURRENT, 1.0/baud_rate);
        }
        else {
            sine_wave(audio, frequency, level_1, CURRENT, 1.0/baud_rate);
        }
    }
}


void fsk (Audio* audio, unsigned int number_of_bits, double baud_rate, double frequency_0, double frequency_1) {
    double duration = (double)number_of_bits * (1.0/baud_rate);

    init_audio(audio, 48000, duration, 1);

    int bit;

    for (int i=0; i<number_of_bits; i++) {
        bit = rand() % 2;

        if (bit == 0) {
            sine_wave(audio, frequency_0, 1.0, CURRENT, 1.0/baud_rate);
        }
        else {
            sine_wave(audio, frequency_1, 1.0, CURRENT, 1.0/baud_rate);
        }
    }
}


void qpsk (Audio* audio, unsigned int number_of_bits, double baud_rate, double frequency) {
    double duration = (number_of_bits/2) * (1.0/baud_rate);
    int bit;

    Audio i_channel, q_channel;
    init_audio(&i_channel, 48000, duration, 1);
    init_audio(&q_channel, 48000, duration, 1);
    init_audio(audio, 48000, duration, 1);

    for (int i=0; i<number_of_bits; i++) {
        bit = rand() % 2;

        if (i%2 == 0) {
            if (bit == 0) {
                cosine_wave(&i_channel, frequency, -0.5, CURRENT, 1.0/baud_rate);
            }
            else {
                cosine_wave(&i_channel, frequency, 0.5, CURRENT, 1.0/baud_rate);
            }
        }
        else {
            if (bit == 0) {
                sine_wave(&q_channel, frequency, -0.5, CURRENT, 1.0/baud_rate);
            }
            else {
                sine_wave(&q_channel, frequency, 0.5, CURRENT, 1.0/baud_rate);
            }
        }
    }

    add_audio(audio, &i_channel);
    add_audio(audio, &q_channel);

    free_audio(&i_channel);
    free_audio(&q_channel);
}


void oqpsk (Audio* audio, unsigned int number_of_bits, double baud_rate, double frequency) {
    double duration = (number_of_bits/2) * (1.0/baud_rate) + 0.5/baud_rate;
    int bit;

    Audio i_channel, q_channel;
    init_audio(&i_channel, 48000, duration, 1);
    init_audio(&q_channel, 48000, duration, 1);
    init_audio(audio, 48000, duration, 1);

    silence(&q_channel, CURRENT, 0.5/baud_rate);

    for (int i=0; i<number_of_bits; i++) {
        bit = rand() % 2;

        if (i%2 == 0) {
            if (bit == 0) {
                cosine_wave(&i_channel, frequency, -0.5, CURRENT, 1.0/baud_rate);
            }
            else {
                cosine_wave(&i_channel, frequency, 0.5, CURRENT, 1.0/baud_rate);
            }
        }
        else {
            if (bit == 0) {
                sine_wave(&q_channel, frequency, -0.5, CURRENT, 1.0/baud_rate);
            }
            else {
                sine_wave(&q_channel, frequency, 0.5, CURRENT, 1.0/baud_rate);
            }
        }
    }

    add_audio(audio, &i_channel);
    add_audio(audio, &q_channel);

    free_audio(&i_channel);
    free_audio(&q_channel);
}
