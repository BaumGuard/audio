// Musical notes and their corresponding frequencies
// Source : https://auditoryneuroscience.com/pitch/fundamental-frequencies-notes-western-music

#define C0 16.35
#define CS0 17.32
#define D0 18.35
#define DS0 19.45
#define E0 20.60
#define F0 21.83
#define FS0 23.12
#define G0 24.50
#define GS0 25.96
#define A0 27.50
#define AS0 29.14
#define B0 30.87

#define C1 32.70
#define CS1 34.65
#define D1 36.71
#define DS1 38.89
#define E1 41.20
#define F1 43.65
#define FS1 46.25
#define G1 49.00
#define GS1 51.91
#define A1 55.00
#define AS1 58.27
#define B1 61.74

#define C2 65.41
#define CS2 69.30
#define D2 73.42
#define DS2 77.78
#define E2 82.41
#define F2 87.31
#define FS2 92.50
#define G2 98.00
#define GS2 103.83
#define A2 110.00
#define AS2 116.54
#define B2 123.47

#define C3 130.81
#define CS3 138.59
#define D3 146.83
#define DS3 155.56
#define E3 164.81
#define F3 174.61
#define FS3 185.00
#define G3 196.00
#define GS3 207.65
#define A3 220.00
#define AS3 233.08
#define B3 246.94

#define C4 261.63
#define CS4 277.18
#define D4 293.66
#define DS4 311.13
#define E4 329.63
#define F4 349.23
#define FS4 369.99
#define G4 392.00
#define GS4 415.30
#define A4 440.00
#define AS4 466.16
#define B4 493.88

#define C5 523.25
#define CS5 554.37
#define D5 587.33
#define DS5 622.25
#define E5 659.26
#define F5 698.46
#define FS5 739.99
#define G5 783.99
#define GS5 830.61
#define A5 880.00
#define AS5 932.33
#define B5 987.77

#define C6 1046.50
#define CS6 1108.73
#define D6 1174.66
#define DS6 1244.51
#define E6 1318.51
#define F6 1396.91
#define FS6 1479.98
#define G6 1567.98
#define GS6 1661.22
#define A6 1760.00
#define AS6 1864.66
#define B6 1975.53

#define C7 2093.00
#define CS7 2217.46
#define D7 2349.32
#define DS7 2489.02
#define E7 2637.02
#define F7 2793.83
#define FS7 2959.96
#define G7 3135.96
#define GS7 3322.44
#define A7 3520.00
#define AS7 3729.31
#define B7 3951.07

#define C8 4186.01
#define CS8 4434.92
#define D8 4698.64
#define DS8 4978.03
#define E8 5120.0
#define F8 5588.48
#define FS8 5918.72
#define G8 6272.0
#define GS8 6645.76
#define A8 7040.0
#define AS8 7459.84
#define B8 7902.72



// Note duration factors
#define WHOLE       1.0
#define HALF        0.5
#define QUARTER     0.25
#define EIGHTH      0.125
#define SIXTEENTH   0.0625


// Calculate the note duration
// bpm:                  Beats per minute
// note_duration_factor: Note duration factor
double note_duration ( double bpm, double note_duration_factor ) {
    return (60.0/bpm)*note_duration_factor;
}
