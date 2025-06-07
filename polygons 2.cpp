// Polygonal Waveform Synth created by Fabian Martinez
// https://github.com/FabianMartinez    
// This code is licensed under the MIT License.

/*
============
MIT License
============

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.







*/






#include <math.h>
#include <new>
#include <distingnt/api.h>

// Fix for M_PI not being defined on some toolchains
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// Utility
inline float clamp(float x, float min, float max) {
    return x < min ? min : (x > max ? max : x);
}
inline float frac(float x) {
    return x - floorf(x);
}

// Standard waveforms
inline float sine_wave(float phase) {
    return sinf(phase);
}
inline float triangle_wave(float phase) {
    return 2.0f * fabsf(2.0f * (phase / (2.0f * (float)M_PI) - floorf(phase / (2.0f * (float)M_PI) + 0.5f))) - 1.0f;
}
inline float saw_wave(float phase) {
    return 2.0f * (phase / (2.0f * (float)M_PI) - floorf(phase / (2.0f * (float)M_PI) + 0.5f));
}
inline float square_wave(float phase) {
    return (sinf(phase) >= 0.0f) ? 1.0f : -1.0f;
}

// --- Wavefolder ---
inline float wavefold(float x, float amount, float offset) {
    x += offset;
    float a = clamp(amount, 0.01f, 2.0f);
    while (x > a || x < -a) {
        if (x > a) x = 2.0f * a - x;
        else if (x < -a) x = -2.0f * a - x;
    }
    x -= offset;
    return x;
}

// Polygonal waveform core
float polygon_amplitude(float phi, float n, float T) {
    float mod_arg = frac((phi * n) / (2.0f * (float)M_PI));
    float angle = 2.0f * (float)M_PI * mod_arg - (float)M_PI + T * n * 2.0f * (float)M_PI / n;
    return cosf(angle);
}

float polygon_waveform(float phi, float n, float T, float phiOffset) {
    float amp = polygon_amplitude(phi, n, T);
    float theta = phi + phiOffset;
    return amp * cosf(theta);
}

// 2D polyhedron (polygon) shape for oscilloscope XY, with tilt and rotation
void polyhedron_xy(float phase, int faces, float tilt, float rotation, float& x, float& y) {
    float angle = phase;
    float step = 2.0f * (float)M_PI / faces;
    int i0 = (int)(angle / step) % faces;
    int i1 = (i0 + 1) % faces;
    float t = (angle - (i0 * step)) / step;

    float a0 = i0 * step;
    float a1 = i1 * step;

    float tiltAmount = clamp(tilt, 0.0f, 1.0f);
    float r = 1.0f - tiltAmount;
    float x0 = r * cosf(a0);
    float y0 = r * sinf(a0);
    float x1 = r * cosf(a1);
    float y1 = r * sinf(a1);

    float xi = (1.0f - t) * x0 + t * x1;
    float yi = (1.0f - t) * y0 + t * y1;

    float rot = rotation;
    float xr = xi * cosf(rot) - yi * sinf(rot);
    float yr = xi * sinf(rot) + yi * cosf(rot);

    x = xr;
    y = yr;
}

// Lissajous figure
void lissajous_xy(float phase, float ratioX, float ratioY, float freq, float& x, float& y) {
    x = sinf(phase * ratioX);
    y = sinf(phase * ratioY + freq);
}

// Star shape (2D) - star with triangles around a polyhedron
void star_xy(float phase, int points, float tilt, float rotation, float& x, float& y) {
    float baseR = 0.7f;
    float spikeR = 1.0f;
    float k = (float)points;
    float starMod = fabsf(triangle_wave(k * phase));
    float r = baseR + (spikeR - baseR) * starMod;
    float rTilt = r * (1.0f - tilt);
    float xr = rTilt * cosf(phase);
    float yr = rTilt * sinf(phase);
    float rot = rotation;
    x = xr * cosf(rot) - yr * sinf(rot);
    y = xr * sinf(rot) + yr * cosf(rot);
}

// Heart shape (2D)
void heart_xy(float phase, float rotation, float& x, float& y) {
    float hx = 16.0f * powf(sinf(phase), 3);
    float hy = 13.0f * cosf(phase) - 5.0f * cosf(2.0f * phase) - 2.0f * cosf(3.0f * phase) - cosf(4.0f * phase);
    hx /= 17.0f;
    hy /= 17.0f;
    x = hx * cosf(rotation) - hy * sinf(rotation);
    y = hx * sinf(rotation) + hy * cosf(rotation);
}

// Biciclette shape (2D) - precision controlled by 'segments' parameter
void biciclette_xy(float phase, float rotation, float& x, float& y, int segments = 8) {
    if (segments < 8) segments = 8;

    float wheelRadius = 0.3f;
    float frameLen = 1.0f;
    float wheel1x = -frameLen / 2.0f;
    float wheel2x = frameLen / 2.0f;
    float wheelY = -0.5f;
    float pedal = sinf(phase * 2.0f);
    float bikerY = 0.2f + 0.1f * pedal;

    float seg = 2.0f * M_PI / segments;
    float bx = 0.0f, by = 0.0f;

    int segIdx = (int)(phase / seg);
    float t = (phase - segIdx * seg) / seg;

    if (segIdx == 0) {
        bx = wheel1x + wheelRadius * cosf(t * 2.0f * M_PI);
        by = wheelY + wheelRadius * sinf(t * 2.0f * M_PI);
    } else if (segIdx == 1) {
        bx = wheel2x + wheelRadius * cosf(t * 2.0f * M_PI);
        by = wheelY + wheelRadius * sinf(t * 2.0f * M_PI);
    } else if (segIdx == 2) {
        bx = wheel1x * (1.0f - t) + wheel2x * t;
        by = wheelY;
    } else if (segIdx == 3) {
        bx = wheel2x * (1.0f - t);
        by = wheelY + t * 0.5f;
    } else if (segIdx == 4) {
        bx = wheel2x + t * 0.2f;
        by = wheelY + 0.5f + t * 0.2f;
    } else if (segIdx == 5) {
        bx = (wheel2x * (1.0f - t) + (wheel2x + 0.2f) * t);
        by = (wheelY + 0.5f) * (1.0f - t) + (wheelY + 0.7f) * t;
    } else if (segIdx == 6) {
        bx = 0.0f;
        by = (wheelY + 0.7f) * (1.0f - t) + (bikerY + 0.7f) * t;
    } else if (segIdx == 7) {
        float headR = 0.08f;
        bx = 0.0f + headR * sinf(t * 2.0f * M_PI);
        by = (bikerY + 0.8f) + headR * cosf(t * 2.0f * M_PI);
    } else {
        float tt = (float)(segIdx - 8 + t) / (segments - 8);
        bx = (0.0f) * (1.0f - tt) + (wheel1x + wheelRadius) * tt;
        by = (bikerY + 0.8f) * (1.0f - tt) + (wheelY) * tt;
    }

    x = bx * cosf(rotation) - by * sinf(rotation);
    y = bx * sinf(rotation) + by * cosf(rotation);
}

// Sun shape (2D): circle with star rays around it, with user-defined circle radius and rotation
void sun_xy(float phase, int points, float circleR, float tilt, float rotation, float& x, float& y) {
    float rayR = 1.0f;
    float k = (float)points;
    float rayMod = 0.5f * (1.0f + triangle_wave(k * phase));
    float r = circleR + (rayR - circleR) * rayMod;
    float rTilt = r * (1.0f - tilt);
    float xr = rTilt * cosf(phase);
    float yr = rTilt * sinf(phase);
    float rot = rotation;
    x = xr * cosf(rot) - yr * sinf(rot);
    y = xr * sinf(rot) + yr * cosf(rot);
}

// --- Nested Circles shape ---
void nested_circles_xy(float phase, int v1, int v2, int v3, float freq, float pwmFreq, float& x, float& y) {
    float t = phase;
    float fv1 = v1 / 1000.0f, fv2 = v2 / 1000.0f, fv3 = v3 / 1000.0f;
    auto S = [](int k, int n, float t, float freq) {
        float duty = 1.0f / n;
        float phase = freq * t - (float)k / n;
        phase -= (int)phase; // approximate fmodf(phase, 1.0f)
        return 0.5f + 0.5f * ((phase < duty) ? 1.0f : -1.0f);
    };
    float s = fv1 * S(0, 3, t, pwmFreq) + fv2 * S(1, 3, t, pwmFreq) + fv3 * S(2, 3, t, pwmFreq);
    x = s * cosf(freq * t);
    y = s * sinf(freq * t);
}

// --- Snail shape ---
void snail_xy(float phase, float freq, int s, float& x, float& y) {
    float t = phase;
    float fs = s / 1000.0f;
    float rampval = (freq / 2.0f) * t;
    rampval -= (int)rampval; // approximate fmodf(rampval, 1.0f)
    float factor = (1.0f + rampval) / 2.0f;
    x = fs * cosf(freq * t) * factor;
    y = fs * sinf(freq * t) * factor;
}

// --- Lemon shape ---
void lemon_xy(float phase, float freq, int s, float& x, float& y) {
    float t = phase;
    float fs = s / 1000.0f;
    float condval = 6.0f * freq * t;
    condval -= (int)condval; // approximate fmodf(condval, 1.0f)
    float cond = condval < 0.8f ? 1.0f : 0.0f;
    x = fs * cosf(freq * t) * cond;
    y = fs * sinf(freq * t) * cond;
}

// --- Butterfly shape ---
void butterfly_xy(float phase, float freq, float& x, float& y) {
    float t = phase;
    float f = freq * t;
    float rho = expf(cosf(f)) - 2.0f * cosf(4.0f * f) - powf(sinf(f / 12.0f), 5.0f);
    x = 0.25f * rho * sinf(f);
    y = 0.25f * rho * cosf(f);
}

// Algorithm struct
struct _polygonal_waveform : public _NT_algorithm
{
    float phase;
    float samplerate;
    _polygonal_waveform() : phase(0.0f), samplerate(48000.0f) {}
    ~_polygonal_waveform() {}
};

// Waveform names
static char const * const enumStringsWaveform[] = {
    "Polygonal",      // 0
    "Sine",           // 1
    "Triangle",       // 2
    "Sawtooth",       // 3
    "Square",         // 4
    "Polyhedron",     // 5
    "Lissajous",      // 6
    "Star",           // 7
    "Heart",          // 8
    "Biciclette",     // 9
    "Sun",            // 10
    "Nested Circles", // 11
    "Snail",          // 12
    "Lemon",          // 13
    "Butterfly"       // 14
};

// Parameter indices
enum {
    PARAM_OUTPUT1 = 0,
    PARAM_OUTPUT1_MODE,
    PARAM_OUTPUT2,
    PARAM_OUTPUT2_MODE,
    PARAM_POLY_FACES,
    PARAM_POLY_ROTATION,
    PARAM_FREQ,
    PARAM_AMP,
    PARAM_PHASEOFFSET,
    PARAM_WAVEFORM,
    PARAM_POLYHEDRON_TILT,
    PARAM_WAVEFOLD_ENABLE,
    PARAM_WAVEFOLD_AMOUNT,
    PARAM_WAVEFOLD_OFFSET,
    PARAM_LISSAJOU_X_RATIO,
    PARAM_LISSAJOU_Y_RATIO,
    PARAM_SUN_RADIUS,

    // --- Nested Circles ---
    PARAM_NC_V1,
    PARAM_NC_V2,
    PARAM_NC_V3,
    PARAM_NC_FREQ,
    PARAM_NC_PWM_FREQ,

    // --- Snail ---
    PARAM_SNAIL_FREQ,
    PARAM_SNAIL_S,

    // --- Lemon ---
    PARAM_LEMON_FREQ,
    PARAM_LEMON_S,

    // --- Butterfly ---
    PARAM_BUTTERFLY_FREQ,

    NUM_PARAMS
};

// Parameter defaults/min/max
static const _NT_parameter paramDefs[NUM_PARAMS] = {
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Output 1", 1, 13)
    NT_PARAMETER_AUDIO_OUTPUT_WITH_MODE("Output 2", 2, 14)
    { "Faces/Vertices", 3, 32, 6, kNT_unitNone, 0, NULL },
    { "Rotation", 0, 360, 0, kNT_unitNone, 0, NULL },
    { "Frequency", 1, 1000, 110, kNT_unitHz, 0, NULL },
    { "Amplitude", 0, 5, 5, kNT_unitVolts, 0, NULL },
    { "Phase Offset", 0, 360, 90, kNT_unitNone, 0, NULL },
    { "Waveform", 0, 14, 0, kNT_unitEnum, 0, enumStringsWaveform },
    { "Polyhedron Tilt", 0, 100, 0, kNT_unitNone, 0, NULL },
    { "Wavefold Enable", 0, 1, 0, kNT_unitNone, 0, NULL },
    { "Wavefold Amount", 1, 200, 100, kNT_unitNone, 0, NULL },
    { "Wavefold Offset", -100, 100, 0, kNT_unitNone, 0, NULL },
    { "Lissajous X Ratio", 1, 20, 1, kNT_unitNone, 0, NULL },
    { "Lissajous Y Ratio", 1, 20, 2, kNT_unitNone, 0, NULL },
    { "Sun Circle Radius", 10, 100, 70, kNT_unitNone, 0, NULL },

    // --- Nested Circles ---
    { "NC v1", 0, 1000, 200, kNT_unitNone, 0, NULL },
    { "NC v2", 0, 1000, 500, kNT_unitNone, 0, NULL },
    { "NC v3", 0, 1000, 800, kNT_unitNone, 0, NULL },
    { "NC Freq", 1, 100, 22, kNT_unitHz, 0, NULL },
    { "NC PWM Freq", 1, 100, 33, kNT_unitHz, 0, NULL },

    // --- Snail ---
    { "Snail Freq", 1, 60, 11, kNT_unitHz, 0, NULL },
    { "Snail s", 0, 1000, 400, kNT_unitNone, 0, NULL },

    // --- Lemon ---
    { "Lemon Freq", 1, 60, 11, kNT_unitHz, 0, NULL },
    { "Lemon s", 0, 1000, 200, kNT_unitNone, 0, NULL },

    // --- Butterfly ---
    { "Butterfly Freq", 1, 30, 4, kNT_unitHz, 0, NULL }
};

void calculateRequirements(_NT_algorithmRequirements& req, const int32_t* specifications)
{
    req.numParameters = NUM_PARAMS;
    req.sram = sizeof(_polygonal_waveform);
    req.dram = 0;
    req.dtc = 0;
    req.itc = 0;
}

_NT_algorithm* construct(const _NT_algorithmMemoryPtrs& ptrs, const _NT_algorithmRequirements& req, const int32_t* specifications)
{
    _polygonal_waveform* alg = new (ptrs.sram) _polygonal_waveform();
    alg->parameters = paramDefs;
    return alg;
}

void step(_NT_algorithm* self_, float* busFrames, int numFramesBy4)
{
    _polygonal_waveform* self = (_polygonal_waveform*)self_;
    int numFrames = numFramesBy4 * 4;
    float* out1 = busFrames + (self->v[PARAM_OUTPUT1] - 1) * numFrames;
    float* out2 = busFrames + (self->v[PARAM_OUTPUT2] - 1) * numFrames;
    bool replace1 = self->v[PARAM_OUTPUT1_MODE];
    bool replace2 = self->v[PARAM_OUTPUT2_MODE];

    int waveform = self->v[PARAM_WAVEFORM];
    int faces = self->v[PARAM_POLY_FACES];
    float rotation = (self->v[PARAM_POLY_ROTATION] / 360.0f) * (2.0f * (float)M_PI);
    float freq = self->v[PARAM_FREQ];
    float amp = self->v[PARAM_AMP];
    float tilt = self->v[PARAM_POLYHEDRON_TILT] / 100.0f;
    float lissajousX = self->v[PARAM_LISSAJOU_X_RATIO];
    float lissajousY = self->v[PARAM_LISSAJOU_Y_RATIO];
    float sunRadius = self->v[PARAM_SUN_RADIUS] / 100.0f;

    // Read new shape parameters from UI
    int nc_v1 = self->v[PARAM_NC_V1];
    int nc_v2 = self->v[PARAM_NC_V2];
    int nc_v3 = self->v[PARAM_NC_V3];
    float nc_freq = self->v[PARAM_NC_FREQ];
    float nc_pwm = self->v[PARAM_NC_PWM_FREQ];
    float snail_freq = self->v[PARAM_SNAIL_FREQ];
    int snail_s = self->v[PARAM_SNAIL_S];
    float lemon_freq = self->v[PARAM_LEMON_FREQ];
    int lemon_s = self->v[PARAM_LEMON_S];
    float butterfly_freq = self->v[PARAM_BUTTERFLY_FREQ];

    float phaseOffsetParam = self->v[PARAM_PHASEOFFSET];
    float phiOffset;
    if (waveform == 5) {
        phiOffset = ((phaseOffsetParam - 90.0f) / 360.0f) * (2.0f * (float)M_PI);
    } else if (waveform == 0) {
        phiOffset = (phaseOffsetParam / 360.0f) * (2.0f * (float)M_PI);
    } else {
        phiOffset = (phaseOffsetParam / 360.0f) * (2.0f * (float)M_PI);
    }

    bool wavefoldEnable = self->v[PARAM_WAVEFOLD_ENABLE] != 0;
    float wavefoldAmount = self->v[PARAM_WAVEFOLD_AMOUNT] / 100.0f;
    float wavefoldOffset = self->v[PARAM_WAVEFOLD_OFFSET] / 100.0f;

    float phase = self->phase; // phase in [0, 1)
    float sr = self->samplerate;
    float phaseInc = freq / sr;

    for (int i = 0; i < numFrames; ++i) {
        phase += phaseInc;
        if (phase >= 1.0f) phase -= 1.0f;

        float phaseRad = phase * 2.0f * (float)M_PI; // use this for all waveform functions

        float val1 = 0.0f, val2 = 0.0f;
        switch (waveform) {
            case 1: val1 = sine_wave(phaseRad + phiOffset); val2 = sine_wave(phaseRad); break;
            case 2: val1 = triangle_wave(phaseRad + phiOffset); val2 = triangle_wave(phaseRad); break;
            case 3: val1 = saw_wave(phaseRad + phiOffset); val2 = saw_wave(phaseRad); break;
            case 4: val1 = square_wave(phaseRad + phiOffset); val2 = square_wave(phaseRad); break;
            case 5: {
                float x, y;
                polyhedron_xy(phaseRad + phiOffset, faces, tilt, rotation, x, y);
                if (wavefoldEnable) {
                    x = wavefold(x, wavefoldAmount, wavefoldOffset);
                    y = wavefold(y, wavefoldAmount, wavefoldOffset);
                }
                val1 = x;
                val2 = y;
                break;
            }
            case 6: {
                float x, y;
                lissajous_xy(phaseRad, lissajousX, lissajousY, freq * 0.01f, x, y);
                float xr = x * cosf(rotation) - y * sinf(rotation);
                float yr = x * sinf(rotation) + y * cosf(rotation);
                val1 = xr;
                val2 = yr;
                break;
            }
            case 7: {
                float x, y;
                star_xy(phaseRad, faces, tilt, rotation, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 8: {
                float x, y;
                heart_xy(phaseRad, rotation, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 9: {
                float x, y;
                biciclette_xy(phaseRad, rotation, x, y, faces);
                val1 = x;
                val2 = y;
                break;
            }
            case 10: {
                float x, y;
                sun_xy(phaseRad, faces, sunRadius, tilt, rotation, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 11: { // Nested Circles
                float x, y;
                nested_circles_xy(phase, nc_v1, nc_v2, nc_v3, nc_freq, nc_pwm, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 12: { // Snail
                float x, y;
                snail_xy(phase, snail_freq, snail_s, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 13: { // Lemon
                float x, y;
                lemon_xy(phase, lemon_freq, lemon_s, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 14: { // Butterfly
                float x, y;
                butterfly_xy(phase, butterfly_freq, x, y);
                val1 = x;
                val2 = y;
                break;
            }
            case 0:
            default: {
                float n = (float)faces;
                float T = rotation / (2.0f * (float)M_PI);
                val1 = polygon_waveform(phaseRad, n, T, phiOffset);
                val2 = polygon_waveform(phaseRad, n, T, 0.0f);
                break;
            }
        }

        float sample1 = clamp(val1 * (amp / 2.0f), -amp / 2.0f, amp / 2.0f);
        float sample2 = clamp(val2 * (amp / 2.0f), -amp / 2.0f, amp / 2.0f);

        if (!replace1)
            out1[i] += sample1;
        else
            out1[i] = sample1;

        if (!replace2)
            out2[i] += sample2;
        else
            out2[i] = sample2;
    }
    self->phase = phase;
}

static const _NT_factory factory = {
    .guid = NT_MULTICHAR('P','O','L','2'),
    .name = "Polygonal Waveform",
    .description = "Polygonal/Polyhedron/Lissajous/Star/Heart/Biciclette/Sun/Nested Circles/Snail/Lemon/Butterfly Waveform Synth (Stereo)",
    .numSpecifications = 0,
    .specifications = NULL,
    .calculateStaticRequirements = NULL,
    .initialise = NULL,
    .calculateRequirements = calculateRequirements,
    .construct = construct,
    .parameterChanged = NULL,
    .step = step,
    .draw = NULL,
    .midiRealtime = NULL,
    .midiMessage = NULL,
    .tags = 0
};

extern "C" {
uintptr_t pluginEntry(_NT_selector selector, uint32_t data)
{
    switch (selector)
    {
    case kNT_selector_version:
        return kNT_apiVersionCurrent;
    case kNT_selector_numFactories:
        return 1;
    case kNT_selector_factoryInfo:
        if (data == 0)
            return (uintptr_t)&factory;
        return 0;
    }
    return 0;
}
}