#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif

/*	A highpass filter made by subtracting the lowpass result from the 
	input sample.
*/ 
static double subnormal = (1.0 / 4294967295.0);
static t_class *highpass_class;

typedef struct _highpass
{
    t_object x_obj; 	
    t_float frequency;    	//the frequency of the filter
	t_float	q;	// the filter quality
	t_float sample_rate;
	t_float pole1;	//the poles
	t_float pole2;
	t_float pole3;
	t_float pole4;
	t_float out1; //variable to store the last output
} t_highpass;


//The method for sample processing
static t_int *highpass_perform(t_int *w)
{
	t_highpass *x = (t_highpass *)(w[1]);
    t_float *in = (t_float *)(w[2]);
    t_float *out = (t_float *)(w[3]);
    int n = (int)(w[4]);

	int sample_number = 0;
	float c, input;
	//The calculation of c that we did in class
	c = exp(-6.283185307179586f * x->frequency / x->sample_rate);

	while (n--)
	{
		//The filtering we did in class, but now the result is subtracted from input
		input = *(in + sample_number) - x->pole4 * x->q;

		if (input > 1.0f) input = 1.0f;
		if (input < -1.0f) input = -1.0f;
		input = 1.5 * input - 0.5 * (input * input * input);

		x->pole1 = input * (1.0f - c) + x->pole1 * c +subnormal;
		x->pole2 = x->pole1 * (1.0f - c) + x->pole2 * c;
		x->pole3 = x->pole2 * (1.0f - c) + x->pole3 * c;
		x->pole4 = x->pole3 * (1.0f - c) + x->pole4 * c;

		//The subtraction causes the filter to do a highpass
		x->out1 = input - x->pole4;
		*(out + sample_number) = x->out1;

		sample_number++;
	}
	return (w + 5);
}


void highpass_q(t_highpass *x, t_floatarg g)
{	
	//check if it is a valid q
	if (g > 0)x->q = g;
	else x->q = 0;
}
    
static void highpass_dsp(t_highpass *x, t_signal **sp)
{
	//Initialize the variables in case they changed while dsp was off.
	x->pole1 = 0;
	x->pole2 = 0;
	x->pole3 = 0;
	x->pole4 = 0;
	x->out1 = 0;
	x->sample_rate = sp[0]->s_sr;
    dsp_add(highpass_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void *highpass_new(void)
{
    t_highpass *x = (t_highpass *)pd_new(highpass_class);
    inlet_new(&x->x_obj, &x->x_obj.ob_pd, gensym("float"), gensym("q"));
    outlet_new(&x->x_obj, gensym("signal"));
	//Initialize the variables when the object is first created.
	x->pole1 = 0;
	x->pole2 = 0;
	x->pole3 = 0;
	x->pole4 = 0;
	x->out1= 0;
    x->q = 0;
    return (x);
}

void highpass_tilde_setup(void)
{
	highpass_class = class_new(gensym("highpass~"), (t_newmethod)highpass_new,0,
		sizeof(t_highpass), 0, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(highpass_class, t_highpass, frequency);
    class_addmethod(highpass_class, (t_method)highpass_dsp, gensym("dsp"), (t_atomtype)0);
    class_addmethod(highpass_class, (t_method)highpass_q, gensym("q"), A_FLOAT, (t_atomtype)0);
}
