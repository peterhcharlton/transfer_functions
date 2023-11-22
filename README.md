# transfer_functions
Transfer functions for photoplethymogram and blood pressure signals

This repository contains Matlab scripts to:
1. Estimate arterial blood pressure waveforms from finger photoplethysmogram waveforms; and
2. Estimate central blood pressure waveforms from peripheral blood pressure waveforms.

## Summary 

Arterial blood pressure (ABP) waveforms contain a wealth of information on the cardiovascular (CV) system, providing potential diagnostic utility [[1]](#1). Furthermore, ABP waveforms can be modelled using physical principles [[2]](#2), allowing one to predict the shape of the ABP pulse under specified CV conditions. However, methods for acquiring ABP signals in vivo either require a skilled operator (_e.g._ applanation tonometry), or are invasive (_e.g._ a pressure catheter). Conversely, the digital photoplethysmogram (PPG) can be easily and non-invasively acquired in a wide range of clinical settings, using a pulse oximetry probe. It is closely related to the BP waveform [[3]](#3), although the exact physiological mechanisms underlying it remain unclear [[4]](#4). A generalised transfer function relating the digital PPG to the ABP waveform has previously been reported [[5]](#5), which can be used to estimate ABP waveforms from the PPG. This transfer function could allow one to conduct ABP-based analyses on estimated ABP signals, derived from easily acquired PPG waveforms. This would expand the utility of ABP-based analyses to clinical settings where it is not practical to acquire ABP waveforms, and would also allow one to use physical modelling to predict the shape of PPG waveforms.

This repository contains two Matlab scripts:
1. `PPG_ABP_TF`, for estimating digital and radial ABP waveforms from the PPG (and vice-versa). It is based on the transfer function reported in [[5]](#5).
2. `per_cen_TF`, for estimating central ABP waveforms from peripheral ABP waveforms.

## Estimating ABP from the PPG

It has been suggested that a transfer function can be used to describe the relationship between PPG and ABP waveforms, which does not differ between normotensive and hypertensive patients, nor during administration of nitroglycerin, a vasodilator which affects the shape of the pulse waveform [[5]](#5). The transfer functions calculated using the ratio of the fast Fourier transforms of waveforms: FFT(ABP) / FFT(PPG), and were described using a Bode plot as shown in Figure 2 of [[5]](#5).

To demonstrate how to use the transfer function, it is helpful to firstly show how to recover the PPG signal from its frequency domain representation. The PPG signal is denoted $x(t)$. The FFT of the PPG,

```math
	X(f) = \mathcal{F}\left(x(t)\right) 
```

is calculated, where $X(f)$ is the frequency domain representation of the PPG (a vector of complex numbers). The original PPG signal, $x(t)$, can be recovered from $X(f)$ using the following sum of cosine and sine waves:

```math
	x[i] = \sum_{k=0}^{N/2} Re \{\overline{X}[k]\} \cos\left(\frac{2\pi k i}{N}\right) \; + \; \sum_{k=0}^{N/2} Im\{\overline{X}[k]\} \sin\left(\frac{2\pi k i}{N}\right) \quad ,
```

where $i$ is the sample number (varying from 0 to $N-1$), $N$ is the length of the FFT, and

```math
	Re\{\overline{X}[k]\} =  \frac{Re\{X[k]\}}{N/2} \quad \mathrm{and} \quad Im\{\overline{X}[k]\} = -\frac{Im\{X[k]\}}{N/2} \quad ,
```

with the exceptions that

```math
Re\{\overline{X}[0]\} =  \frac{Re\{X[0]\}}{N} \quad \mathrm{and} \quad Re\{\overline{X}[N/2]\} = -\frac{Im\{X[N/2]\}}{N} \quad .
```

A similar approach can be used to estimate an ABP signal, $y(t)$, from the PPG, $x(t)$. Firstly, $X(f)$ is calculated from $x(t)$ using the FFT. Secondly, $y(t)$ is calculated using a sum of cosine and sine waves, where the amplitude and phase offset of each wave are adjusted according to the transfer function. For instance, 

```math
y[i] = \sum_{k=0}^{N/2} A[k] Re\{\overline{X}[k]\} \cos\left(\frac{2\pi k i}{N} + \phi[k] \right) \; + \; \sum_{k=0}^{N/2} A[k] Im\{\overline{X}[k]\} \sin\left(\frac{2\pi k i}{N} + \phi[k] \right) \; ,
```
where the amplitude and phase of the transfer function are denoted by $A(f)$ and $\phi(f)$ respectively, and $\phi$ is expressed in radians.

The reported transfer functions can also be used to estimate a PPG signal from the ABP. This can be achieved by using the reciprocal of the amplitude, and the negative of the phase, of the reported transfer function.


## References

<a id="1">[1]</a> 
M. O'Rourke and D. Gallagher, 'Pulse wave analysis', _Journal of Hypertension_, vol. 14, no. suppl 5, pp. S147-57, 1996.

<a id="2">[2]</a>
J. Alastruey, K. Parker, and S. Sherwin, 'Arterial pulse wave haemodynamics', in _Proc. BHR Group's 11th International Conference on Pressure Surges_, Lisbon, Portugal, 2012, pp. 401-43.

<a id="3">[3]</a>
A. A. Awad et al., How does the plethysmogram derived from the pulse oximeter relate to arterial blood pressure in coronary artery bypass graft patients?" _Anesthesia and Analgesia_, vol. 93, no. 6, pp. 1466-71, 2001.

<a id="4">[4]</a>
A. A. Kamshilin et al., 'A new look at the essence of the imaging photoplethysmography'. Scientific Reports, vol. 5, p. 10494, 2015.

<a id="5">[5]</a>
S. C. Millasseau et al., 'Noninvasive assessment of the digital volume pulse. Comparison with the peripheral pressure pulse'. _Hypertension_, vol. 36, no. 6, pp. 952-6, 2000.

<a id="6">[6]</a>
M. Karamanoglu et al., 'An analysis of the relationship between central aortic and peripheral upper limb pressure waves in man', _European Heart Journal_, vol. 14, no. 2, pp. 160-7, 1993.
