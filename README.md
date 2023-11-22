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
