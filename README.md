# Phase-Response-Curves

This repository contains code for three Neurospora crassa circadian clock models, and functions for doing basic phase response analysis. The three models coded are:
1.	An original simple model (1)
2.	A model by Dovzhenok et al. (2)
3.	A model by Tseng et al. (3)
The three initialize files (InitializeModDovzhenok, InitializeModSimple, and InitializeModTseng) include the ode models, the parameters, and initial values of each model necessary to run on the limit cycle. After the model is initialized in each of these files (using the ODE class), a simple test case is run, where the model is run for 100 hours on the limit cycle starting at frq mRNA max, and then frq mRNA is plotted over time.
The fourth file in this repository, ODE, defines the ODE class that has been built specifically for running circadian oscillators. Various methods in this class allow for solving the model, finding the period of the limit cycle, evaluating PRCs, etc. The methods are defined in alphabetic order, and the most crucial method in this class is dirPRC, which evaluates the PRC of a model using the direct method given the parameters of a pulse (the times to apply the pulses, which parameters are pulsed, the pulse amplitudes, and the duration of the pulse).

References
1.	Bellman, J. 2016. Phase Response Optimization of the Circadian Clock in Neurospora crassa. University of Cincinnati.
2.	Dovzhenok, A. A., M. Baek, S. Lim, and C. I. Hong. 2015. Mathematical modeling and validation of glucose compensation of the neurospora circadian clock. Biophys J 108:1830-1839.
3.	Tseng, Y. Y., S. M. Hunt, C. Heintzen, S. K. Crosthwaite, and J. M. Schwartz. 2012. Comprehensive modelling of the Neurospora circadian clock and its temperature compensation. PLoS Comput Biol 8:e1002437.
