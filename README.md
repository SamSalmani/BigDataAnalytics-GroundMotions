# Big Data Analytics

In ehVSR.m, earthquake Horizintal to Vertical Spectral ratio can be computed using ground motion data. Intitiall all ground motion data in the desired location are derived using the code that is in Geound_Motion_Derivation repository. Then, eHVSRs are computed using this code.
Then, there is another code, HVSR_Inversion.m, that can be used to use Ensemble Kalman Inversion to compute Vs and Vp. For more information regarding the methodology, you can refer to "Albers, D. J., Blancquart, P. A., Levine, M. E., Seylabi, E. E., & Stuart, A. (2019). Ensemble Kalman methods with constraints. Inverse Problems, 35(9), 095007."
