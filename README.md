## W3C Model: The Wong-Wang-Wilson-Cowan hybrid model.

This model takes the form of
```
   dS_i,E/dt = - S_i,E/tau_E + (1-S_i,E) * gamma_E * H_E(w_EE.*S_i,E - w_IE.*S_i,I + I_E + I_i,G)
   dS_i,I/dt = - S_i,I/tau_I + (1-S_i,I) * gamma_I * H_I(w_EI.*S_i,E - w_II.*S_i,I + I_I)
 
 where the global input
   I_i,G = G * sum_j [C_ij*S_i,E] where C_ij = 0 for i=j
 
 and the transfer function
   H_p(x) = (r_max + (a_p*x-b_p-r_max) / (1-exp(d_p * (a_p*x - b_p - r_max))) )/ (1-exp(-d_p*(a_p*x-b_p)))
 
 for p=E or I.
 
 ```

The model can be thought of as a high-dimensional generalization of
 Wong-Wang model which is more consistent with the form of the
 Wilson-Cowan model. Importantly, the model use a sigmoidal transfer
 function that matches the linear threshold transfer function as used in
 Wong-Wang (2006)/Deco et al (2014JON) asymptotically for lower level of
 of activity. Note that sigmoidal transfer function will cap the firing
 rate such that the model is more realistic. 

You can explore the global bifurcation using the file bifurcation_global_model.m, 
which systematically study the bifurcation of the global W3C model. However, please note that this program could take more than a week to run. If you
 simply want to get the gist of the computation of bifurcation diagrams, please see the computation for the local model. 


Please cite this paper if you use this code:
Zhang, Mengsen, Yinming Sun, and Manish Saggar. "Cross-attractor repertoire provides new perspective on structure-function relationship in the brain." NeuroImage 259 (2022): 119401.
https://doi.org/10.1016/j.neuroimage.2022.119401
