
mainCLV3activation.m is the main file for simulations of CLV3 expression
 with mutant or WT CRM. 
 str_mutant.m determines which mutant to simulate and it includes mutants of
'WT', '950M','970M','DM', '970Mi' or '1060i'. 

distance_matrix.m generates matrices for the distance between different cis-elements.

stochastic_events_probs.m determines the probabilities of possible events at
any state. Cooperativity is also determined in this file. In the default setup, k_on 
is increased and k_off is decreased if there is an interaction
between cis-elements. The amount of this change is determined by
distance_matrix and f_d.

plot_all.m plots the results generated for multiple cis-elements after running 
mainCLV3activation for each mutant with multiple cis-elements ('WT', '950M','970M','DM').

plot_single.m plots the results generated for single cis-elements (970Mi','1060i').

