"""
gammadifference_sample.
A code example for using gammadifference.
"""

import gammadifference
import time


alpha_1, beta_1 = 50, 10000
alpha_2, beta_2 = 100, 20000
print "alpha_1, beta_1, alpha_2, beta_2 =",alpha_1, beta_1, alpha_2, beta_2

center = float(alpha_1 ) / float(beta_1) - float(alpha_2) / float(beta_2) 
standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5 



##########################################################
'''
    The probability density function of Gamma_1 - Gamma_2
    
    Args:
        z: A float. x1 - x2.
        alpha_1: An int. Alpha of Gamma1.
        beta_1: An int. Beta of Gamma1.
        alpha_2: An int. Alpha of Gamma2.
        beta_2: An int. Beta of Gamma2.
        
    Returns:
        The probability density of Gamma1- Gamma2 at z.
'''

start_time = time.time()
z = center
pdf = gammadifference.pdf_of_gamma_difference(z, alpha_1, beta_1, alpha_2, beta_2)
print "pdf:",pdf

print("--- %s seconds ---" % (time.time() - start_time))


##########################################################

'''
    The cumulative distribution function of Gamma_1 - Gamma_2 
    using mpmath.quad integrating pdf directly.
    
    Args:
        z: A float. x1 - x2.
        alpha_1: An int. Alpha of Gamma1.
        beta_1: An int. Beta of Gamma1.
        alpha_2: An int. Alpha of Gamma2.
        beta_2: An int. Beta of Gamma2.
        
    Returns:
        The cumulative probability of Gamma1- Gamma2 at z.
'''

start_time = time.time()
cdf = gammadifference.cdf_of_gamma_difference(z + 0.5 * standard_err, alpha_1, beta_1, alpha_2, beta_2)
print "cdf:",cdf

print("--- %s seconds ---" % (time.time() - start_time))

##########################################################
'''The inverse function of cdf of Gamma_1 - Gamma_2
    Given a cumulative probability y, compute z.
    
    Args:
        alpha_1: An int. Alpha of Gamma1.
        beta_1: An int. Beta of Gamma1.
        alpha_2: An int. Alpha of Gamma2.
        beta_2: An int. Beta of Gamma2.
        n: An int. The number of segments.
        
    Returns:
        The z corresponding to the cumulative probability y of Gamma1- Gamma2.
'''

start_time = time.time()

cdf_seg = gammadifference.cdf_segment(alpha_1, beta_1, alpha_2, beta_2, 200)
print cdf_seg

print("--- %s seconds ---" % (time.time() - start_time))
