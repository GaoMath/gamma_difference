
"""gammadifference
Is a module used to compute the pdf, cdf, inverse_cdf, cdf_segment of 
the gamma difference distribution. A gamma difference distribution can 
be defined as the difference of two gamma distribution.
"""

import mpmath
import math
import numpy
from scipy.integrate import quad

mpmath.mp.prec = 24
mpmath.mp.dps = 8
mpmath.mp.pretty = True

    
def sum_of_ln_fac(alpha_1):
    '''ln(n!)
    Is the sum of ln(i), where i start from 1 to n.
    
    Args:
        alpha_1: An int.
        
    Returns:
        A float, as the value of ln(n!)
    '''
    f_sum = 0
    for i in range(1, alpha_1):
        f_sum += math.log(i)
    return f_sum

def integrand_plus(x, alpha_1, beta_1, alpha_2, beta_2, z, c):
    standard_err = (float(alpha_1) / beta_1 ** 2 + float(alpha_2) / beta_2 ** 2 ) ** 0.5
    if 10000 * x - 10000 * z <= standard_err:
        return 0
    else:
        return math.exp( c + beta_2 * z + (alpha_1 - 1) * math.log(x) + (alpha_2 - 1) * math.log(x - z) - (beta_1 + beta_2) * x)

def integrand_minus(x, alpha_1, beta_1, alpha_2, beta_2, z, c):
    standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5
    if 10000 * x + 10000 * z <= standard_err:
        return 0
    else:
        return math.exp( c - beta_1 * z + (alpha_2 - 1) * math.log(x) + (alpha_1 - 1) * math.log(x + z) - (beta_1 + beta_2) * x)
    
def pdf_of_gamma_difference(z, alpha_1, beta_1, alpha_2, beta_2):
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
    center = float(alpha_1 ) / float(beta_1) - float(alpha_2) / float(beta_2) 
    standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5
    c = alpha_1 * math.log(beta_1) + alpha_2 * math.log(beta_2) - sum_of_ln_fac(alpha_1) - sum_of_ln_fac(alpha_2)
    if z >= 0:   
        f = quad(lambda x: integrand_plus(x, alpha_1, beta_1, alpha_2, beta_2, z, c), z, center + 100 * standard_err)[0]
    else:
        f = quad(lambda x: integrand_minus(x, alpha_1, beta_1, alpha_2, beta_2, z, c), -z, center + 100 * standard_err)[0]
    return f
    
def cdf_of_gamma_difference(z, alpha_1, beta_1, alpha_2, beta_2):
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
    center = float(alpha_1 ) / float(beta_1) - float(alpha_2) / float(beta_2) 
    standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5
    if z > center + 30 * standard_err:
        z = center + 30 * standard_err
    return mpmath.quad(lambda x: pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2), [center - 200 * standard_err, z])
    
def artificial_cdf_of_gamma_difference(z, alpha_1, beta_1, alpha_2, beta_2):
    '''
    The cumulative distribution function of Gamma_1 - Gamma_2 
    without using mpmath.
    
    Args:
        z: A float. x1 - x2.
        alpha_1: An int. Alpha of Gamma1.
        beta_1: An int. Beta of Gamma1.
        alpha_2: An int. Alpha of Gamma2.
        beta_2: An int. Beta of Gamma2.
        
    Returns:
        The cumulative probability of Gamma1- Gamma2 at z.
    '''
        
    center = float(alpha_1 ) / float(beta_1) - float(alpha_2) / float(beta_2) 
    standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5
    if z > center + 30 * standard_err:
        z = center + 30 * standard_err
    x = center - 300 * standard_err
    cdf = 0
    pdf_1 = pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2)
    delta = standard_err / 50.
    epsilon = standard_err / 10000.
    while pdf_1 * delta < 1e-10:
        x += delta
        pdf_2 = pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2)
        pdf_1 = pdf_2
        if math.fabs(pdf_2 + pdf_1) * delta < 1e-8:
            pass
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) < 0.2:
            delta = delta * 2
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) > 0.6:
            delta = delta / 2.
        if delta < epsilon:
            print 'Error'
            return None
    while x <= z:
        pdf_2 = pdf_of_gamma_difference(x + delta, alpha_1, beta_1, alpha_2, beta_2)
        cdf += (pdf_1 + pdf_2) / 2 * delta
        x += delta
        if delta < epsilon:
                print 'Error'
                return None
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) < 0.001:
            #print delta, '*2'
            delta = delta * 2
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) > 0.003:
            #print delta, '/2'
            delta = delta / 2.
        pdf_1 = pdf_2
    return cdf
    
def inverse_cdf(y, alpha_1, beta_1, alpha_2, beta_2):
    '''The inverse function of cdf of Gamma_1 - Gamma_2
    Given a cumulative probability y, compute z.
    
    Args:
        z: A float. x1 - x2.
        alpha_1: An int. Alpha of Gamma1.
        beta_1: An int. Beta of Gamma1.
        alpha_2: An int. Alpha of Gamma2.
        beta_2: An int. Beta of Gamma2.
        
    Returns:
        The z corresponding to the cumulative probability y of Gamma1- Gamma2.
    '''
    center = float(alpha_1 ) / float(beta_1) - float(alpha_2) / float(beta_2) 
    standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5
    x = center - 300 * standard_err
    cdf = 0
    pdf_1 = pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2)
    pdf_2 = pdf_1
    delta = standard_err / 30.
    epsilon = standard_err / 10000.
    if y < 1e-8 or y >= 1 - 1e-8:
        return numpy.inf
    while pdf_1 * delta  < 1e-10:
        x += delta
        pdf_2 = pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2)
        pdf_1 = pdf_2
        if math.fabs(pdf_2 + pdf_1) < 1e-8:
            pass
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) < 0.2:
            delta = delta * 2
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) > 0.5:
            delta = delta / 2.
        if delta < epsilon:
            print 'Error'
            return None
    while cdf < y:
        pdf_1 = pdf_2
        pdf_2 = pdf_of_gamma_difference(x + delta, alpha_1, beta_1, alpha_2, beta_2)
        cdf += (pdf_1 + pdf_2) / 2 * delta
        x += delta
        if delta < epsilon:
            print 'Error'
            return None
        if mpmath.fabs(pdf_2 + pdf_1) * delta < 1e-8:
            pass
        elif mpmath.fabs(pdf_2 - pdf_1) / mpmath.fabs(pdf_2 + pdf_1) < 0.001:
            delta = delta * 2.
        elif mpmath.fabs(pdf_2 - pdf_1) / mpmath.fabs(pdf_2 + pdf_1) > 0.003:
            delta = delta / 2.
    x = x - delta * (cdf - y) / ((pdf_1 + pdf_2) / 2 * delta)
    return x

def cdf_segment(alpha_1, beta_1, alpha_2, beta_2, n=200):
    '''The segment of cdf of Gamma_1 - Gamma_2
    Compute the z variables of n(=200) segments of cdf. 
    For example, cdf = 1/n, 2/n, 3/n
    
    Args:
        alpha_1: An int. Alpha of Gamma1.
        beta_1: An int. Beta of Gamma1.
        alpha_2: An int. Alpha of Gamma2.
        beta_2: An int. Beta of Gamma2.
        n: An int. The number of segments.
        
    Returns:
        A list of z. The z corresponding to the segmented cumulative probability of Gamma1- Gamma2.
    '''
    center = float(alpha_1 ) / float(beta_1) - float(alpha_2) / float(beta_2) 
    standard_err = (float(alpha_1 ) / beta_1 ** 2 + float(alpha_2 ) / beta_2 ** 2 ) ** 0.5
    x = center - 300 * standard_err
    cdf_n = []
    cdf = 0
    pdf_1 = pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2)
    pdf_2 = pdf_1
    delta = standard_err / 30.
    epsilon = standard_err / 10000.
    while pdf_1 * delta < 1e-10:
        x += delta
        pdf_2 = pdf_of_gamma_difference(x, alpha_1, beta_1, alpha_2, beta_2)
        pdf_1 = pdf_2
        if math.fabs(pdf_2 + pdf_1) * delta < 1e-8:
            pass
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) < 0.2:
            delta = delta * 2
        elif math.fabs(pdf_2 - pdf_1) / math.fabs(pdf_2 + pdf_1) > 0.5:
            delta = delta / 2.
        if delta < epsilon:
            print 'Error'
            return None
    for i in range(1,n):
        while cdf < float(i) / n:
            pdf_1 = pdf_2
            pdf_2 = pdf_of_gamma_difference(x + delta, alpha_1, beta_1, alpha_2, beta_2)
            cdf += (pdf_1 + pdf_2) / 2 * delta
            x += delta
            if delta < epsilon:
                print 'Error'
                return None
            if mpmath.fabs(pdf_2 + pdf_1) * delta < 1e-8:
                pass
            elif mpmath.fabs(pdf_2 - pdf_1) / mpmath.fabs(pdf_2 + pdf_1) < 0.001:
                delta = delta * 2
            elif mpmath.fabs(pdf_2 - pdf_1) / mpmath.fabs(pdf_2 + pdf_1) > 0.003:
                delta = delta / 2.
        cdf_n.append(float(x - delta * (cdf - float(i) / n) / ((pdf_1 + pdf_2) / 2 * delta)))
    return cdf_n
