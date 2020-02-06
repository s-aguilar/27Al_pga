import numpy as np
import scipy
from scipy import special

"""
Compute analytically the legendre parameters that minimizes chi2 a la
Bevington:

DETERMINANT Method:
-compute partial derivatives of chi2 w.r.t each 'a_i' parameter for i in [0,5]
 set equal to zero and solve for each parameter

MATRIX Method:
-solve matrix equation by computing its inverse
"""

###############################################################################
###############################################################################
###############################################################################
"""
CALCULATE P-values
"""
def pValue(numOfData,numOfPars,x2):
    # return 1 - scipy.special.gammainc( (numOfData - numOfPars)*.5 , x2*.5 )
    return scipy.special.gammaincc( (numOfData - numOfPars)*.5 , x2*.5 )

###############################################################################
###############################################################################

"""
MATRIX METHOD, INCLUDING ERRORS OF PARAMETERS
Derivation and formalism found in Bevington's Data reduction book CH 7
"""
def chi2_mat(data,data_unc,angle,order):

    x = np.cos(angle)
    weights = 1/data_unc**2

    # Needed for creating the legendre coefficient array 'legCoef' with proper
    # length
    dim = 2*order+1
    legCoef = np.zeros(dim)

    rank = order+1

    # Store array objects in the indices
    legDeriv_a = np.zeros(rank,dtype=object)

    # The derivatives w.r.t 'a_i' of the Legendre polynomial extract the value
    # of i'th term in the expansion and are stored in the 'legDeriv_a' array
    for ind in range(0,dim,2):

        # Select out only the i'th term (even terms of Legendre)
        legCoef[ind] = 1

        # Compute value of Legendre polynomial at that order (for all angles)
        legDeriv_a[ind//2] = np.polynomial.legendre.legval(x,legCoef)

        # reset legCoef so we are only Selecting the i'th term (even)
        legCoef[ind] = 0


    # beta matrix is a row vector
    beta_matrix = np.zeros(rank)

    # alpha matrix is a 'rank x rank' matrix
    alpha_matrix = np.zeros((rank,rank))

    # a_i coefficient matrix is a row vector
    aCoef_matrix = np.zeros(rank)

    # a_i error matrix is a 'rank x rank' matrix
    aCoef_err_matrix = np.zeros((rank,rank))

    a_errs = np.zeros(rank)


    # Fill out the terms in the alpha and beta matrix
    for row in range(rank):
        for col in range(rank):
            alpha_temp = np.multiply(legDeriv_a[row],legDeriv_a[col])
            alpha_matrix[row][col] = np.sum( np.multiply(alpha_temp,weights) )

        beta_temp = np.multiply(data,legDeriv_a[row])
        beta_matrix[row] = np.sum( np.multiply(beta_temp,weights) )


    # Compute the inverse of the alpha matrix which is also just the \
    # error/covariance matrix
    # alpha_inv = np.linalg.inv(alpha_matrix)
    alpha_pinv = np.linalg.pinv(alpha_matrix)


    # Matrix multiply to get the aCoef_matrix for both inverse calculations
    # aCoef_matrix = np.matmul(beta_matrix,alpha_inv)
    # print('\nResults from Bevington (Matrix method inv):\n',aCoef_matrix)
    # print('\nErr:\n',alpha_inv)

    aCoef_matrix = np.matmul(beta_matrix,alpha_pinv)
    # print('\nResults from Bevington (Matrix method pinv):\n',aCoef_matrix)
    # print('\nErr:\n',alpha_pinv)


    # Calculate chi2 per degree of freedom
    chi2 = 0
    chi2ndf = 0


    # Prepare a new list of coefficients (Read the 'convert' function
    # definition in 'legendre.py')
    new = []
    for _ in aCoef_matrix:
        new.append(_)
        new.append(0)

    # Pop the last 0 that is unnecessary
    new.pop()

    temp1 = data-np.polynomial.legendre.legval(x,new)
    temp2 = np.multiply(temp1**2,weights)
    chi2 = np.sum(temp2)


    chi2ndf = chi2/(len(data)-rank)


    """
    # Calculate errors for each coefficient from error matrix and store in list
    # Variances are diagonal terms of covariance (error matrix)
    # Errors are sqrt of that
    """
    for ind in range(rank):
        a_errs[ind] = alpha_pinv[ind][ind]  # Store the variances
    a_errs = np.sqrt(a_errs)

    # SUSPECTED WRONG, errors for parameters are diagonal terms
    # for row in range(rank):
    #     # a_var_temp = 0
    #     for col in range(rank):
    #         # Get the variances
    #         if row == col:
    #             a_var_temp += alpha_pinv[row][col]
    #         # else:
    #         #     a_var_temp += 2*alpha_pinv[row][col]
    #
    #     # Used if encountering a negative value. investigate this some more
    #     try:
    #         a_errs[row] = np.sqrt(a_var_temp)
    #     # except RuntimeWarning as err:
    #     except Exception:
    #         print('A problem was encountered:')
    #         a_errs[row] = 0

    # print('\nErr:\n',a_errs)

    # Calculate p-value
    pVal = pValue(len(angle),rank,chi2)

    results = [aCoef_matrix,a_errs,chi2,chi2ndf,pVal]

    return results


"""
Test Suite
"""
if __name__ == "__main__":
    # import sys
    # args = sys.argv[1]
    test_x = [1.1,.9,1.1,.9,1.1,.9,1.1]
    test_x_err = [.1,.1,.1,.1,.1,.1,.1]
    trialOrders = [0,1,2]

    for blah in trialOrders:
        test_results = chi2_mat(test_x,test_x_err,blah)
        print('\nFor order {0}:'.format(blah*2))
        print('coef:',test_results[0])
        print('errs:',test_results[1])
        print('chi2:',test_results[2])
        print('chiN:',test_results[3])
        print('pval:',test_results[4])
