#! /usr/bin/env python3

#
# Author: Martin Schreiber
# Email: schreiberx@gmail.com
# Date: 2017-06-18
#

import sys
import math
import mule.rexi.EFloat as ef



#
# Supported Functions to approximate
#
class Functions:

    def phiNDirect(
        self,
        n: int,
        z: float
    ):
        """
        Direct formulation of phiN functions.
        WARNING: There's a singularity close to 0!
        """
        if n == 0:
            return self.efloat.exp(z)
        elif n == 1:
            return (self.efloat.exp(z)-1)/z
        elif n == 2:
            return (self.efloat.exp(z)-1-z)/(z*z)
        elif n == 3:
            return (2*self.efloat.exp(z)-2-2*z-z*z)/(2*z*z*z)
        elif n == 4:
            return (6*self.efloat.exp(z)-6-6*z-3*z*z-z*z*z)/(6*z*z*z*z)
        elif n == 5:
            return (24*self.efloat.exp(z) -24 - 24*z - 12*z*z - 4*z*z*z - z*z*z*z)/(24*z*z*z*z*z)
        else:
            raise Exception("Not yet implemented")



    def factorial(self, N):
        """
        Helper function to support (N-1) factorials
        """
        if N < 0:
            return 1
        return math.factorial(N)



    def phiNDirectFormula(
        self,
        n: int,
        z
    ):
        """
        retval = self.factorial(n-1)*self.efloat.exp(z)
        for i in range(n):
            retval -= (self.factorial(n-1)/self.factorial(i))*self.efloat.pow(z, i)
        retval /= self.factorial(n-1)*self.efloat.pow(z, n)
        """

        # self.factorial(n-1) is cancelled out
        retval = self.efloat.exp(z)

        for i in range(n):
            retval -= self.efloat.pow(z, i)/self.factorial(i)

        retval /= self.efloat.pow(z, n)
        return retval



    def phiNRec(
        self,
        n: int,
        z
    ):
        """
        Recursive calculation of phiN functions.
        """

        if n == 0:
            return self.efloat.exp(z)

        return (self.phiN(n-1, z) - self.efloat.to(1.0)/self.efloat.to(math.factorial(n-1)))/z;



    def phiNSeries(
        self,
        n: int,
        z
    ):
        """
        It takes less than 20 iterations for cases (abs(z) < 0.5) to converge
        """

        niters = 20
        #for i in range(niters):
        #    retval += self.efloat.pow(z, i)/math.factorial(i+n)

        # Avoid repeated factorial and pow computations
        powz = self.efloat.to(1.0)
        facn = math.factorial(n)

        retval = powz/facn
        for i in range(1, niters):
            powz *= z
            facn *= (n+i)
            retval += powz/facn

        return retval



    def phiN(
        self,
        n: int,
        z
    ):
        # Use Series if z < 0.2 since this converges relatively fast
        if self.efloat.abs(z) < 0.2:
            return self.phiNSeries(n, z)

        return self.phiNRec(n, z)



    def upsNDirect(
        self,
        n: int,
        z
    ):
        if n == 1:
            return (-4-z+self.efloat.exp(z)*(4-3*z+z*z)) / (z*z*z)

        if n == 2:
            return (2+z+self.efloat.exp(z)*(-2+z)) / (z*z*z)

        if n == 3:
            return (-4-3*z-z*z+self.efloat.exp(z)*(4-z)) / (z*z*z)

        raise Exception("ups number "+str(n)+" is not supported!")



    def upsNSeries(
        self,
        n: int,
        z
    ):
        """
        It takes less than 20 iterations for cases (abs(z) < 0.5) to converge
        """

        niters = 20
        if n == 1:

            #retval = 0
            #for l in range(niters):
            #    retval += self.efloat.pow(z, l)*(l+1)*(l+1)/math.factorial(l+3)
            #return retval

            # avoid repeated pow and factorial computations
            powz = self.efloat.to(1.0)
            facn = math.factorial(3)

            retval = powz/facn
            for l in range(1, niters):
                powz *= z
                facn *= (l+3)
                retval += powz*(l+1)*(l+1)/facn

            return retval

        if n == 2:
            retval = self.efloat.to(1.0)/self.efloat.to(2.0)

            #for l in range(niters):
            #    retval += (z-2)*self.efloat.pow(z, l)/math.factorial(l+3)

            powz = self.efloat.to(1.0)
            facn = math.factorial(3)

            retval += (z-2)*powz/facn
            for l in range(1, niters):
                powz *= z
                facn *= (l+3)
                retval += (z-2)*powz/facn

            return retval

        if n == 3:
            retval = -self.efloat.to(1.0)/self.efloat.to(2.0)

            #for l in range(niters):
            #    retval += (4-z)*self.efloat.pow(z, l)/math.factorial(l+3)
            #return retval

            powz = self.efloat.to(1.0)
            facn = math.factorial(3)

            retval += (4-z)*powz/facn
            for l in range(1, niters):
                powz *= z
                facn *= (l+3)
                retval += (4-z)*powz/facn

            return retval



        raise Exception("ups number "+str(n)+" is not supported!")



    def upsN(
        self,
        n: int,
        z
    ):
        # Use Series if z < 0.2 since this converges relatively fast
        if self.efloat.abs(z) < 0.2:
            return self.upsNSeries(n, z)

        return self.upsNDirect(n, z)





    def __init__(
        self,
        function_name = "phi0",
        efloat_mode = None
    ):
        self.efloat = ef.EFloat(efloat_mode)

        self.function_name = function_name

        self.function_complex = True

        if self.efloat.floatmode == 'mpfloat':
            import mpmath as mp
            # Set numerical threshold to half of precision
            self.epsthreshold = 1e-15
        else:
            self.epsthreshold = 1e-10


        # Exponential integrator: phi0
        if self.function_name[0:3] == 'phi':

            N = int(self.function_name[3:])

            def fun(x):
                return self.phiN(N, x)

            self.eval = fun

            if self.function_complex:
                self.is_real_symmetric = True
                self.is_complex_conjugate_symmetric = True
            else:
                self.is_real_symmetric = True
                self.is_complex_conjugate_symmetric = False

        elif self.function_name[0:3] == 'ups':

            N = int(self.function_name[3:])

            if self.efloat.floatmode == 'mpfloat':
                import mpmath as mp
                # Set numerical threshold to half of precision
                self.epsthreshold = 1e-10
            else:
                self.epsthreshold = 1e-10


            if N == 1:
                #
                # Setup \upsilon_1 for EDTRK4
                # See document notes_on_time_splitting_methods.lyx
                #
                def fun(x):
                    K = x
                    if abs(x) < self.epsthreshold:
                        return self.efloat.to(1.0)/self.efloat.to(2.0*3.0)
                    else:
                        return (-self.efloat.to(4.0)-K+self.efloat.exp(K)*(self.efloat.to(4.0)-self.efloat.to(3.0)*K+K*K))/(K*K*K)
                self.eval = fun

                if self.function_complex:
                    self.is_real_symmetric = True
                    self.is_complex_conjugate_symmetric = True
                else:
                    self.is_real_symmetric = True
                    self.is_complex_conjugate_symmetric = False


            elif N == 2:
                #
                # Setup \upsilon_2 for EDTRK4
                # See document notes_on_time_splitting_methods.lyx
                #
                def fun(x):
                    K = x
                    if abs(x) < self.epsthreshold:
                        return self.efloat.to(1.0)/self.efloat.to(2.0*3.0)
                    else:
                        return (self.efloat.to(2.0)+1.0*K+self.efloat.exp(K)*(self.efloat.to(-2.0)+K))/(K*K*K)
                self.eval = fun

                if self.function_complex:
                    self.is_real_symmetric = True
                    self.is_complex_conjugate_symmetric = True
                else:
                    self.is_real_symmetric = True
                    self.is_complex_conjugate_symmetric = False


            elif N == 3:
                #
                # Setup \upsilon_3 for EDTRK4
                # See document notes_on_time_splitting_methods.lyx
                #
                def fun(x):
                    K = x
                    if abs(x) < self.epsthreshold:
                        return self.efloat.to(1.0)/self.efloat.to(2.0*3.0)
                    else:
                        return (-self.efloat.to(4.0) - 3.0*K - K*K + self.efloat.exp(K)*(self.efloat.to(4.0)-K))/(K*K*K)
                self.eval = fun

                if self.function_complex:
                    self.is_real_symmetric = True
                    self.is_complex_conjugate_symmetric = True
                else:
                    self.is_real_symmetric = True
                    self.is_complex_conjugate_symmetric = False

            else:
                print("Unknown ups function "+str(N))
                sys.exit(1)

        else:
            print("Unknown basis function "+str(self.function_name))
            sys.exit(1)


