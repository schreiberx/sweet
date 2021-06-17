

import sys
import numpy as np
import mule_local.rexi.quadrature as quadrature



def lagrange_interpol(x, y, eval_x):
    N_x = len(x)
    N_eval = len(eval_x)

    l = np.ones((N_eval, N_x), dtype=y.dtype)
    for j in range(N_x):
        for i in range(N_x):
            if i != j:
                l[:,j] *= (eval_x - x[i])/(x[j] - x[i])

    return l.dot(y)




class rk_co:
    def __init__(self, M, collocation_type="gauss_legendre"):
        """
        _M: Number of quadrature points
        _collocation_type: e.g. "gauss_lobatto"
        """

        self._M = M
        self._collocation_type = collocation_type

        ref_x, w = quadrature.quad_coeffs(M, collocation_type)
        ref_x = (1.0 + ref_x)*0.5

        if np.isclose(ref_x[0], 0):
            print("WARNING:")
            for i in range(5):
                print("WARNING: Removing first point since this would create a singular matrix!")
            print("WARNING:")

            # Regenerate coefficients with one more coefficient
            ref_x, w = quadrature.quad_coeffs(M+1, collocation_type)
            ref_x = (1.0 + ref_x)*0.5

            ref_x = ref_x[1:]
            w = w[1:]

        self.ref_x = ref_x

        #self.A, self.b, self.c = rkanalysis.irk_collocation(ref_x)
        #return

        self.c = ref_x
        self.A = np.zeros((M, M), dtype=float)
        self.b = np.zeros(M, dtype=float)

        if 0:
            """
            Don't use this! Polyfit doesn't really work!!!
            It has very large numerical round-off errors, see assertion below
            """

            for j in range(M):
                y = np.zeros(M, dtype=float)
                y[j] = 1.0

                p = np.polyfit(self.c, y, len(y)-1)

                for i in range(M):
                    if i == j:
                        assert np.isclose(np.polyval(p, self.c[i]), 1)
                    else:
                        print(np.polyval(p, self.c[i]))
                        #
                        # This fails!!!!
                        #
                        assert np.isclose(np.polyval(p, self.c[i]), 0)

                pint = np.polyint(p)
                self.b[j] = np.polyval(pint, 1) - np.polyval(pint, 0)

                for i in range(M):
                    self.A[i,j] = np.polyval(pint, self.c[i]) - np.polyval(pint, 0)

        else:
            """
            Use simple quadrature to solve integrals
            """

            """
            Solve polynomial integrals
            """
            quad_x, quad_w = quadrature.quad_coeffs(M, "gauss_legendre")
            quad_x = (1.0 + quad_x)*0.5
            quad_w *= 0.5


            def quad(x, y, start=0, end=1):
                dx = end - start

                if np.abs(dx) < 1e-12:
                    z = np.zeros_like(x)

                else:
                    x = (x - start)/dx
                    z = lagrange_interpol(x, y, quad_x)

                return z.dot(quad_w)*dx
            
            for j in range(M):
                y = np.zeros(M, dtype=float)
                y[j] = 1.0

                self.b[j] = quad(self.c, y, 0, 1)
                for i in range(M):
                    self.A[i,j] = quad(self.c, y, 0, self.c[i])


    def integrate(self, lam, u0, dt):
        """
        u0: initial condition
        _M: Number of collocation points
        lam: ode coefficient   d_dt y(t) = lam * y(t)
        dt: time step size
        """
        k = np.linalg.solve(np.eye(self._M) - dt*lam*self.A, lam*np.array([u0 for i in range(self._M)]))
        return u0 + dt*self.b.dot(k)



    def integrate_ver2(self, lam, u0, dt):
        """
        u0: initial condition
        _M: Number of collocation points
        lam: ode coefficient   d_dt y(t) = lam * y(t)
        dt: time step size
        """
        k = np.linalg.solve((1.0/(dt*lam))*np.eye(self._M) - self.A, np.array([u0 for i in range(self._M)]))
        return u0 + self.b.dot(k)



if __name__ == "__main__":
    """
    Simple test
    """
    x = np.array([-5., -2., -1., 1., 2., 6., 8.])
    y = np.array(np.sin(x))

    z = lagrange_interpol(x, y, x)
    assert np.all(np.isclose(z, y))

    r = rk_co(16, "gauss_legendre")

