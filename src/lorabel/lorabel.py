import builtins
import time
import math
from typing import List, NoReturn, Union

from . import types as ty

# grab all core routines (zeros, absolute, matmul, multiply…)  
from unipy import *  

# and pull in linear‐algebra routines  
import unipy.linalg as linalg


class LORABEL:
    """
    LOw-RAnk Background ELimination (LORABEL) for stable principal component pursuit

    Decomposes input data A into:
        A = B + C + D
    where
        B = low-rank background
        C = sparse point sources
        D = residual noise

    Parameters
    ----------
    maxiter : int, default=50
        Maximal number of iterations

    sv_factor : float [0-1], default=0.05
        Singular value factor to add in order to use truncated svd (see _update_svd)

    lambda_factor : float, default=1
        Setting factor to be multiplied with standard lambda: 1/sqrt(n) (see [1])

    lambda_opt : float, default=None
        Setting lambda automatically to preferred value

    verbose : bool, default=True
        Intermediate results printing

    Attributes
    ----------
    b :
        Output reconstructed low-rank array

    c :
        Output sparse error array

    d :
        Output residual error array

    run : array
        Run principal component pursuit

    References
    ----------
    [1] Zhou, Zihan, et al. "Stable principal component pursuit."
    2010 IEEE international symposium on information theory. IEEE, 2010.

    [2] Aybat, Necdet S., and Garud Iyengar. "An alternating direction
    method with increasing penalty for stable principal component pursuit."
    Computational Optimization and Applications 61.3 (2015): 635-668.

    """
    
    def __init__(
        self,
        maxiter=50,
        stdev=1,
        dirac=None,
        sv_factor=0.05,
        lambda_factor=1,
        lambda_opt=None,
        stopping=1e-3,
        verbose=True,
        **kwargs,
    ) -> NoReturn:
        self.maxiter = maxiter
        self.sv_factor = sv_factor
        self.stdev = stdev
        self.dirac = dirac
        self.lambda_factor = lambda_factor
        self.lambda_opt = lambda_opt
        self.stopping = stopping
        self.verbose = verbose
        self.kwargs = kwargs

        pass
    
    @property
    def b(self) -> ty.Array:
        """Reconstruct b (low-rank component) from svd components

        Returns
        -------
            array : reconstructed low-rank component

        """
        return matmul(multiply(self._b[0], self._b[1]), self._b[2])
    
    @property
    def b_old(self) -> ty.Array:
        """Reconstruct b_old (low-rank component) from svd components

        Returns
        -------
            array : reconstructed low-rank component

        """
        return matmul(multiply(self._b_old[0], self._b_old[1]), self._b_old[2])
    
    @property
    def d(self) -> ty.Array:
        """
        Reconstruct optimization residual (e)

        """
        return self.a - self.b - self.c
    
    @b.setter
    def b(self, a: tuple[ty.Array]) -> None:
        """Sets b from tuple

        Parameters
        ----------
        a : tuple[array]
            Tuple containing matrices of svd to set b

        """
        self._b[0] = a[0]
        self._b[1] = a[1]
        self._b[2] = a[2]

        return None
    
    def _initialize_lambda(self) -> None:
        """Initialize lambda depending on its setting, lambda_opt is overwriting"""
        if self.lambda_opt is None:
            self.lambda_opt = (
                self.lambda_factor * 1 / math.sqrt(builtins.max(self._m, self._n))
            )
        else:
            self.lambda_factor = (
                1 / math.sqrt(builtins.max(self._m, self._n)) / self.lambda_opt
            )

        return None
    
    def _initialize_dirac(self):
        """Initialize dirac depending on its setting, dirac_opt is overwriting"""
        if self.dirac is None:
            self.dirac = (
                math.sqrt(
                    builtins.min(self._m, self._n)
                    + sqrt(8 * builtins.min(self._m, self._n))
                )
                * self.stdev
            )
        else:
            self.stdev = (
                1
                / math.sqrt(
                    builtins.min(self._m, self._n)
                    + math.sqrt(8 * builtins.min(self._m, self._n))
                )
                * self.dirac
            )

        return None
    
    def _initialize_b(self) -> None:
        """Initialize b"""
        self._b = [0, 0, 0]
        self._b[0] = zeros(
            (self._m, 1), atype=find_package(self.a)[1], dtype=self.datatype
        )
        self._b[1] = zeros((1, 1), atype=find_package(self.a)[1], dtype=self.datatype)
        self._b[2] = zeros(
            (1, self._n), atype=find_package(self.a)[1], dtype=self.datatype
        )

        return None
    
    def _initialize_b_old(self) -> None:
        """Initialize b_old"""
        self._b_old = [0, 0, 0]
        self._b_old[0] = zeros(
            (self._m, 1), atype=find_package(self.a)[1], dtype=self.datatype
        )
        self._b_old[1] = zeros(
            (1, 1), atype=find_package(self.a)[1], dtype=self.datatype
        )
        self._b_old[2] = zeros(
            (1, self._n), atype=find_package(self.a)[1], dtype=self.datatype
        )

        return None

    def _initialize_czy(self, c: float) -> None:
        """Initialize c, z, and y

        Parameters
        ----------
        c : float
            Threshold value

        """
        self.c = sign(self.a) * maximum(absolute(self.a) - c, 0)
        self._z = self.theta / (self.theta + self.mu) * (self.a - self.c)
        self._y = -self.mu * self._z

        return None

    def _initialize(self, a: ty.Array) -> None:
        """Initialize other help variables

        Parameters
        ----------
        a : array
            Input array

        """
        self.a = a
        self.datatype = a.dtype
        self._m, self._n = self.a.shape
        # Initialize variables
        self.atype = find_package(self.a)[1]
        self.c = zeros((self._m, self._n), atype=self.atype, dtype=a.dtype)
        self.c_hat = zeros((self._m, self._n), atype=self.atype, dtype=a.dtype)
        self._initialize_b()
        self._initialize_b_old()
        # Initialize Attributes
        self._a_norm_2 = linalg.norm(self.a, 2)
        self.mu = 1.25 / self._a_norm_2
        self.rho = 1.65  # Kappa in the paper
        self.mu_bar = self.mu * 1e8
        self._initialize_dirac()
        self._initialize_lambda()
        self._theta_search(self.a.flatten())
        c = (
            0
            if self.theta == 0
            else self.lambda_opt * (self.mu + self.theta) / (self.mu * self.theta)
        )
        self._initialize_czy(c)
        # Initialize Class Attributes
        self._k = 0  # Set iteration counter
        self._sv = 0  # set sv
        self.b_rank = [0]  # Rank change over iterations
        self.c_card = [0]  # Cardinality change over iterations
        self._x = {
            "a": self.a,
            "b": self.b,
            "b_old": self.b_old,
            "b_old_sv": self._b_old[1],
            "c": self.c,
            "c_hat": self.c_hat,
            "y": self._y,
            "a_2": self._a_norm_2,
            "mu": self.mu,
            "k": self._k,
            "svp": 0,
        }  # Dictionary with references for criteria
        self._tic = 0  # Initialize counter

        return None

    def _check_convergence(self) -> bool:
        """
        Default convergence: ||B - B_old||_F and ||C - C_hat||_F small relative to previous singular values.
        """
        # numerator: sqrt(||b - b_old||_F^2 + ||c - c_hat||_F^2)
        num = math.sqrt(
            linalg.norm(self.b - self.b_old, 'fro')**2
            + linalg.norm(self.c - self.c_hat, 'fro')**2
        )
        # denominator: sqrt(||b_old_sv||_2^2 + ||c_hat||_F^2) + 1
        den = math.sqrt(
            sum(self._b_old[1]**2) + linalg.norm(self.c_hat, 'fro')**2
        ) + 1
        self.crit_num = (num / den)
        return self.crit_num < self.stopping
    
    def decompose(self, a: ty.Array) -> tuple[ty.Array]:
        """Run algorithm

        Parameters
        ----------
        a : array
            Input array

        """
        self._initialize(a)  # set other variables
        self._tic = time.time()  # Start timer

        while True:
            self._update_b()  # Update L
            self._update_sv()  # Update sv
            self._theta_search(
                (absolute(self.a - self.b - (1 / self.mu) * self._y)).flatten()
            )  # Update Theta
            self._update_c_hat()  # Update c_hat
            self._update_c()  # Update c
            self._update_z()  # Update z
            self.check = self._check_convergence()
            self._update_y()  # Update Y
            self._update_mu()  # Update mu
            self._update_b_old()  # Update b_old
            self._update_parameters()  # Update Parameters
            self._update_k()  # Update k
            if self.verbose:  # Print option
                self._print_out()
                
            if self.check or self._k >= self.maxiter:
                break

        self._tic = time.time() - self._tic  # Finish timer

        return self.b, self.c, self.d
    
    def _update_y(self):
        """Update Y"""
        self._y += self.mu * (self.b - self._z)

        return None
    
    def _update_c(self):
        """Update c"""
        c = (
            0
            if self.theta == 0
            else self.lambda_opt * (self.mu + self.theta) / (self.mu * self.theta)
        )
        self.c = sign(self.a - self.b - (1 / self.mu) * self._y) * maximum(absolute(self.a - self.b - (1 / self.mu) * self._y) - c, 0)

        return None
    
    def _update_z(self):
        """Update z"""
        self._z = (self.theta / (self.mu + self.theta)) * (self.a - self.c) + (
            self.mu / (self.mu + self.theta)
        ) * (self.b + (1 / self.mu) * self._y)

        return None
    
    def _update_c_hat(self):
        """Update c_hat (estimate)"""
        self.c_hat = astype(self.c, dtype=self.c.dtype, copy=True)

        return None
    
    def _update_b(self) -> None:
        """Update b by performing truncated svd"""
        if "method" in self.kwargs:
            sv = (
                self._sv + 10 if linalg.requires_sv(self.kwargs["method"]) is True else 0
            )
        else:
            sv = 0

        self.kwargs.update({"sv": sv})
        self.kwargs.update(
            {"v0": self._b[2] if sum(self._b[1]) > 1e-2 else None}
        )

        self.b = linalg.svd(self._z - (1 / self.mu) * self._y, self.kwargs)
        self._svp = int(sum(self._b[1] > 1 / self.mu))
        if self._svp != 0:
            self._threshold(self._svp, 1 / self.mu)
        else:  # if svp is 0, just reset the svd result
            self._initialize_b()

        return None
    
    def _threshold(self, num: int = 1, value: float = 0.0) -> None:
        """Threshold b to certain value

        Parameters
        ----------
        num: int
            Truncate to this number of values

        value : float
            Threshold to this value

        """
        self._b[0] = self._b[0][:, :num]
        self._b[1] = self._b[1][:num] - value
        self._b[2] = self._b[2][:num, :]

        return None
    
    def _update_b_old(self):
        """Update b_old"""
        self._b_old[0] = astype(self._b[0], dtype=self._b[0].dtype, copy=True)
        self._b_old[1] = astype(self._b[1], dtype=self._b[1].dtype, copy=True)
        self._b_old[2] = astype(self._b[2], dtype=self._b[2].dtype, copy=True)

        return None
    
    def _update_mu(self):
        """Update mu (small update: self.mu_bar+self._k removed last part)"""
        self.mu = min(self.rho * self.mu, self.mu_bar)

        return None
    
    def _print_out(self) -> None:
        """Output iteration steps during optimization"""
        text = "#{} - Rank = {}% ({}) - Card = {}% ({}) - sv = {} - time = {}s - crit = {}".format(
            self._k,
            round(100 * self.b_rank[-1] / min(self._m, self._n), 1),
            self.b_rank[-1],
            round(100 * self.c_card[-1] / (self._m * self._n), 1),
            self.c_card[-1],
            self.kwargs['sv'],
            int(time.time() - self._tic),
            self.crit_num,
        ) 
        print(text)

        return None
    
    def _update_parameters(self) -> None:
        """Update the low-rank matrix rank and cardinality of c (each iteration)"""
        self.b_rank.append(self._svp)
        self.c_card.append(count_nonzero(self.c))

        return None
    
    def _update_k(self) -> None:
        """Update k (iteration number)"""
        self._k += 1

        return None
    
    def _update_sv(self) -> None:
        """Update (predict) the number of singular values to be calculated by svd"""
        n = min(self._m, self._n)
        if self._svp < self._sv:
            self._sv = int(min(self._svp + 1, n))
        else:
            self._sv = int(min(self._svp + round(self.sv_factor * n), n))
            
        return None
    
    def _theta_search(self, a: ty.Array) -> None:
        """Return new theta value

        Parameters
        ----------
        a : array
            Input matrix

        """
        norm_a = linalg.norm(a, 2)
        if norm_a <= self.dirac:
            self.theta = 0
        else:
            absolute(a, out=a)  # this should be changed

            # added to the algorithm to speed up the process of sorting (a[a>=self.lambda_opt/self.mu])
            sort_a = sort(a[a >= self.lambda_opt / self.mu])[::-1]
            del a
            flag = 0
            if len(sort_a) == 0:  # Fix problem if no a>=self.lambda_opt/self.mu
                sort_a = 1e-10 * ones(
                    (1,)
                )  # this should be defined properly towards gpu/cpu
                flag = 1
            if (
                norm_a * (1 - self.lambda_opt / (sort_a[0] * self.mu)) < self.dirac
                or flag == 1
            ):
                self.theta = self.mu * ((1 / self.dirac) * norm_a - 1)
            else:
                flag = 0
                qq = sort_a.shape[0]
                qq_k = linspace(1, qq - 1, qq - 1, atype=find_package(self.a)[1])
                qq_ratio = reshape(
                    sort_a[1:] - (1 / self.mu) * self.lambda_opt, (qq - 1,)
                )
                qq_tail_sum = qq_k * (qq_ratio**2)  # (qq-1)
                qq_sort_a = cumsum(sort_a**2)  # (qq,)
                qq_sum_a = norm_a**2 - qq_sort_a  # (qq,)
                qq_test_sum = sqrt(
                    qq_tail_sum
                    + ((1 - (1 / (sort_a[1:] * self.mu)) * self.lambda_opt) ** 2)
                    * qq_sum_a[1:]
                )
                if sum(qq_test_sum < self.dirac) == 0:
                    pp = qq
                else:
                    pp = int(qq + 1 - sum(qq_test_sum < self.dirac))
            
                qq_coeff = numpy.array([
                        self.dirac**2,
                        (2 * (self.dirac**2) * self.mu),
                        (
                            (self.dirac**2 - tonumpy(qq_sum_a[pp - 1])) * self.mu**2
                            - (pp * self.lambda_opt**2)
                        ),
                        (-2 * self.mu * pp * self.lambda_opt**2),
                        (-pp * (self.lambda_opt**2) * self.mu**2),
                    ])

                self.theta = roots(qq_coeff.flatten())
                ind = (absolute(imag(self.theta)) < 1e-12) * (self.theta > 0)
                self.theta = matmul(self.theta, ind).real
        return None
