"""
    Discrete Truncated Wigner Approximation for 1d Ising model with
    long range interactions and time-periodic drive

    * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    * Copyright (c) 2015 Analabha Roy (daneel@utexas.edu)
    *
    *This is free software: you can redistribute it and/or modify it under
    *the terms of version 2 of the GNU Lesser General Public License
    *as published by the Free Software Foundation.
    *Notes:
    *1. The initial state is currently hard coded to be the classical ground
    *    state
    *2. Primary references are
    *   Anatoli: Ann. Phys 325 (2010) 1790-1852
    *   Mauritz: arXiv:1209.3697
    *   Schachenmayer: arXiv:1408.4441
    * - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
from __future__ import division, print_function
__all__ = ["dtwa_ising_longrange", "reductions","redirect_stdout"]
from dtwa_ising_longrange import *
