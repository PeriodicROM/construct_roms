# construct_roms
Construct reduced order models (ROMs) for Rayleigh-Bénard convection. Each ROM is a system or ordinatry differential equations designed to approximate the following governing  equations:
$$
    \dot{\psi}_{mn} &= -\sigma \rho_{mn} \psi_{mn} + (-1)^{m+n} \sigma \R \frac{mk}{\rho_{mn}} \theta_{mn} + \frac{k}{\rho_{mn}} Q_{mn}, \label{eq:psiTrunc} \\
     \dot{\theta}_{mn} &= -\rho_{mn} \theta_{mn} + (-1)^{m+n} (mk) \psi_{mn} + k \widetilde{Q}_{mn},
$$
