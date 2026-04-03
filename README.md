# Machine Learning the 2D Ising Model

This repository reproduces the main ideas of the paper by Carrasquilla & Melko (2017), using a simple neural network to distinguish ordered and disordered phases of the 2D Ising model.

The workflow is:

1. Generate spin configurations with a Metropolis Monte Carlo sampler written in Fortran for faster sampling.
2. Save configurations for several temperatures and lattice sizes
3. Train a neural network to classify ordered vs disordered states
4. Evaluate the neural-network output as a function of temperature
5. Estimate the critical temperature using finite-size scaling and crossing analysis

---



## Physics Background

The 2D Ising model consists of spins:

```math
s_i = \pm 1
```

on a square lattice with nearest-neighbor interactions.

The Hamiltonian is:

```math
H = -J \sum_{\langle i,j \rangle} s_i s_j
```

At low temperature, spins align and the system is ordered.
At high temperature, thermal fluctuations dominate and the system becomes disordered.

The exact critical temperature is:

```math
T_c = \frac{2}{\ln(1 + \sqrt{2})} \approx 2.269
```

---

## Monte Carlo Sampler

The Fortran code generates spin configurations using the Metropolis algorithm.

For each temperature:

* The lattice is initialized with a random hot start
* Several thermalization sweeps are performed
* Decorrelated configurations are collected
* Each configuration is saved in binary format
* Magnetization density is also computed and stored

Output files:

```text
configs_L20_T2.0000.bin
configs_index.dat
```

---

## Binary Data Format

Each `.bin` file contains:

```text
N_CONFIGS × (L × L)
```

signed 8-bit integers with values:

```text
-1, +1
```

Configurations are written in Fortran column-major order.

Example Python loading:

```python
raw = np.fromfile(filename, dtype=np.int8)
raw = raw.reshape(n_configs, L * L)
```


## Neural Network



The neural network is trained to classify:

* Ordered configurations → label = 1
* Disordered configurations → label = 0

Configurations near the critical region are excluded from training but still used during testing.


## Finite-Size Scaling

After training, the network output is averaged as a function of temperature:

```math
P(\mathrm{ordered})
```

The critical temperature can be estimated from the crossing points between different lattice sizes.

The code also produces a finite-size scaling collapse using:

```math
t = \frac{T - T_c}{T_c}
```

and:

```math
tL^{1/\nu}
```

with:

```math
\nu = 1
```

for the 2D Ising universality class.

---

## Running the Project

### 1. Generate configurations

Compile and run the Fortran sampler:

```bash
gfortran -O3 ising_sampler.f90 -o ising_sampler
./ising_sampler
```

### 2. Train on one dataset

```python
main()
```



## References

* Juan Carrasquilla and Roger G. Melko, “Machine learning phases of matter”, Nature Physics 13, 431–434 (2017)


