
---

### **1. Maximum Likelihood Sequence Estimator (MLSE)**

This is the optimal receiver for channels with Inter-Symbol Interference (ISI). Instead of making decisions one symbol at a time, it looks at the entire received sequence to find the most probable transmitted path.

#### **Mathematical Formulation**
The estimator searches for the sequence $\{I_n\}$ that minimizes the Euclidean distance between the actual received signal and the reconstructed signal:

$$\hat{I}_n = \arg \min_{I_n} \int_{-\infty}^{\infty} \left| r(t) - \sum_{n} I_n h(t-nT) \right|^2 dt$$


**Symbol Definitions:**
* **$\hat{I}_n$**: The estimated sequence of transmitted symbols (the final output).
* **$\arg \min_{I_n}$**: The operator that finds the specific sequence $I_n$ that minimizes the expression.
* **$I_n$**: A candidate sequence of transmitted symbols (a hypothesis being tested).
* **$\int_{-\infty}^{\infty} ... dt$**: Integral over continuous time $t$, representing the total energy of the error.
* **$|\cdot|^2$**: The squared magnitude of the complex value inside.
* **$r(t)$**: The actual continuous-time signal received from the channel.
* **$\sum_{n}$**: Summation over all symbol indices $n$.
* **$h(t)$**: The composite channel impulse response (combining the transmitter filter, physical multipath channel, and receiver filter).
* **$t$**: Continuous time variable.
* **$T$**: The symbol period (time duration of one symbol).

#### **Implementation**
* **Viterbi Algorithm:** Because the complexity of testing every possible sequence is too high, the Viterbi algorithm is used to efficiently find the optimal path through a trellis diagram. The states of the trellis are defined by the channel memory length $L$.

#### **Trade-offs**
* **Advantages:** It is the **optimal** receiver; it minimizes the probability of sequence error.
* **Disadvantages:** **High Complexity**. The computational cost grows exponentially with the channel memory $L$ (complexity $\propto M^L$, where $M$ is the constellation size). It is impractical for long channels.

---

### **2. Linear Equalizers (Time Domain / Burst Model)**

To reduce complexity, these equalizers use linear filters. They rely on a discrete-time **matrix model** of the transmission burst.

**The System Matrix Model:**
$$\underline{r} = \underline{\underline{H}} \cdot \underline{I} + \underline{z}$$


**Symbol Definitions:**
* **$\underline{r}$**: The received signal vector (containing samples $r_0, r_1, ...$).
* **$\underline{\underline{H}}$**: The channel convolution matrix (Toeplitz structure representing the multipath delay).
* **$\cdot$**: Matrix multiplication operator.
* **$\underline{I}$**: The vector of transmitted symbols ($I_0, I_1, ...$).
* **$\underline{z}$**: The noise vector (Gaussian noise samples).

---

#### **A. Matched Filter (MF)**

The Matched Filter maximizes the Signal-to-Noise Ratio (SNR) at the sampling instant but ignores interference.

**Equation:**
$$\underline{\hat{I}}_{MF} = \underline{\underline{H}}^H \cdot \underline{r}$$


**Symbol Definitions:**
* **$\underline{\hat{I}}_{MF}$**: The estimated symbol vector output by the Matched Filter.
* **$\underline{\underline{H}}^H$**: The Hermitian transpose (conjugate transpose) of the channel matrix $\underline{\underline{H}}$.
* **$\underline{r}$**: The received signal vector.

**Trade-offs:**
* **Advantages:** Maximizes SNR in the presence of noise.
* **Disadvantages:** Does not remove Inter-Symbol Interference (ISI). Performance is poor in multipath environments.

---

#### **B. Zero-Forcing (ZF) Equalizer**

The ZF equalizer completely inverts the channel matrix to force the ISI to zero.

**Filter Equation:**
$$\underline{\underline{F}}_{ZF} = (\underline{\underline{H}}^H \cdot \underline{\underline{H}})^{-1} \cdot \underline{\underline{H}}^H$$


**Estimation Equation:**
$$\underline{\hat{I}}_{ZF} = \underline{\underline{F}}_{ZF} \cdot \underline{r} = \underline{I} + (\underline{\underline{H}}^H \cdot \underline{\underline{H}})^{-1} \underline{\underline{H}}^H \underline{z}$$


**Symbol Definitions:**
* **$\underline{\underline{F}}_{ZF}$**: The Zero-Forcing equalizer matrix.
* **$(\cdot)^{-1}$**: Matrix inversion.
* **$\underline{\underline{H}}$**: The channel matrix.
* **$\underline{\underline{H}}^H$**: The Hermitian transpose of the channel matrix.
* **$\underline{\hat{I}}_{ZF}$**: The estimated symbol vector.
* **$\underline{I}$**: The original transmitted symbols (perfectly recovered).
* **$\underline{z}$**: The noise vector.

**Trade-offs:**
* **Advantages:** Completely eliminates ISI ($\underline{\hat{I}}$ is unbiased).
* **Disadvantages:** **Noise Enhancement**. If the channel has deep fades (small eigenvalues), the inverse matrix contains very large numbers, which amplifies the noise term ($\underline{z}$) significantly.

---

#### **C. Minimum Mean Square Error (MMSE) Equalizer**

The MMSE equalizer balances the removal of ISI with the suppression of noise. It minimizes the total error variance.

**Equation:**
$$\underline{\underline{F}}_{MMSE} = \left( \frac{2N_0}{\sigma_I^2}\underline{\underline{I}}_N + \underline{\underline{H}}^H \cdot \underline{\underline{H}} \right)^{-1} \cdot \underline{\underline{H}}^H$$


**Symbol Definitions:**
* **$\underline{\underline{F}}_{MMSE}$**: The MMSE equalizer matrix.
* **$2N_0$**: The variance (power) of the noise samples.
* **$\sigma_I^2$**: The variance (average power) of the transmitted symbols.
* **$\underline{\underline{I}}_N$**: The identity matrix of size $N \times N$.
* **$\underline{\underline{H}}$**: The channel matrix.
* **$\underline{\underline{H}}^H$**: The Hermitian transpose of the channel matrix.
* **$(\cdot)^{-1}$**: Matrix inversion.

**Trade-offs:**
* **Advantages:** Generally outperforms ZF because it prevents noise enhancement. At high SNR, it behaves like ZF; at low SNR, it behaves like a Matched Filter.
* **Disadvantages:** Does not perfectly remove ISI (leaves some residual interference).

---

### **3. Decision-Feedback Equalizer (DFE)**

The DFE is a non-linear equalizer. It uses a linear filter for the current symbol and subtracts the interference from **previously detected** symbols.

**Equation:**
$$\underline{\hat{I}}^{DFE} = \underline{\hat{I}}' - (\underline{\underline{L}}^{-1} - \underline{\underline{I}}_N) \cdot \tilde{\underline{I}}$$


**Symbol Definitions:**
* **$\underline{\hat{I}}^{DFE}$**: The final estimated symbol vector from the DFE.
* **$\underline{\hat{I}}'$**: An intermediate linear estimate (typically from an MMSE filter before feedback).
* **$\underline{\underline{L}}^{-1}$**: The inverse of the lower triangular matrix $\underline{\underline{L}}$ (derived from the Cholesky decomposition of the error matrix).
* **$\underline{\underline{I}}_N$**: The identity matrix of size $N \times N$.
* **$\tilde{\underline{I}}$**: The vector of **previously detected** symbols (hard decisions, e.g., rounded to the nearest QAM point).

**Trade-offs:**
* **Advantages:** Outperforms linear equalizers by canceling known interference without amplifying noise.
* **Disadvantages:** **Error Propagation**. If the previous decision $\tilde{\underline{I}}$ is incorrect, the feedback loop adds more error instead of fixing it, corrupting future symbols.

---

### **4. OFDM Equalizers (Frequency Domain)**

In OFDM, the Cyclic Prefix (CP) and FFT transform the channel matrix into a **diagonal** matrix. This decouples the system into independent subcarriers, allowing for simple "One-Tap" equalization.

**Model per subcarrier $q$:**
$$r_q^F = h_q^F \cdot I_q^F + z_q^F$$


#### **A. One-Tap Zero-Forcing (ZF)**

**Equation:**
$$\hat{I}_q^F = \frac{r_q^F}{h_q^F}$$


**Symbol Definitions:**
* **$\hat{I}_q^F$**: The estimated symbol on subcarrier $q$ (Frequency Domain).
* **$r_q^F$**: The received signal on subcarrier $q$ (after FFT).
* **$h_q^F$**: The channel frequency response (complex coefficient) at subcarrier $q$.

**Trade-offs:**
* **Advantages:** Extremely simple implementation (one division). Perfectly removes channel distortion.
* **Disadvantages:** Massive noise amplification if $|h_q^F| \approx 0$ (deep fade).

#### **B. One-Tap MMSE**

**Equation:**
$$\hat{I}_q^F = \frac{(h_q^F)^*}{|h_q^F|^2 + \frac{2N_0}{\sigma_I^2}} \cdot r_q^F$$


**Symbol Definitions:**
* **$\hat{I}_q^F$**: The estimated symbol on subcarrier $q$.
* **$(h_q^F)^*$**: The complex conjugate of the channel coefficient at subcarrier $q$.
* **$r_q^F$**: The received signal on subcarrier $q$.
* **$|h_q^F|^2$**: The squared magnitude (power) of the channel gain at subcarrier $q$.
* **$2N_0$**: The noise variance (power).
* **$\sigma_I^2$**: The signal variance (transmitted power).

**Trade-offs:**
* **Advantages:** Solves the noise enhancement problem of ZF by adding a regularization term ($\frac{2N_0}{\sigma_I^2}$) in the denominator.
* **Disadvantages:** Requires knowledge of the SNR ($2N_0$ and $\sigma_I^2$) to compute accurately.