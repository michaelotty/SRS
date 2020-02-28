from scipy import signal
import numpy as np
import matplotlib.pyplot as plt


class SRS:
    def __init__(self, waveform, fs, q=10, f=12):
        """Shock Response Spectrum"""

        self.fs = fs
        self.spectrum = [i for i in range(100)]

        f = 500.0
        dt = 1/fs
        phi = 0.05
        omega_n = f * 2 * np.pi
        omega_d = omega_n*np.sqrt(1-phi**2)

        a = 2 * np.exp(-phi * omega_n * dt) * np.cos(omega_d * dt)
        b = -np.exp(-2 * phi * omega_n * dt)
        c = 2 * phi * omega_n * dt
        d = omega_n * dt * np.exp(-phi * omega_n * dt)*((omega_n/omega_d * (1-2*phi**2))* np.sin(omega_d*dt)-2*phi*np.cos(omega_d*dt))
        e = 0

        num = [c, d, e]
        den = [1, -a, -b]

        filt = signal.lfilter(num, den, waveform)
        plt.loglog(filt)

        plt.show()


    def plot(self):
        pass

