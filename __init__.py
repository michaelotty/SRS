from scipy import signal
import numpy as np
import matplotlib.pyplot as plt


class SRS:
    def __init__(self, waveform, fs, q=10, f=12):
        """Shock Response Spectrum"""

        iunit = 1  # input(' Enter acceleration unit:   1= G   2= m/sec**2  ')

        t = np.linspace(0, len(waveform)/fs, len(waveform))
        y = waveform

        tmx = t.max()
        tmi = t.min()
        n = len(y)
            
        dt = (tmx-tmi) / (n-1)
        sr = 1 / dt

        fn = np.array([1.0])  # input(' Enter the starting frequency (Hz)  ')
        if fn[0] > sr / 30:
            fn[0] = sr / 30

        idamp = 2  # input(' Enter damping format:  1= damping ratio   2= Q  ')

        if idamp == 1:
            damp = 0.05  # input(' Enter damping ratio (typically 0.05) ')
        else:
            Q = 10  # input(' Enter the amplification factor (typically Q=10) ')
            damp = 1.0 / (2.0 * Q)

        # 1=Kelly-Richman
        ialgorithm = 1

        tmax = (tmx-tmi) + 1./fn[0]
        limit = round( tmax/dt )
        n = limit
        yy = np.zeros(int(limit))
        for i, pt in enumerate(y):
            yy[i] = pt

        print(' Calculating response..... ')

        # SRS engine
        x_pos = np.array([])
        x_neg = np.array([])
	
        for j in range(1000):
            omega = 2 * np.pi * fn[j]
            omegad = omega * np.sqrt(1 - (damp**2))
            cosd = np.cos(omegad * dt)
            sind = np.sin(omegad*dt)
            domegadt = damp * omega * dt

            if ialgorithm == 1:
                a1 = 2 * np.exp(-domegadt) * cosd
                a2 = -np.exp(-2*domegadt)
                b1 = 2 * domegadt
                b2 = omega * dt * np.exp(-domegadt)
                b2 = b2 * ((omega / omegad) * (1 - 2 * (damp**2)) * sind - 2 * damp * cosd)
                b3 = 0
            else:
                E = np.exp(-damp * omega * dt)
                K = omegad * dt
                C = E * np.cos(K)
                S = E * np.sin(K)
                Sp = S / K

                a1 = 2 * C
                a2 = -E**2
                b1 = 1.0 - Sp
                b2 = 2.0 * (Sp - C)
                b3 = E**2 - Sp

            forward = [b1, b2,  b3]
            back =    [1, -a1, -a2]

            resp = signal.lfilter(forward,back,yy)

            x_pos = np.append(x_pos, resp.max())
            x_neg = np.append(x_neg, resp.min())

            jnum = j
            if  fn[j] > sr / 8.0:
                return

            fn = np.r_[fn, fn[0] * (2.0 ** (j * (1.0 / 12.0)))]

        print(' plotting output..... ')

        #  Find limits for plot
        srs_max = x_pos.max()
        if np.abs(x_neg).max() > srs_max:
            srs_max = np.abs(x_neg ).max()

        srs_min = x_pos.min()
        if np.abs(x_neg).min() < srs_min:
            srs_min = np.abs(x_neg).min()
          
        plt.figure(1)
        plt.plot(fn, x_pos, fn, np.abs(x_neg), '-.')

        if iunit==1:
            plt.ylabel('Peak Accel (G)')
        else:
            plt.ylabel('Peak Accel (m/sec**2)')

        plt.xlabel('Natural Frequency (Hz)')
        Q = 1.0 / (2.0 * damp)
        out5 = f'Acceleration Shock Response Spectrum Q={Q}'
        plt.title(out5)
        plt.grid(true)

        # set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log')
        # leg ('positive','negative',2)

        ymax= 10 ** (np.round(np.log10(srs_max)+0.8))
        ymin= 10 ** (np.round(np.log10(srs_min)-0.6))

        fmax = fn.max()
        fmin = fmax / 10.0

        fmax = 10 ** (np.round(np.log10(fmax)+0.5))

        if  fn(1) >= 0.1:
            fmin=0.1

        if  fn(1) >= 1:
            fmin=1

        if  fn(1) >= 10:
            fmin=10

        if  fn(1) >= 100:
            fmin=100

        plt.axis([fmin, fmax, ymin, ymax])

        print('Plot pseudo velocity?')
        vchoice = int(input('1=yes   2=no'))

        if vchoice == 1:
            plt.figure(2)

        # Convert to pseudo velocity
        #
        for j in range(1, jnum):  
            if iunit==1:
               x_pos[j] = 386.0 * x_pos[j] / (2.0 * np.pi * fn[j])
               x_neg[j] = 386.0 * x_neg[j] / (2.0 * np.pi * fn[j])
            else:
               x_pos[j] = x_pos[j] / (2.0 * np.pi * fn[j])
               x_neg[j] = x_neg[j] / (2.0 * np.pi * fn[j])

        srs_max = x_pos.max()
        if np.abs(x_neg).max() > srs_max:
            srs_max = np.abs(x_neg).max()

        srs_min = x_pos.min()
        if np.abs(x_neg).min() < srs_min:
            srs_min = np.abs(x_neg).min()

        plt.plot(fn, x_pos, fn, np.abs(x_neg), '-.')

        if iunit==1:
            plt.ylabel('Velocity (in/sec)')
        else:
            plt.ylabel('Velocity (m/sec)')   

        plt.xlabel('Natural Frequency (Hz)')
        Q = 1.0 / (2.0 * damp)
        out5 = f'Pseudo Velocity Shock Response Spectrum Q={Q}'
        plt.title(out5)
        plt.grid(True)

        # set(gca,'MinorGridLineStyle','none','GridLineStyle',':','XScale','log','YScale','log')
        # leg ('positive','negative',2)

        ymax= 10 ** (np.round(np.log10(srs_max)+0.8))
        ymin= 10 ** (np.round(np.log10(srs_min)-0.6))

        plt.axis([fmin, fmax, ymin, ymax])
        plt.show()


    def plot(self):
        pass

