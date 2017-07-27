import numpy as np
from matplotlib.pyplot import figure, show

class SimpleIF(object):
    ''' Simple Integrated-and-fire model of the form:
    Cm*dV/dt + gl*(V-El) = Iapp
    
    Cm is the membrane capacitante
    gl is the leak conductance
    El is the resting membrane potential
    Ipp is the applied current
    '''

    def __init__(self, Cm, gl, El, Vinit=-70.):
        self._Cm = float(Cm) # in ms
        self._gl = float(gl) # in uS
        self._El = float(El) # in mV
        self._Vinit = float(Vinit) # in mV

        # fixed threshold
        self._Vthr = -50 # in mV

    def dVdt(self, Vmb, Iapp):
        ''' returns the instantaneous voltage change
        as a function of Cm, gl, El and Iapp
        '''
        Cm = self._Cm
        gl = self._gl
        El = self._El

        return (Iapp - gl*(Vmb - El) )/Cm

    def timecourse(self, current):
        ''' returns the time course of the voltage as
        a function of the current injected when solved
        by the implicit Euler method:
        f(x) = f(x-1) + dt*f'(x)
        
        we assume dt is the sampling interval of the current
        '''


        voltage = np.empty(len(current))

        # initial condition
        voltage[0] = self._Vinit
        dt = 0.02 # step size integration

        # Eulers solver
        for i in range(1, len(current)):
            dVdt = self.dVdt(Vmb = voltage[i-1], Iapp = current[i-1])
            voltage[i] = voltage[i-1] + dt*dVdt
            # action potential threshold
            if voltage[i] > self._Vthr:
                voltage[i-1] = 0.0 # overshooting
                voltage[i] = self._El # reset voltage

        return voltage

if __name__ == '__main__':
    # current injection
    dt = 0.02 # in ms
    voltage_scaler = 40
    T = int(X.shape[0]/dt)
    I = np.interp(np.linspace(0,T*dt,T),np.arange(len(yhat)),yhat)*voltage_scaler
    mycell = SimpleIF(Cm =4.9, gl = .16, El=-65., Vinit = -75) # spike threshold = -50 mV!!!
    voltage = mycell.timecourse(I)

    # figure
    fig1 = figure()
    ax1 = fig1.add_subplot(111)
    time = np.linspace(0, len(I)*(dt/1000), len(I)) # transform in sec

    ax1.plot(time, voltage)
    ax1.set_ylabel('Voltage (mV)')
    ax1.set_xlabel('Time (ms)')

    show()