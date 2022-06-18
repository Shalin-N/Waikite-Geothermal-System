import math
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize as op


def interp(t0, t1, dt, flag=True):
    """ Return the interpolated values of pressure for a given range of time.

            Parameters:
            -----------
            t0 : float
                Initial time.
            t1 : float
                Final time.
            dt : float
                Time step size.
            flag : bool
                falg to use data approximation

            Returns:
            --------
            t : array-like
                array of times
            pr : array-like
                array of interpolated and converted pressure values

            Notes:
            ------
            for the specfic case where dt = 1, we have approximated the uncertain data from 1950-1984 using a circle
            function, this estimation only works for a specfic step size hence the checks.
    """
    n = int(np.ceil((t1 - t0) / dt))    # number of steps
    t = t0 + np.arange(n + 1) * dt      # time array

    # read and interpolate water level
    time, water_level = np.genfromtxt('gr_p.txt', delimiter=',', skip_header=1).T  # water level (gr_p data)
    water_interp = np.interp(t, time, water_level)

    # Calculation of water lvl from 1850 - 1875
    # Ratouis2017, figure 16, gives us historical water data
    # Can be approximated to a circle function (x + 2.0211)²  +  (y - 11.518)²  =  8.1733e+4
    # Calculated using http://www.1728.org/circle2.htm, parameters {0, 297.4; 30, 296.9; 69,295} respectively
    #  (x + 2.0211)²  +  (y - 11.518)²  =  8.1733e+4  for 1 step?
    if (dt == 1 and flag):
        for i in range(0, 34):
            water_interp[i] = math.sqrt(8.1733e+4 - (i + 2.0211) ** 2) + 11.518

    # Conversion of water level to pressure
    pr = ((water_interp - 297.4) * 997 * 9.81) + 5000

    return t, pr


def interpolate_q_total(t):
    """ interplate two extraction rates to find the total extraction rate, q.

        Parameters:
        -----------
        t : float
            Independent variable.

        Returns:
        --------
        q : float
            total extraction rate.

        Notes:
        ------
        tq1 is the total extraction rate over the rotrua region including the rhyolite formation
        tq2 is only the rhyolite formation

        instead of adding a new term to our ode we opted to simply minus the reinjection rate from the extraction rate
    """
    # read extraction rate data
    tq1, pr1 = np.genfromtxt('gr_q1.txt', delimiter=',', skip_header=1).T  # production rate 1 (gr_q1 data)
    tq2, pr2 = np.genfromtxt('gr_q2.txt', delimiter=',', skip_header=1).T  # production rate 2 (gr_q2 data)

    # interpolate extraction data
    ex1 = np.interp(t, tq1, pr1)
    ex2 = np.interp(t, tq2, pr2)

    # calculation of reinjection rate of water
    ex_final = ex1
    ex_final[34:] = ex_final[34:] - 1500  # 1500    @1985
    ex_final[41:] = ex_final[41:] - 3800  # 5300    @1992
    ex_final[50:] = ex_final[50:] - 2200  # 7500    @2001

    return ex_final / 86.4


def analytical(t, q, a, b, c, p0):
    """ Return the pressure value for a given time

    Args:
        t : float
            Independent variable.
        ap : float
            extraction strength parameter.
        bp : float
            recharge strength parameter.
        cp : float
            slow drainage strength parameter.
        p0 : float
            hydrostatic pressure value of recharge source.

    Returns:
        P : float
            Pressure value for given time and parameters

    Notes:
            This function assumes a constant extraction rate therefore the slow drainage term is equal to 0
    """
    return p0 - (a * q) / b * (1 - math.exp(-b * t))


def ode_model(t, pr, q, dqdt, a, b, c, p0):
    """ Return the derivative dP/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable.
        pr : float
            Dependent variable.
        q : float
            relative extraction rate.
        dqdt : float
            rate of change of relative extraction rate.
        a : float
            extraction strength parameter.
        b : float
            recharge strength parameter.
        c : float
            slow drainage strength parameter.
        p0 : float
            hydrostatic pressure value of recharge source.

        Returns:
        --------
        dpdt : float
            Derivative of dependent variable with respect to independent variable.

        Notes:
        ------
        q = {qtotal,qrhyolite,qnotrhyolite}
        q is found by using interpolate_q_total(t):
        - 1. interpolating extraction rate 1 and extraction rate 2 to independent variables values that corresponds to input variable t.
        - 2. summing the two interpolated data.

        Examples:
        ---------
        >>> pressure_ode_model(0, 1, 2, 3, 4, 5, 6, 7)
        = -12
    """

    # the first derivative returns dP/dt = -a*q-b*(P-P0)-c*dqdt where all the parameters are provided as inputs
    return -a * q - b * (pr - p0) - c * dqdt


def solve_ode(f, t0, t1, dt, x0, indicator, pars):
    """ Solve the pressure ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dxdt given variable and parameter inputs.
        t0 : float
            Time at solution start.
        t1 : float
            Time at solution end.
        dt : float
            Time step length.
        x0 : float
            Initial value of solution.
        indicator: string
            string that describes the future operation of production. PLEASE REFER TO NOTES TO MORE INFORMATION.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        ts : array-like
            Independent variable solution vector.
        ys : array-like
            Dependent variable solution vector.

        Notes:
        ------
        ODE should be solved using the Improved Euler Method.

        Function q(t) should be hard coded within this method. Create duplicates of
        solve_ode for models with different q(t).

        Assume that ODE function f takes the following inputs, in order:
            1. independent variable
            2. dependent variable
            3. forcing term, q
            4. all other parameters

        if indicator is:
            'SAME' => No extrapolation is wanted or production is maintained from year 2014; there is no change in production rate from 2014
            'STOP' => production is stopped from year 2014; production rate = 0 from q[65]
            'DOUBLE' => production is doubled from year 2014
            'HALF' => production is halved from year 2014
    """
    n = int(np.ceil((t1 - t0) / dt))    # number of steps
    ts = t0 + np.arange(n + 1) * dt     # time array
    ys = 0. * ts                        # initialise solution array
    ys[0] = x0                          # set initial value of solution array

    # total extraction rate and gradient arrays
    q = interpolate_q_total(ts)
    dqdt = np.gradient(q)

    # different value of q for extrapolation scenarios
    # written so the change is non linear and not instant
    if indicator == 'SAME':
        temp = 0

    elif indicator == 'STOP':
        a = q[64] / 36
        for i in range(65, 101):
            q[i] = q[64] - a * (i - 65)

    elif indicator == 'DOUBLE':
        a = q[64] / 36
        for i in range(65, 101):
            q[i] = q[64] + a * (i - 65)

    elif indicator == 'HALF':
        a = q[64] / 72
        for i in range(65, 101):
            q[i] = q[64] - a * (i - 65)

    # calculate solution values using Improved Euler
    for i in range(n):
        fk = f(ts[i], ys[i], q[i], dqdt[i], *pars)
        fk1 = f(ts[i] + dt, ys[i] + dt * fk, q[i], dqdt[i], *pars)
        ys[i + 1] = ys[i] + dt * ((fk + fk1) / 2)

    # Return both arrays contained calculated values
    return ts, ys


def helper(t, dt, x0, indicator, a, b, c, p0):
    """ A helper method for curve_fit.

        Parameters:
        -----------
        t : float
            Independent variable.
        dt : float
            Time step length.
        x0 : float
            Initial value of the dependent variable.
        indicator: string
            String that describes the future operation of production. 
        a : float
            Extraction strength parameter.
        b : float
            Recharge strength parameter.
        c : float
            Slow drainage strength parameter.
        p0 : float
            hydrostatic pressure of recharge source.

        Returns:
        --------
        dxdt : float
            Derivative of dependent variable with respect to independent variable.

        Notes:
        ------
            this method when used in combination with the lambda function allows us to make quick change to our other
            parameters and resolve for the a, b, c.

            It also break up the time array into the correct format.

    """

    # correct time format
    t0 = t[0]
    t1 = t[-1]

    return solve_ode(ode_model, t0, t1, dt, x0, indicator, pars=[a, b, c, p0])[1]


def fit(t, wp, dt, x0, p0, flag=False):
    """ A helper method for curve_fit.

            Parameters:
            -----------
            t : float
                Independent variable.
            wp : float
                Dependent variable.
            dt : float
                Time step length.
            x0 : float
                Initial pressure value.
            p0 : float
                Hydrostatic pressure value of recharge source.
            flag : bool
                indicator whether to use slow drainage term or not.

            Returns:
            --------
            ap : float
                Extraction strength parameter.
            bp : float
                Recharge strength parameter.
            cp : float
                Slow drainage strength parameter.

            Notes:
            ------

            the noise of the data is the how rainfall effects to the water levels which in turn effects to
            pressure readings. The covariance array from curve_fit can be used to model uncertainty

    """

    # implementation of uncertainty generated by the varying annual rainfall
    tqr, prr = np.genfromtxt('gr_rainfall.txt', delimiter=',', skip_header=1).T  # rainfall data
    rf = np.interp(t, tqr, prr)
    sigma = rf*997*9.81/1000

    # check for slow drainage
    if (flag):
        para, cov = op.curve_fit(lambda t, a, b, c: helper(t, dt, x0, 'SAME', a, b, c, p0),
                                 xdata=t,
                                 ydata=wp, p0=[0.15, 0.12, 0.6],
                                 bounds=((0, 0, -np.inf), (np.inf, np.inf, np.inf)),
                                 sigma=sigma)
    else:
        c = 0
        para, cov = op.curve_fit(lambda t, a, b: helper(t, dt, x0, 'SAME', a, b, c, p0),
                                 xdata=t,
                                 ydata=wp, p0=[0.15, 0.12],
                                 bounds=((0, 0), (np.inf, np.inf)),
                                 sigma=sigma)

    return para, cov
