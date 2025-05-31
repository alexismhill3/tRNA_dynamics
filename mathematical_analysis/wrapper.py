# Wrapper around the nsolve function from sympy 
# solves two-codon system (see eqs below) with different initial conditions

from copy import copy
from sympy import symbols, Eq, solve, nsolve
import pandas as pd

Tt1, Tt2, Tc1, Tc2, Rb1, Rb2, fopt = symbols('Tt1 Tt2 Tc1 Tc2 Rb1 Rb2 fopt')
Rt, Rb, Tt, Tc, N, L, Ksp, Kbd, Kch = symbols('Rt Rb Tt Tc N L Ksp Kbd Kch')

EQ1 = Eq((Rt-Rb1-Rb2)*Kbd*N, (Tc1*Rb1+Tc2*Rb2)*Ksp*(1/L)) 
EQ2 = Eq(Tc1*Rb1*(1-fopt), Tc2*Rb2*fopt)
EQ3 = Eq((Tt1-Tc1)*Kch, Tc1*Rb1*Ksp)
EQ4 = Eq((Tt2-Tc2)*Kch, Tc2*Rb2*Ksp)

def solve_numeric(fixed_params, 
                target, 
                vals, 
                name="param",
                time=100,
                verbose=True,
                solver="mnewton",
                guess=(15,15,15,15)):
    
    # Solve numerically for Rb1, Tc1, Rb2, and Tc2 using an mnewton solver.
    # fix_params: model parameters with fixed values
    # target: model parameter to vary (a sympy Symbol object)
    # vals: values for the target parameter
    # name: name to use for the target parameter
    # time: corresponds to seconds of simulation time
    
    fopts = [x/100 for x in range(0, 101)]
    df = None
    for idx, val in enumerate(vals):
        fixed_params[target] = val
        if verbose:
            print(fixed_params)
        solutions = []
        xvals = []
        for f in fopts:
            cp = copy(fixed_params)
            cp[fopt] = f
            Eq1 = EQ1.evalf(subs=cp)
            Eq2 = EQ2.evalf(subs=cp)
            Eq3 = EQ3.evalf(subs=cp)
            Eq4 = EQ4.evalf(subs=cp)
            try:
                solution = nsolve([Eq1, Eq2, Eq3, Eq4], 
                                  [Rb1, Rb2, Tc1, Tc2],
                                  guess,
                                  dict=True, 
                                  solver=solver)
                solutions.append(solution)
                xvals.append(f)
            except Exception:
                continue
        rb1 = [float(x[0][Rb1]) for x in solutions]
        rb2 = [float(x[0][Rb2]) for x in solutions]
        tc1 = [float(x[0][Tc1]) for x in solutions]
        tc2 = [float(x[0][Tc2]) for x in solutions]
        tmp = pd.DataFrame({"fopt": xvals,
                             "rb1": rb1,
                             "rb2": rb2,
                             "tc1": tc1,
                             "tc2": tc2})
        tmp[name] = val
        tmp["tu1"] = float(fixed_params[Tt1]) - tmp["tc1"]
        tmp["tu2"] = float(fixed_params[Tt2]) - tmp["tc2"]
        # formula for calculating total protein from Tu1
        tmp["protein"] = (time*float(fixed_params[Kch])*tmp["tu1"]) / (tmp["fopt"] * float(fixed_params[L]))
        tmp["Rb"] = tmp["rb1"] + tmp["rb2"]
        tmp["Tc"] = tmp["tc1"] + tmp["tc2"]
        
        if df is not None:
            df = df.append(tmp, ignore_index=True)
        else:
            df = tmp
    
    return df


if __name__ == "__main__":
    # test that the module is working correctly
    fixed_params = {N: 100, Kbd: 0.005, Rt: 500, Tt1: 1250, Tt2: 1250, Ksp: 0.05,L: 300}
    vals = [100]
    df = param_sweep(fixed_params, Kch, vals, name="kch")
    print(df)
