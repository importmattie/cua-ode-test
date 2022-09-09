# Constants/Tunable Fields
import numpy as np

import matplotlib.pyplot as plt

N = 5812069 # ...I didn't see a dNdt so I guess we'll assume it's constant?
beta = 0.85 # transmission rate
gamma = 0.1 # recovery rate
alpha = 4   # days from onset to infectious
TC = 0.49   # transmission control to scale beta
pS = 0.75   # percent symptomatic
hosp = 0.21 # hospitalization rate among infected
los = 10    # 10 day hospital stay

total_days = 161  # match the plot given in the example, start on day 0
days = np.arange(0,total_days)

# Modeled variables
# Store the results in 1D numpy arrays; we can initialize to zeros since we
# know only I and S have non-zero initial values and we can set those
I = np.zeros([total_days])  # I counts infectious individuals
I[0] = 1
S = np.zeros([total_days])  # S counts people susceptible to COVID
S[0] = N - I[0]
R = np.zeros([total_days])  # R counts recovered, immune people
E = np.zeros([total_days])  # E counts exposed people
A = np.zeros([total_days])  # A counts asymptomatic people
Ih = np.zeros([total_days]) # Ih counts hospitalized people

def main():
    # Modeling portion
    for t in range(0, total_days-1):
        # Simplify some of the complex/reused components of the ODEs
        beta_S_1mTC_over_N = (beta*S[t]*(1-TC))/N
        E_over_alpha = E[t]/alpha
        I_gamma_hosp = I[t]*gamma*hosp
        A_gamma = A[t]*gamma
        Ih_over_los = Ih[t]/los

        # Always difficult to decide on following the Python convention for
        # spaces between operators or not!
        dSdt = (-1*beta_S_1mTC_over_N)*(I[t]+A[t])
        dEdt = (-1*E_over_alpha)+(-1*dSdt) # the second two terms in dEdt are -1*dSdt
        dIdt = (E_over_alpha*pS)-I_gamma_hosp
        dAdt = (E_over_alpha*(1-pS))-A_gamma
        dIhdt = I_gamma_hosp-Ih_over_los
        dRdt = Ih_over_los+A_gamma

        # Calcuate the next timestep for our variables
        S[t+1] = S[t] + dSdt
        E[t+1] = E[t] + dEdt
        I[t+1] = I[t] + dIdt
        A[t+1] = A[t] + dAdt
        Ih[t+1] = Ih[t] + dIhdt
        R[t+1] = R[t] + dRdt

    #Plotting portion
    plt.plot(days, S/1e6, label='Susceptible', color='blue')
    plt.plot(days, I/1e6, label='Infected', color='red')
    plt.plot(days, R/1e6, label='Recovered', color='green')
    plt.plot(days, E/1e6, label='Exposed', color='purple')
    plt.plot(days, A/1e6, label='Asymptomatic', color='pink')
    plt.plot(days, Ih/1e6, label='Hospitalized', color='orange')
    plt.xlim([0, 160])
    plt.xlabel('Days (since initial infection)')
    plt.ylim([0, 6])
    plt.ylabel('People (in millions)')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
