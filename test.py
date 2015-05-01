import numpy as np
from pprint import pprint
from copy import copy

class HMM():
    def __init__(self, pi=None, T=None, E=None):
        """
        All the variable are numpy arrays.

        """
        self.pi = pi
        self.T = T
        self.E = E
        self.Tlabels = []
        self.Elabels = []
        self.states = None
        self.alphabet = None
        
    def setLabels(self, T, E):
        self.Tlabels = T
        self.Elabels = E

    def setStates(self,states):
        """
        Fills in the object's state/sample count when it can.
        """
        self.states = self.T.shape[0] # Gets the number of hidden states
        
    def setAlphabet(self, obs):
        """
        Sets alphabet given an observation sequence (or list of obs).
        """
        pass
    
    def forward(self, n, obs):
        """
        Run forward algorithm up to t = n, where t > 0.
        """
        if not self.states:
            self.setStates(states)
        alpha = np.zeros((states, samples)) # Creates entire grid    
        pass

    def backward(self, n):
        """
        Run backward algorithm up to t = n, where t < T.
        """
        pass

    def forwardbackward(self, s, n):
        """
        Returns probability of being at state s at time = n given
        observations and parameters.
        """
        pass

    def baumwelch(self, obs, threshold, log = False):
        states = self.T.shape[0] # Gets the number of hidden states
        samples = len(obs)

        A = self.T
        B = self.E
        pi = copy(self.pi)

        done = False
        while not done:
            alpha = np.zeros((states, samples)) # Creates entire grid
            #c = np.zeros(samples) # Scale factor for normalizing columns later

            # Initializing first column
            #pprint("Transposed Pi : %s" % pi.T)
            #pprint("Initial Observations : %s" % self.E[:, obs[0]])

            # Transposed initial states times respective prob of obs
            alpha[:,0] = pi.T * self.E[:,obs[0]]

            # First scale factor
            #c[0] = 1.0/np.sum(alpha[:,0])
            #alpha[:,0] = c[0] * alpha[:,0]

            for t in range(1, samples):
                # Next step/column
                alpha[:,t] = np.dot(alpha[:, t-1].T, self.T).T * self.E[:,obs[t]]
                # Calculate and modify with scale factor
                #c[t] = 1.0/np.sum(alpha[:,t])
                #alpha[:,t] = c[t] * alpha[:,t]
            
            pprint(alpha)
            beta = np.zeros((states, samples))
            beta[:,samples - 1] = 1
            #beta[:, samples - 1] = c[samples - 1] * beta[:, samples - 1]

            for t in range(len(obs) - 1, 0, -1):
                beta[:, t - 1] = np.dot(self.T, (self.E[:, obs[t]] * beta[:,t]))
                #beta[:, t - 1] = c[t-1] * beta[:, t - 1]
           
            # Gamma is our forward backward portion
            gamma = np.zeros((states,samples))
            for t in range(samples):
                for i in range(states):
                    gamma[i,t] = (alpha[i,t] * beta[i,t]) / np.dot(alpha[:,t].T,beta[:,t])

            pprint(gamma)
            # Transition probabilities
            eps = np.zeros((states, states, samples))
            for t in range(samples - 1):
                for i in range(states):
                    for j in range(states):
                        eps[i,j,t] = alpha[i,t] * self.T[i,j] * beta[j,t + 1] * self.E[j,obs[t+1]]
                        eps[i,j,t] = eps[i,j,t] / np.dot(alpha[:,t].T, beta[:,t])

            print(eps[:,:,:])
            
            return 0

    def viterbi(self, obs):
        pass

    def simulate(self, steps, log=False):
        """
        Generates an observed sequence using the current parameters.
        """
        def draw(prob_vector):
            """
            Returns the column which 'succeeded' in the multinomial roll.
            """
            return np.where(np.random.multinomial(1, prob_vector) == 1)[0][0]


        print("Initializing required arrays.")
        states = np.zeros(steps)
        states[0] = draw(self.pi)
        obs = np.zeros(steps)
        obs[0] = draw(self.E[states[0],:])
        for t in range(1 ,steps):
            states[t] = draw(self.T[states[t-1],:])
            obs[t] = draw(self.E[states[t],:])
        if log:
            print("".join([str(int(x)) for x in states]))
            print("".join([str(int(x)) for x in obs]))
        return obs, states


    def __str__(self):
        return "Pi\n%s\nT\n%s\nE\n%s" % (self.pi, self.T, self.E)


if __name__ == "__main__":
    test = HMM()
    
    # Uses test on wikipedia
    obs = ["11", "11", "11", "11", "11",
            "10", "00", "01", "11", "11"] 
    pi = np.array([.2, .8])
    a = np.array([
        [0.5, 0.5],
        [.3, .7]
        ])
    b = np.array([
        [0.3, .7],
        [0.8, .2]
        ])
    test = HMM(pi, a, b)
    test.baumwelch("00", .01, True)
