import pandas as pd
import numpy as np
from math import log


### For a better understanding, view the Wikipedia page for Absorbing Markov chain: https://en.wikipedia.org/wiki/Absorbing_Markov_chain ###
# In our case, there are 3 absorbing states and 3^n-3 transient states for a total of m = 3^n states
n = 4


def digits_ternary(x):
    global n
    if x == 0:
        return [0]*n
    nums = []
    while x:
        x, r = divmod(x, 3)
        nums.append(r)
    while len(nums) < n:
        nums.append(0)
    nums.reverse()
    return nums


def value_digits_ternary(digits):
    global n
    result = 0
    for i in range(n):
        result += digits[i] * (3**(n-i-1))
    return result


def value_all_ones():
    global n
    result = 0
    for i in range(n):
        result += 3**i
    return result


def next_states(digits, position_1, position_2):
    result = []
    original_1 = digits[position_1]

    digits[position_1] = 0
    value = value_digits_ternary(digits)
    result.append(value)

    digits[position_1] = 1
    value = value_digits_ternary(digits)
    result.append(value)

    digits[position_1] = 2
    value = value_digits_ternary(digits)
    result.append(value)

    digits[position_1] = original_1

    digits[position_2] = 0
    value = value_digits_ternary(digits)
    result.append(value)

    digits[position_2] = 1
    value = value_digits_ternary(digits)
    result.append(value)

    digits[position_2] = 2
    value = value_digits_ternary(digits)
    result.append(value)

    return result


def transition_matrix(sequence_length: int):
    global n
    n = sequence_length
    """ Compute the transition matrix (related to our algorithm) of all ternary number of a given length
    Args:
        n (int): The given length of the ternary sequence
    Returns:
        P (np.array): The corresponding transition matrix of dimension m*m with m = 3^n
    """
    m = 3**n
    P = np.zeros(shape=(m, m))

    for x in range(m):
        digits = digits_ternary(x)
        position_1 = -1
        position_2 = -1

        if digits[n-1] > digits[0]:
            position_1 = n-1
            position_2 = 0
        else:
            for i in range(n-1):
                if digits[i] > digits[i+1]:
                    position_1 = i
                    position_2 = i+1

        if position_1 == -1 | position_2 == -1:  # Absorbing case
            P[x][x] = 1
        else:
            transitions = next_states(digits, position_1, position_2)
            for y in transitions:
                P[x][y] += 1/6
    return P


def transient_matrix(P: np.array):
    """ Extract the transient matrix from a transition matrix
    Args:
        P (np.array): A transition matrix of dimension m*m
    Returns:
        Q (np.array): The corresponding transient matrix of dimension (m-3)*(m-3)
    """
    global n
    ones = value_all_ones()
    twos = 3**n - 1

    Q = np.delete(P, twos, 0)
    Q = np.delete(Q, twos, 1)

    Q = np.delete(Q, ones, 0)
    Q = np.delete(Q, ones, 1)

    Q = np.delete(Q, 0, 0)
    Q = np.delete(Q, 0, 1)

    return Q


def fundamental_matrix(Q: np.array):
    """ Compute the fundamental array of a transient array
    Args:
        Q (np.array): A transient matrix of dimension (m-3)*(m-3)
    Returns:
        N (np.array): The corresponding fundamental matrix of dimension (m-3)*(m-3)
    """
    I = np.identity(Q.shape[0])
    N = np.linalg.inv(I-Q)
    return N


def absorption_time_vector(Q: np.array):
    """ Calculate the absorption time vector from a transient matrix
    Args:
        Q (np.array): A transient matrix of dimension (m-3)*(m-3)
    Returns:
        t (np.array): The absorption time vector of dimension (m-3)*1
    """
    I = np.identity(Q.shape[0])
    one = np.ones(Q.shape[0])
    return np.linalg.solve(I-Q, one)


def absorbing_probabilities_matrix(N: np.array, R: np.array):
    """ Calculate the absorption time vector from a
    Args:
        N (np.array): A fundamental matrix of dimension (m-3)*(m-3)
        R (np.array): A right absorbtion matrix of dimension (m-3)*3
    Returns:
        B (np.array): The corresponding absorbing probabilities matrix of dimension (m-3)*3
    """
    B = np.zero(1)
    return B


if __name__ == '__main__':
    for length in range(2, 8):
        P = transition_matrix(length)
        # df_transition = pd.DataFrame(P)
        # filepath = "transition_matrix.xlsx"
        # df_transition.to_excel(filepath)

        Q = transient_matrix(P)
        # df_transient = pd.DataFrame(Q)
        # filepath = "transient_matrix.xlsx"
        # df_transient.to_excel(filepath)

        t = absorption_time_vector(Q)
        print(length, "-->  mean:", round(t.mean(), 1), " n*n:", length *
              length, "n*n*ln(n):", round(length*length*log(length), 1))
