from markov import digits_ternary
from random import randint

n = 4
m = 3**n

zeros = [0]*n
ones = [1]*n
twos = [2]*n

sample_size = 10
accumulator = 0
for sample in range(sample_size):
    x = randint(0, m)
    digits = digits_ternary(x)
    count = 0
    while digits not in [zeros, ones, twos]:
        if digits[n-1] > digits[0]:
            position_1 = n-1
            position_2 = 0
        else:
            for i in range(n-1):
                if digits[i] > digits[i+1]:
                    position_1 = i
                    position_2 = i+1
        i = randint(1, 6)
        if i == 1:
            digits[position_1] = 0
        elif i == 2:
            digits[position_1] = 1
        elif i == 3:
            digits[position_1] = 2
        elif i == 4:
            digits[position_2] = 0
        elif i == 5:
            digits[position_2] = 1
        else:
            digits[position_2] = 2
        count += 1

    accumulator += count

print(accumulator/sample_size)
