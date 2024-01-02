#%%
import numpy as np 


def c_wrap(i, spin): 
    # s : spin is 0 for down and 1 for up 
    # i : site where c acts from 0 to (2*s+1 - 1)
    # (i, s) gives the location of the target 
    # n : state on which c acts 
    # sign : coefficient of the state

    def c(n, sign):
    
        if sign == 0: return n, sign 
         # if the state is the null state, there is nothing to do

        n_bin = bin(n)[2:][::-1]
        # binary expression of n, reversed to keep operator action convention

        if 2 * i + spin < len(n_bin): filled = n_bin[2 * i + spin]
        else: filled = 0
        # if the target is out of binary expansion, return null
    
        if not int(filled): return n, 0
    
        else:
            # count the number of "1" before the target
            anti_c = n_bin[:2 * i + spin].count('1')
            # Removes an electron at the filled target
            new_n = n - 2 ** (2 * i + spin)
            return new_n, sign * (-1) ** anti_c
    
    def cdag(n, sign):
    
        if sign == 0: return n, sign 

        n_bin = bin(n)[2:][::-1]

        if 2 * i + spin < len(n_bin): filled = n_bin[2 * i + spin]
        else: filled = 0
    
        if int(filled): return n, 0
        else:
            anti_c = n_bin[:2 * i + spin].count('1')
            new_n = n + 2 ** (2 * i + spin)
            return new_n, sign * (-1) ** anti_c
        
    return c, cdag


testn = 3
N = 4

def format_state(testn, N):
    n = bin(testn)[2:][::-1]
    n = n + '0' * (2*N - len(n))
    spl = list(n)
    sites = []
    for i in range(N):
        sites.append(''.join(spl[2*i:2*i + 2]))
    return '__'.join(sites)

        

N = 4
testn = 4
print(format_state(testn, N))
print('-'*10)

for i in range(N):
    aup, cup = c_wrap(i, 1) 
    adown, cdown = c_wrap(i, 0)

    n_aup, s_aup = aup(testn, 1)
    n_cup, s_cup = cup(testn, 1)
    n_adown, s_adown = adown(testn, 1)
    n_cdown, s_cdown = cdown(testn, 1)


    n1 = format_state(n_aup, N)
    n2 = format_state(n_cup, N)
    n3 = format_state(n_adown, N)
    n4 = format_state(n_cdown, N)

    print(str(s_aup) + "*" + n1, 'a up {}'.format(i))
    print(str(s_cup) + "*" + n2, 'c up {}'.format(i))
    print(str(s_adown) + "*" + n3, 'a down {}'.format(i))
    print(str(s_cdown) + "*" + n4, 'c down {}'.format(i))
    print('-'*10)


# %%

# %%
