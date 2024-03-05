# Differencial forms 

## Why $1/n!$ outside of "components"
Differencial forms are usually written as $da = \frac{1}{n!} a_{\mu_1 \cdots \mu_n}  dx^{\mu_1} \wedge \cdots \wedge dx^{\mu_n}$. This notation is mysterious because the $n!$ factor seems to be unjustified. To see why it is important in defining components $a_{\mu_1 \cdots \mu_n}$, we need to rewrite the sum (initially over unordered multiplets $\mu_1 \cdots \mu_n$) as a sum over ordered multiplets denoted $\text{ord}(\mu_1 \cdots \mu_n)$ essentially a sum on multiplets where $\mu_i$ are arranged in increasing order. We have 
$$
da = \frac{1}{n!} \sum_{\mu_1 \cdots \mu_n} a_{\mu_1 \cdots \mu_n}  dx^{\mu_1} \wedge \cdots \wedge dx^{\mu_1} = \frac{1}{n!} \sum_{\text{ord}(\mu_1 \cdots \mu_n)} \sum_{p \in P_n} a_{p(\text{ord}(\mu_1 \cdots \mu_n))}  dx^{p(\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_n))}
$$
where $P_n$ is the group of permutation on $n$ objects, $p$ is a specific permutation and $\sigma(p)$ is the associated signature. Using the anticommutation of the wedge product, we can interchange the $1$-forms to place them in increasing $\mu_i$ order in each term to get (einstein summation is not used when a sum is written explicitly)
$$
da = \frac{1}{n!} \sum_{\text{ord}(\mu_1 \cdots \mu_n)} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_n)} \sum_{p \in P_n} a_{p(\text{ord}(\mu_1 \cdots \mu_n))} (-1)^{\sigma(p)} = \sum_{\text{ord}(\mu_1 \cdots \mu_n)} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_n)} a_{[\mu_1 \cdots \mu_n]}
$$
so the factorial fraction removes the redundance in the summation due to the multiple representation of the same tensor $dx^1 \wedge \cdots \wedge dx^n$ with wedge products. 
## Why $(l + k)!/(l! k!)$ inside of the components of a wedge product
When multiplying two forms $da = \frac{1}{l!} a_{\mu_1 \cdots \mu_l}  dx^{\mu_1} \wedge \cdots \wedge dx^{\mu_l}$ and  $db = \frac{1}{k!} b_{\nu_1 \cdots \nu_k}  dx^{\nu_1} \wedge \cdots \wedge dx^{\nu_k}$ with a wedge product we get a new form 
$$
dc = \frac{1}{l!} a_{\mu_1 \cdots \mu_l} dx^1 \wedge \cdots \wedge dx^l \wedge \frac{1}{k!} b_{\nu_1 \cdots \nu_k}  dx^{\nu_1} \wedge \cdots \wedge dx^{\nu_k} = \frac{1}{(l + k)!} c_{\mu_1 \cdots \mu_{l + k}}  dx^{\mu_1} \wedge \cdots \wedge dx^{\mu_{l+k}}. 
$$
Using the notation established before, we have 
$$
dc 
= \frac{1}{l!k!} a_{\mu_1 \cdots \mu_l} b_{\nu_1 \cdots \nu_k} dx^{\mu_1} \wedge \cdots \wedge dx^{\mu_l} \wedge dx^{\nu_1} \wedge \cdots \wedge dx^{\nu_k}\\
= \sum_{\text{Ord}(\mu_1 \cdots \mu_l)}\sum_{\text{ord}(\nu_1 \cdots \nu_k)} a_{\text{Ord}(\mu_1 \cdots \mu_l)} b_{\text{Ord}(\nu_1 \cdots \nu_k)} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_l)} \wedge dx^{\text{Ord}(\nu_1} \wedge \cdots \wedge dx^{\nu_k)}\\
= \sum_{\text{Ord}(\mu_1 \cdots \mu_l)}\sum_{\text{ord}(\nu_1 \cdots \nu_k)} a_{\text{Ord}(\mu_1 \cdots \mu_l)} b_{\text{Ord}(\nu_1 \cdots \nu_k)} (-1)^{\sigma(\text{Ord}, \text{Ord} \to \text{Ord})} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_l} \wedge dx^{\nu_1} \wedge \cdots \wedge dx^{\nu_k)}\\
= \sum_{\text{Ord}(\mu_1 \cdots \mu_{k+l})} \sum_{\text{Split}(l, k)} a_{\text{Split}_l} b_{\text{Split}_k} (-1)^{\sigma(\text{Split}_l, \text{Split}_k \to \text{Ord})} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_{l+k})}.
$$
Where $\text{Split}(l, k)$ is a pair of ordered multiplets of length $l$ ($\text{Split}_l$) and $k$ ($\text{Split}_k$) of values chosen in the $k+l$ multiplet $\text{Ord}(\mu_1 \cdots \mu_{k+l})$. The concatenation $\text{Split}_l \text{Split}_k$ is not ordered in general and requires transpositions from $\text{Ord}(\mu_1 \cdots \mu_{k+l})$. The signature of this numbur of transposition is $\sigma(\text{Split}_l, \text{Split}_k \to \text{Ord}) = \sigma(\text{Ord} \to \text{Split}_l, \text{Split}_k)$. We note that any permutation is the combination of a spliting and internal permutations in each block. Expanding again the internal permutation giving $a,b$ antisymmetric structure, we get 
$$
dc = \sum_{\text{Ord}(\mu_1 \cdots \mu_{k+l})} \sum_{\text{Split}(l, k)} \frac{1}{l!k!}\sum_{p \in P_k, q \in P_l} a_{p(\text{Ord}(\mu_1 \cdots \mu_l))} b_{q(\text{Ord}(\mu_{l+1} \cdots \mu_{l+k}))} (-1)^{\sigma(p) + \sigma(q) +\sigma(\text{Split}_l, \text{Split}_k \to \text{Ord})} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_{l+k})}\\
= \sum_{\text{Ord}(\mu_1 \cdots \mu_{k+l})} \frac{1}{l!k!}\sum_{p \in P_{k+l}} a_{p(\text{Ord}(\mu_1 \cdots \mu_l} b_{\mu_{l+1} \cdots \mu_{l+k}))} (-1)^{\sigma(p)} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_{l+k})}\\
= \sum_{\text{Ord}(\mu_1 \cdots \mu_{k+l})} \frac{(l+k)!}{l!k!} a_{[\mu_1 \cdots \mu_l} b_{\mu_{l+1} \cdots \mu_{l+k}]} dx^{\text{Ord}(\mu_1} \wedge \cdots \wedge dx^{\mu_{l+k})}
$$
So the new factor inside the components is there because of the disagreement of antisymmetrisation of individual multiplets of components of $a$ and $b$ against the antisymmetrisation of the antisymmetrisation of the combined multiplets. The combinatorial coefficient is the number of ways to choose $k$ from $l+k$ which is the number of splittings with parts with fixed order (increasing order). It is easier to see everything all at once with 
$$
(l+k)! a_{[\mu_1 \cdots \mu_l} b_{\mu_{l+1} \cdots \mu_{l+k}]} = l!k! \sum_{\text{Split}(l, k)} a_{\text{Split}_l} b_{\text{Split}_k} (-1)^{\sigma(\text{Split}_l, \text{Split}_k \to \text{Ord})} 
$$
