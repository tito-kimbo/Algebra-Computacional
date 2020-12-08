def shifting_root(n, e, b):
    """
        Finds the nth root of e in base b using the shifting algorithm
        e must be represented as a list of digits in base b
    """

    def sum_base(e,b):
        return sum([b**i * x for i,x in enumerate(e[::-1])], b.ring.zero)

    R = b.ring

    x = R.zero
    y = R.zero
    r = R.zero

    mod = len(e) % n
    if mod > 0:
        alpha = sum_base(e[0:mod], b)
        i = mod
    else:
        alpha = sum_base(e[0:n], b)
        i = n

    while i <= len(e):

        x = b**n * x + alpha
        
        beta = R.zero
        while (b*y + beta)**n < x:
            beta = beta + R.one

        y = b*y + beta
        r = x - y**n

        i += n
        alpha = sum_base(e[i-n:i], b)

    return y
