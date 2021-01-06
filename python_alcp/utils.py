from collections import defaultdict
from itertools import chain, combinations
from functools import reduce

######## Error checking #########

"""
    Alternative for asserts
    Usage:
        def myfun(myarg):
            assuming(myarg > 0, "MyArg must be >0")
"""
def assuming(val: bool, msg:str = None):
    if not val:
        if msg is not None:
            raise AssertionError(f"Unfulfilled assumption: {msg}")
        else:
            raise AssertionError(f"Unfulfilled assumption")


def try_op(op_name = None):
    def op_decorator(op):
        nonlocal op_name
        if op_name == None:
            op_name = op.__name__.strip("_")
        def wrapped_op(*args, **kwargs):
            try:
                return op(*args, **kwargs)
            except (AssertionError, TypeError) as e:
                raise ValueError(f"Can't {op_name} {' and '.join([str(a) for a in args])}: {e}")
        return wrapped_op

    return op_decorator



OPERAND_ERROR = "All operands must be within the following list "
def op_typecheck(operand,allowed):
    if not type(operand) in allowed:
        raise TypeError(OPERAND_ERROR+" "+str(allowed))


######## Pretty printing #########

def print_superscript(n: int):
    uni = [
        "\U00002070",
        "\U000000B9",
        "\U000000B2",
        "\U000000B3",
        "\U00002074",
        "\U00002075",
        "\U00002076",
        "\U00002077",
        "\U00002078",
        "\U00002079"
    ]
    digits = []
    while n > 0:
        digits.append(uni[n % 10])
        n //= 10

    return ''.join(digits[::-1])


######## Primes and int factorization #########

curprimes = [2,3]
def primes():

    def _p():
        # Generator that returns primes
        
        count = 0
        while count < len(curprimes):
            yield curprimes[count]
            count += 1

        i = curprimes[-1] + 2
        while True:
            prime = True
            for p in curprimes:
                if i % p == 0:
                    prime = False
                    break

            if prime:
                curprimes.append(i)
                yield i

            i += 2

    return _p()



def prime_factors(n: int):
    # Crude factorization
    res = defaultdict(lambda: 0)
    for p in primes():
        if n == 1:
            return res
        elif p**2 > n:
            res[n] += 1
            return res
        while n % p == 0:
            res[p] += 1
            n //= p

def all_factors(n: int):
    facts = sum([[f]*mul for f,mul in prime_factors(n).items()], [])
    def mul(nums):
        return reduce(int.__mul__, nums, 1)
    return set((mul(c) for r in range(1,len(facts)) for c in combinations(facts, r)))



######## Double-add algorithm #########

def double_and_add(element, n: int, op, neutral):
    """ Computes op^n(element) in O(log(n)) calls to op """

    res = neutral
    acc = element
    while n > 0:
        if n % 2 == 1:
            res = op(res,acc)
        n //= 2
        acc = op(acc,acc)

    return res


def repeated_addition(element, n: int):
    """ Computes a+a+a+a... (n times) in O(log(n)) time """

    return double_and_add(element, n, type(element).__add__, type(element).zero)

def fast_exponentiation(element, n: int):
    """ Computes a^n in O(log(n)) time """

    return double_and_add(element, n, type(element).__mul__, type(element).one)





######## Externals #########

"""
    External operations. Avoids circular imports
    Usage:
        If a module needs to use a function but importing it would
        cause circular imports, decorate the function with @external
        and then use it in the requesting module as externals.<func>

    Example:
        rings.Ring requires ideals.GetQuotient to implement its operator
        __truediv__, but GetQuotient is defined on ideals.py, which
        depends on rings.Ring.
        To solve this, ideals imports utils.external and decorates GetQuotient
        Then rings.Ring imports externals, and can call externals.GetQuotient
        without issues
        
"""

class Externals():
    pass

externals = Externals()

def external(elem):
    externals.__setattr__(elem.__name__, elem)
    return elem
