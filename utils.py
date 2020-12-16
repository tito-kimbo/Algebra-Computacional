class AssumptionnError(Exception):
    pass


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
